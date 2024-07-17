// Variant calling pipeline, from FASTQ to BAM to gVCF
// (no gVCF merging)


// ########### SOFTWARE AND CONTAINERS ############
// gatk_docker = "public.ecr.aws/biocontainers/gatk4:4.1.9.0--py39_0"
gatk_invoc = params.gatk_invoc          // "gatk"
picard_invoc = params.picard_invoc     // "java -jar /usr/gitc/picard.jar "
bwa_invoc = params.bwa_invoc           //    "/usr/gitc/bwa "
samtools_invoc = params.samtools_invoc // "samtools"

// ############################################################### //
// ##############   Prepare inputs ############################### //
// ############################################################### //
// Parsing sample sheet, sending to channel "samples_ch"
Channel.fromPath(params.input)
  .splitCsv(header: true)
  .map { row -> tuple( row.sample
                       ,row.readgroup
                       ,file(row.fastq_1)
                       ,file(row.fastq_2)
                       ,row.num_rg.toInteger()
                       ) }
  .dump(tag: "samples_ch")
  .set { samples_ch }

// channels for reference
ref_rel_files = Channel.value( params.ref )
            .map { f -> tuple(
                 file(f),
                 file("${f}.amb"),
                 file("${f}.ann"),
                 file("${f}.bwt"),
                 file("${f}.pac"),
                 file("${f}.sa")
                 )  }

ref_files_for_varcall = Channel.value( tuple(
               file(params.ref),
               file(params.ref_fai),
               file(params.ref_dict)
               ) )


// ####################################################### //
//                 PROCESS - ALIGN
// ####################################################### //
process ALIGN {
  container 'public.ecr.aws/p1p2i7f9/samtools-picard-bwa:v1'

 input:
 tuple path(ref), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa)
 tuple val(sid), val(rg), path(read1), path(read2), val(num_rg)

 output:
 tuple val(sid), val(rg), val(num_rg),  path("${sid}*${rg}.bam"), emit: bam
 stdout  emit: verbiage

 script:
 """
 echo "cpu: ${task.cpus}  mem: ${task.memory}"
 echo "${workflow.containerEngine}"
 echo "##   "  ${sid} ':' ${read1} ${read2}
 $bwa_invoc mem -t ${task.cpus}  ${ref} ${read1} ${read2} > out.sam
 samtools view -Obam,level=5 out.sam > ${sid}_${rg}.bam
 rm *.sam # if conversion successful, remove SAM
 """
}


// ####################################################### //
// ##### Processing BAMs
// ####################################################### //

process bamProcess {
 container 'public.ecr.aws/p1p2i7f9/samtools-picard-bwa:v1'

 publishDir 'stats', pattern: "*.metrics", mode: 'move'
 storeDir 'processed-bams'  // mode: 'copy'

 input:
 tuple val(sid), val(rg), val(num_rg),  path(bam)

 output:
 tuple val(sid), val(rg), val(num_rg),  path("${sid}.${rg}.proc.bam"), path("${sid}.${rg}.proc.bai"),  emit: proc_bam
 path "*.metrics",  emit: mkdup_metrics_file, optional: true

 script:
 """
 set -euo pipefail

 echo $bam
  # sort sam
    $picard_invoc  SortSam -I $bam -O sorted.bam \
    --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE \
    --SORT_ORDER coordinate \
    ${params.sort_extra_params} --TMP_DIR ${params.tmpdir}

  # fix mate
  $picard_invoc  FixMateInformation -I  sorted.bam  -O fxmt.bam \
      --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE \
    ${params.fixmate_extra_params} --TMP_DIR ${params.tmpdir}

  rm sorted.bam

  # mark dup
  $picard_invoc  MarkDuplicates -I  fxmt.bam  -O markdup.bam \
      --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE \
    --METRICS_FILE ${sid}.markdup.metrics \
    --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
    --TMP_DIR ${params.tmpdir} \
    ${params.markdup_extra_params}

  rm fxmt.bam

  # add/replace read groups
  $picard_invoc  AddOrReplaceReadGroups -I  markdup.bam  -O addrep.bam \
      --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE \
    --RGLB lib1 \
    --RGSM ${sid} \
    --RGID ${sid} \
    --RGPU PU1 \
    --RGPL ${params.PL} \
    --RGCN ${params.CN} \
	--TMP_DIR ${params.tmpdir} \
    ${params.addrep_extra_params}

  rm markdup.bam

  # set the name of final file
  mv addrep.bam ${sid}.${rg}.proc.bam
  mv addrep.bai ${sid}.${rg}.proc.bai
 """

}

// ################################################################ //
// ##### Merge BAM files from different read groups of the same sample
// ################################################################ //

process mergeBam {
 container 'public.ecr.aws/p1p2i7f9/samtools-picard-bwa:v1'

 input:
 tuple val(sample), path(input_files, stageAs: "?/*"), path(input_bai, stageAs: "?/*")

 output:
 tuple val(sample), path("${prefix}.bam"), path("${prefix}.bai"),  emit: bam


 script:
 prefix = task.ext.prefix ?: "${sample}"
 def args = task.ext.args   ?: ''
 def files_str = input_files instanceof List ? input_files.join(" ") : input_files
 """
 samtools merge  --threads ${task.cpus-1} \
   ${prefix}.bam $files_str

 samtools index ${prefix}.bam ${prefix}.bai
 """

}


// ######################################################### //
// ######        Variant Calling  ##
// ######################################################## //

process variantCall {
 container 'public.ecr.aws/biocontainers/gatk4:4.1.9.0--py39_0'

 publishDir "${params.outdir}/gvcf/", pattern: "*vcf*", mode: 'copy'
 publishDir "${params.outdir}/bamout/", pattern: "*bamOut.ba*", mode: 'copy'

 input:
 tuple path(ref), path(ref_fai), path(ref_dict)
 tuple val(sid), path(bam), path(bai)

 output:
 tuple val(sid), path("*.vcf.gz"), emit: gvcf
 path("*.bamOut.ba*"),  emit: bamOut

 script:

  def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

 """
 mkdir ./TMP

 $gatk_invoc --java-options "-Xmx${avail_mem}M -XX:-UsePerfData"  HaplotypeCaller \
         -R $ref \
         -ERC  GVCF \
         -I $bam \
         -O ${sid}.g.vcf.gz  \
         -mbq 20 \
          --showHidden true \
         -bamout ${sid}.bamOut.bam \
       --tmp-dir ./TMP

 """

}


// MAIN WORKFLOW
workflow {

  main:
     ALIGN( ref_rel_files, samples_ch)

    bamProcess( ALIGN.out.bam )

    bam_mapped = bamProcess.out.proc_bam . map{
           sample, rg, num_rg, bam, bai ->
            [ groupKey( sample, num_rg )
               , bam
               , bai
               ]
     } . groupTuple() .dump(tag: "mapped" )

     // Merge BAMs
    bam_to_merge = bam_mapped.branch{ sample, bam, bai
         ->
        // bam is a list, so use bam.size()
        single:   bam.size() <= 1
            return [  sample
                     , bam[0]
                     , bai[0]
            ]
        multiple: bam.size() > 1
    }

    // Apply merging on the multiple read-group samples
    mergeBam( bam_to_merge.multiple )

    // Add the single read group BAM files, creating new channel bam_all
    bam_all = mergeBam.out.bam.mix(bam_to_merge.single)
    bam_all.dump(tag: "bam_all")


     variantCall( ref_files_for_varcall, bam_all)

 emit:
   variantCall.out.gvcf
}
