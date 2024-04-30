// Variant calling pipeline, from FASTQ to BAM to gVCF
// (no gVCF merging)


// ########### SOFTWARE AND CONTAINERS ############
// gatk_docker = "public.ecr.aws/biocontainers/gatk4:4.1.9.0--py39_0"
gatk_invoc = params.gatk_invoc          // "gatk"
// picard_docker = gatk_docker // since gatk v4
picard_invoc = params.picard_invoc     // "java -jar /usr/gitc/picard.jar "
bwa_invoc = params.bwa_invoc           //    "/usr/gitc/bwa "
samtools_invoc = params.samtools_invoc // "samtools"

// ############################################################### //
// ##############   Prepare inputs ############################### //
// ############################################################### //
// Parsing sample sheet, sending to channel "samples_ch"
Channel.fromPath(params.samplesheet)
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

 // container "${picard_docker}" // container settings in config file

 input:
 tuple val(sid), val(rg), val(num_rg),  path(bam)

 output:
 tuple val(sid), val(rg), val(num_rg),  path("${sid}.${rg}.proc.bam"), path("${sid}.${rg}.proc.bai"),  emit: proc_bam
 path "*.metrics",  emit: mkdup_metrics_file, optional: true
 //path("${sid}*.bai"), optional: true

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
 // storeDir 'store-gvcf'
 publishDir "${params.outdir}/bamout/", pattern: "*bamOut.ba*", mode: 'copy'

 containerOptions " -v ${params.tmpdir}:/tmpdir "

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
 ## for older GATK version
 ##optionalT=" -T "

 bands="--gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30  --gvcf-gq-bands 33 --gvcf-gq-bands 35 --gvcf-gq-bands 37  --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 42  --gvcf-gq-bands 45 --gvcf-gq-bands 46  --gvcf-gq-bands 48  --gvcf-gq-bands 50 --gvcf-gq-bands 52  --gvcf-gq-bands 55  --gvcf-gq-bands 57  --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99"

 # consider adding --gvcf-bands


 echo "#----------- Working directory ------------#"
 pwd
 echo "#------------------------------------------#"

 $gatk_invoc --java-options "-Xmx${avail_mem}M -XX:-UsePerfData"  HaplotypeCaller \
         -R $ref \
         -ERC  GVCF \
         -I $bam \
         -O ${sid}.g.vcf.gz  \
         -mbq 20 \
          --showHidden true \
         -bamout ${sid}.bamOut.bam \
       --tmp-dir tmpdir

  #
  #      \$bands

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
   // variantCall( ref_files_for_varcall, bamProcess.out.proc_bam )

 emit:
   variantCall.out.gvcf
}
