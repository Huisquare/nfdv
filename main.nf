/*----------------------------------------------------------------------
  DeepVariant as a Nextflow pipeline for whole genome sequencing data
----------------------------------------------------------------------*/


//  Fasta, indexed fasta and zipped input files

assert (params.fasta != true) && (params.fasta != null) : "please specify --fasta path/to/fasta_file"

fasta=file(params.fasta)

//  Bam and indexed bam input files

assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder path/to/bam_folder"

Channel.fromPath("${params.bam_folder}/*.bam").map{ file -> tuple(file.name, file) }.set{bamChannel}

//  output directory

params.resultdir = "results";


//  generate indexed files and zipped files from the input fasta file 
//  file types: .fai, .gz, .gz.fai, .gz.gzi

process preprocessFASTA{

  container 'huisquare/htslib-and-samtools'
  publishDir "$baseDir/intermediate/preprocessFASTA"


  input:
  file fasta from fasta

  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi") into fastaChannel
  script:
  """
  samtools faidx $fasta ;
  bgzip -c ${fasta} > ${fasta}.gz ;
  bgzip -c -i ${fasta} > ${fasta}.gz ;
  samtools faidx "${fasta}.gz" ;
  """

}

//  Params for the Read Group Line to be added in case it is needed.

params.rgid=4;
params.rglb="lib1";
params.rgpl="illumina";
params.rgpu="unit1";
params.rgsm=20;


//  Produces indexed bam files if user did not provide them

process preprocessBAM{

  tag "${bam[0]}"
  container 'huisquare/samtools-picard-config'
  publishDir "$baseDir/intermediate/preprocessBAM"

  input:
  set val(prefix), file(bam) from bamChannel
  output:
  set file("final/${bam[0]}"), file("final/${bam[0]}.bai") into completeChannel
  script:
  """
  mkdir final
  [[ `samtools view -H ${bam[0]} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam final;}|| \
  { picard AddOrReplaceReadGroups \
    I=${bam[0]} \
    O=final/${bam[0]} \
    RGID=${params.rgid} \
    RGLB=${params.rglb} \
    RGPL=${params.rgpl} \
    RGPU=${params.rgpu} \
    RGSM=${params.rgsm};}
    cd final;
    samtools index ${bam[0]};
  """
}

fastaChannel.map{file -> tuple (1,file[0],file[1],file[2],file[3],file[4])}
            .set{all_fa};

completeChannel.map { file -> tuple(1,file[0],file[1]) }
               .set{all_bam};

// below code will create 
// [[1, ref.fa, ref.fa.fai, ref.fa.gz, ref.fa.gz.gzi, ref.gz.fai], [1, bam_file_0.bam, bam_file_0.bai]]
// [[1, ref.fa, ref.fa.fai, ref.fa.gz, ref.fa.gz.gzi, ref.gz.fai], [1, bam_file_1.bam, bam_file_1.bai]]

all_fa.cross(all_bam)
      .set{all_fa_bam};

//  Number of cores to be used for makeExamples

int cores = Runtime.getRuntime().availableProcessors();
params.numCores=cores
numCoresMinusOne=params.numCores-1;

//  Getting bam files and converting them to images (named examples)

process makeExamples{

    tag "${bam[1]}"
    cpus params.numCores

    input:
      set file(fasta), file(bam) from all_fa_bam
    output:
      set file("${fasta[1]}"),file("${fasta[1]}.fai"),file("${fasta[1]}.gz"),file("${fasta[1]}.gz.fai"), file("${fasta[1]}.gz.gzi"),val("${bam[1]}"), file("shardedExamples") into examples
    shell:
    '''
    mkdir shardedExamples
    time seq 0 !{numCoresMinusOne} | \
    parallel --eta --halt 2 \
      python /opt/deepvariant/bin/make_examples.zip \
      --mode calling \
      --ref !{fasta[1]}.gz\
      --reads !{bam[1]} \
      --examples shardedExamples/examples.tfrecord@!{params.numCores}.gz\
      --task {}
    '''
}

//  Doing the variant calling based on the ML trained model.

process call_variants{

  tag "${bam}"
  cpus params.numCores

  input:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"),val(bam), file("shardedExamples") from examples
  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"), val(bam), file('call_variants_output.tfrecord') into called_variants
  script:
  """
  /opt/deepvariant/bin/call_variants \
    --outfile call_variants_output.tfrecord \
    --examples shardedExamples/examples.tfrecord@${params.numCores}.gz \
    --checkpoint /opt/models/wgs/model.ckpt \
    --num_readers ${params.numCores}
  """
}


//  Transforming the variant calling output (tfrecord file) into a standard vcf file.

process postprocess_variants{

  tag "${bam}"
  cpus params.numCores

  publishDir params.resultdir, mode: 'copy'
  input:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"), val(bam),file('call_variants_output.tfrecord') from called_variants
  output:
   set val(bam),file("${bam}.vcf") into postout
  script:
  """
    /opt/deepvariant/bin/postprocess_variants \
    --ref "${fasta}.gz" \
    --infile call_variants_output.tfrecord \
    --outfile "${bam}.vcf"
  """
}

//  completion of workflow

workflow.onComplete {
  if (workflow.success){
    System.out.println("Congrats! The job was successful. Please find your results in $baseDir/${params.resultdir}");
  } else {
    System.out.println("Job was not successful. Please look at the .nextflow.log file for more information on errors.");
  }
}
