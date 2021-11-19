/*----------------------------------------------------------------------
  DeepVariant as a Nextflow pipeline for whole genome sequencing data
----------------------------------------------------------------------*/


//  Number of cores to be used for makeExamples

int cores = Runtime.getRuntime().availableProcessors();
params.numCores=cores
numberShardsMinusOne=params.numCores-1;


//  Fasta, indexed fasta and zipped input files

params.test="";

params.fasta="nofasta";
params.fai="nofai";
params.fastagz="nofastagz";
params.gzfai="nogzfai";
params.gzi="nogzi";

if(!("nofasta").equals(params.fasta)){
  fasta=file(params.fasta)
  fai=file(params.fai);
  fastagz=file(params.fastagz);
  gzfai=file(params.gzfai);
  gzi=file(params.gzi);
}

else if(params.test){
  fasta=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta");
  fai=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.fai");
  fastagz=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.gz");
  gzfai=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.gz.fai");
  gzi=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.gz.gzi");
}

else{
  System.out.println("please input your fasta file using --fasta \"/path/to/your/genome\" ");
  System.exit(1);
}


//  Bam and indexed bam input files

params.getBai="false";

if(params.test){
    params.bam_folder="$baseDir/testdata"
}

assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder path/to/bam_folder"

params.bam_file_prefix="*"

if( !("false").equals(params.getBai)){
  Channel.fromFilePairs("${params.bam_folder}/${params.bam_file_prefix}*.{bam,bam.bai}").set{bamChannel}
}else{
  Channel.fromPath("${params.bam_folder}/${params.bam_file_prefix}*.bam").map{ file -> tuple(file.name, file) }.set{bamChannel}
}

//  output directory

params.resultdir = "results";


//  generate indexed files and zipped files from the input fasta file if user did not provide them in input
//  file types: .fai, .gz, .gz.fai, .gz.gzi

process preprocessFASTA{

  container 'huisquare/htslib-and-samtools'
  publishDir "$baseDir/sampleDerivatives"


  input:
  file fasta from fasta
  file fai from fai
  file fastagz from fastagz
  file gzfai from gzfai
  file gzi from gzi
  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi") into fastaChannel
  script:
  """
  [[ "${params.fai}"=="nofai" ]] &&  samtools faidx $fasta || echo " fai file of user is used, not created"
  [[ "${params.fastagz}"=="nofastagz" ]]  && bgzip -c ${fasta} > ${fasta}.gz || echo "fasta.gz file of user is used, not created "
  [[ "${params.gzi}"=="nogzi" ]] && bgzip -c -i ${fasta} > ${fasta}.gz || echo "gzi file of user is used, not created"
  [[ "${params.gzfai}"=="nogzfai" ]] && samtools faidx "${fasta}.gz" || echo "gz.fai file of user is used, not created"
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
  container 'huisquare/samtools-config'
  publishDir "$baseDir/sampleDerivatives"

  input:
  set val(prefix), file(bam) from bamChannel
  output:
  set file("ready/${bam[0]}"), file("ready/${bam[0]}.bai") into completeChannel, completeStats
  script:
  """
	mkdir ready
  [[ `samtools view -H ${bam[0]} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam ready;}|| { picard AddOrReplaceReadGroups \
    I=${bam[0]} \
    O=ready/${bam[0]} \
    RGID=${params.rgid} \
    RGLB=${params.rglb} \
    RGPL=${params.rgpl} \
    RGPU=${params.rgpu} \
    RGSM=${params.rgsm};}
    cd ready ;samtools index ${bam[0]};
  """
}


//  Use samtools to collect statistical information of the alignments

process BAMstats{

  tag "${bam[0]}"
  container 'huisquare/samtools-config'

  input:
  set file(bam), file(bai) from completeStats
  output:
  file("*") into bam_multiqc
  script:
  """
  samtools stats $bam > stats.txt
  samtools flagstat $bam > flagstat.txt
  samtools idxstats $bam > idxstats.txt
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
    time seq 0 !{numberShardsMinusOne} | \
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

  tag "$bam"
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

//  Use vcftools to collect data on transitions and transversions.

process vcftools{
  tag "$vcf"

  container 'huisquare/vcftools-config'

  input:
  set val(bam),file(vcf) from postout
  output:
  file("*") into vcfout

  script:
  """
  vcftools --vcf $vcf --TsTv-summary
  vcftools --vcf $vcf --TsTv-by-count
  vcftools --vcf $vcf --TsTv-by-qual
  # remove rows containing 'inf' which breaks multiqc report
  sed -i '/inf/d' out.TsTv.qual
  """
}

//  Use multiqc to generate a summary report.

process multiqc{
  tag "multiqc_report.html"

  publishDir "${params.resultdir}/MultiQC", mode: 'copy'
  container 'ewels/multiqc:latest'

  input:
  file(vcfout) from vcfout
  file(bamout) from bam_multiqc
  output:
  file("*") into multiqc

  script:
  """
  multiqc . -m vcftools -m samtools
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
