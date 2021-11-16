# Incorporating Google's DeepVariant variant caller as a Nextflow pipeline on Google Cloud Platform

A ZB4171 Project, incorporating Google DeepVariant with Nextflow pipeline. 
## Our Project 
1. Incorporating DeepVariant as a Nextflow pipeline

2. Hosting the pipeline on Google Cloud Platform (GCP)

3. Perform pipeline on a breast cancer dataset against clinical studies, and cross validation
    a. Identify and analyse gene variants from the result
    b. Supporting evidences for previously identified breast cancer
    c. Discover possible new genes associated with breast cancer

## Introduction 
### “Google’s DeepVariant: Highly Accurate Genomes With Deep Neural Networks” 
DeepVariant transforms the task of variant calling, a reconstruction problem in genomics, into an image classification problem. It reconstructs the true genome sequence from high throughput sequencer data with significantly greater accuracy than previous classical methods

### “Nextflow: Workflow framework that eases data-intensive computational pipelines” 
Benefits of Nextflow: 
  Fast prototyping of computational pipeline
  Support containerisation (e.g. Kubernetes, Docker)
  Portable
  Parallelisation
  Continuous checkpointing
  All the intermediate results produced during the pipeline execution are automatically tracked

### Google Cloud Platform (GCP) 
Firstly, new customers are offered $300 free credits to run, test and deploy programs on Google Cloud. As DeepVariant is the brainchild of Google Life Sciences, its workflows can be run easily on Google Cloud Platform with optimised configurations. There is also support for running Nextflow pipelines with Google Life Sciences API.

## Motivation
DeepVariant specifies that it only takes in bam or cram files and its indexed bam file for aligned reads. It requires fasta files and its index file for reference genomes. There is a need to provide input files in such format to run DeepVariant. In addition, when multiple bam files are passed into DeepVariant, it runs them one by one. Also, DeepVariant requires high computing power, therefore running it on local machine may not be an ideal choice.

### Preprocessing & Parallelisation 
The nextflow pipeline aims to automatically handle the creation of some extra needed index files such as the fai and the bai files, which are needed as inputs for DeepVariant. This files were normally produced by users. The support for Docker in Nextflow allows us to bundle each preprocessing step with an individual Docker container. In addition, incorporating DeepVariant into a Nextflow pipeline allows DeepVariant to leverage on the parallelisation that Nextflow offers. This allows variant calling to be performed at the same time on multiple bam files. Finally, the nextflow pipeline can be deployed onto cloud to use the high computing power offered by cloud computing companies.

## Dependencies

[Nextflow](https://www.nextflow.io/)
[Docker](https://www.docker.com/)

## Quick Start

A typical run on **whole genome data** looks like this: 
```
git clone https://github.com/Huisquare/nfdv
cd nfdv
nextflow run main.nf --fasta path/to/fastaFile --bam_folder path/to/bamFolder
```
The h38 version of the reference genome is used.
Two vcf files are produced and can be found in the folder "results" in the home directory

## More about the pipeline 

The workflow **accepts one reference genome and multiple BAM files as input**. The variant calling for the several input BAM files will be processed completely indipendently and will produce indipendent VCF result files. The advantage of this approach is that the variant calling of the different BAM files can be parallelized internally by Nextflow and take advantage of all the cores of the machine in order to get the results at the fastest.


## INPUT PARAMETERS

### About preprocessing

DeepVariant, in order to run at its fastest, requires some indexed and compressed versions of both the reference genome and the BAM files. With DeepVariant in Nextflow, if you wish, you can only use as an input the fasta and the BAM file and let us do the work for you in a clean and standarized way (standard tools like [samtools](http://samtools.sourceforge.net/) are used for indexing and every step is run inside of  a Docker container).

<!-- This is how the list of the needed input files looks like. If these are passed all as input parameters, the preprocessing steps will be skipped. 
```
NA12878_S1.chr20.10_10p1mb.bam   NA12878_S1.chr20.10_10p1mb.bam.bai	
ucsc.hg19.chr20.unittest.fasta   ucsc.hg19.chr20.unittest.fasta.fai 
ucsc.hg19.chr20.unittest.fasta.gz  ucsc.hg19.chr20.unittest.fasta.gz.fai   ucsc.hg19.chr20.unittest.fasta.gz.gzi
```
If you do not have all of them, these are the file you can give as input to the Nextflow pipeline, and the rest will be automatically  produced for you .
```
NA12878_S1.chr20.10_10p1b.bam  
ucsc.hg19.chr20.unittest.fasta
``` -->

### Parameters definition 

- ### BAM FILES 

```
--bam_folder "/path/to/folder/where/bam/files/are"            REQUIRED
--getBai "true"                                               OPTIONAL  (default: "false")
```
In case only some specific files inside the BAM folder should be used as input, a file prefix can be defined by: 
```
--bam_file_prefix MYPREFIX
```

All the BAM files on which the variant calling should be performed should be all stored in the same folder. If you already have the index files (BAI) they should be stored in the same folder and called with the same prefix as the correspoding BAM file ( e.g. file.bam and file.bam.bai ). 

**! TIP** 

- ### REFERENCE GENOME

 By default the h38 version of the reference genome is used. 
 
 Alternatively, a user can use an own reference genome version, by using the following parameters:

  ```
  --fasta "/path/to/myGenome.fa"                OPTIONAL
  --fai   "/path/to/myGenome.fa.fai"            OPTIONAL
  --fastagz "/path/to/myGenome.fa.gz"           OPTIONAL
  --gzfai  "/path/to/myGenome.fa.gz.fai"         OPTIONAL
  --gzi  "/path/to/myGenome.fa"                  OPTIONAL
  ```
If the optional parameters are not passed, they will be automatically be produced for you and you will be able to find them in the "preprocessingOUTPUT" folder.

### Advanced parameters options

- ### CPUS 

The **make_example** process can be internally parallelized and it can be defined how many cpus should be assigned to this process.
By default all the cpus of the machine are used.

```
-- numCores 2          OPTIONAL (default: all)
```
- ### MODEL 

The trained model which is used by the **call_variants** process can be changed.
The default one is the 0.6.0 Version for the whole genome. So if that is what you want to use too, nothing needs to be changed.

## More about the dataset
We have chosen to look at the HCC1143 cell line, which is a publicly available illumina whole genome sequencing data. The cell line was generated from a 52 year old caucasian woman with breast cancer tumor. Fastq files of both matched normal and tumor were preprocessed, subjected to GATK best practices. The bam files containing the reads for the cancer cell line and the matched normal consists of these 2 bam files. We chose to only look at reads from chromosome 17 as we wanted to start with a smaller dataset to test our pipeline. The genome sequence reads were aligned to the Human GRCh38 reference genome.

## More about Docker Containers 
To ensure that such tools used in the Nextflow pipeline can run on any machine without running into errors of uninstalled dependencies, Docker containers are used. We are able to encapsulate each process in a Docker container. The configurations needed in the container can be specified in a Docker image, and the image is used like a piece of instruction to build the Docker container. When a process is being run, it would be running in a container, which we can think of as an environment that is specially configured for that process. Therefore, there are no worries about the environment in the deployed machine affecting the execution of each process.
The Docker files that we wrote for each process can be found in the **Dockerfiles** folder. The images of the Docker files can be found on Dockerhub. 
[htslib-and-samtools](https://hub.docker.com/repository/docker/huisquare/htslib-and-samtools)
[samtools-config](https://hub.docker.com/repository/docker/huisquare/samtools-config)
[vcftools-config](https://hub.docker.com/repository/docker/huisquare/vcftools-config)


### Acknowledgements
We referenced similar pipelines developed by [lifebit.ai] (https://github.com/lifebit-ai/DeepVariant) and [nf-core] (https://github.com/nf-core/deepvariant) when building our pipeline.

We were able to run our pipeline on Google Cloud due to the USD$300 free credits that we received as new users in Google Cloud. 

### Improvements from current similar pipelines
The DeepVariant model we used is the 1.2.0 version, which is a huge advancement from v0.5.1/v0.6.1 used by [lifebit.ai] (https://github.com/lifebit-ai/DeepVariant) and v1.0.0 used by [nf-core] (https://github.com/nf-core/deepvariant).

The model used is within the DeepVariant docker container, instead of using a model stored on cloud. Our pipeline is therefore more efficient in this aspect as it does not need to download the trained DeepVariant model from cloud storage.
