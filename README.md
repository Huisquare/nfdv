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
git clone https://github.com/huisquare/nfdv
cd nfdv
nextflow run main.nf --fasta path/to/fastaFile --bam_folder path/to/bamFolder
```
The h38 version of the reference genome is used.
Two vcf files are produced and can be found in the folder "results" in the home directory

## More about Google Cloud Platform 

Note: Instructions are adapted from GCP. 

DeepVariant doesn't require GCP, but if you want to use it, these are some instructions that we found to be useful when getting started.

### Set up a Google Cloud account

To get started using DeepVariant on Google Cloud Platform (GCP), you first need
to set up an account and a project to contain your cloud resources.

*   If you do not have an account yet, you should create one at
    [cloud.google.com](https://cloud.google.com). You should then [enable
    billing for your
    account](https://support.google.com/cloud/answer/6288653?hl=en) but note
    that if your account is new, [you receive $300 of free
    credit](https://cloud.google.com/free/). Once your cloud account is set up,
    you should be able to log in to the [Cloud
    Console](https://console.cloud.google.com) to view or administer your cloud
    resources.

*   From the Cloud Console, [set up a
    project](https://cloud.google.com/resource-manager/docs/creating-managing-projects)
    to house all of the cloud resources (storage, compute, services) that you
    will associate with your use of DeepVariant. For example, if your
    organization is AcmeCorp, you might call your project
    `acmecorp-deepvariant`.

*   Finally, please visit the ["Compute Engine" page on Cloud
    Console](https://console.cloud.google.com/compute). You don't need to create
    Compute Engine instances at this time, but simply visiting this page will
    initialize your compute engine "service account" so that we can authorize
    it.

### Install the Google Cloud SDK

The Google Cloud SDK comes with two very useful command line utilities that you
can use on your local workstation---`gcloud`, which lets you administer your
cloud resources, and `gsutil`, which lets you manage and transfer data to Google
Cloud Storage buckets. We will make use of these tools in the following
instructions. To install the Cloud SDK, [follow the installation instructions
here](https://cloud.google.com/sdk/downloads).

The final step in the installation process (`gcloud init`) will have you
authenticate via your web browser and select a default [zone and
region](https://cloud.google.com/compute/docs/regions-zones/regions-zones) for
your cloud resources, which you can choose based on your location and regional
hardware availability.

NOTE: Not all zones are equipped with GPUs, so if you want to use GPUs for your
project, please take note of the availability listing
[here](https://cloud.google.com/compute/docs/gpus/).

To verify that the installation and authentication succeeded, run

```shell
gcloud auth list
```

and verify that your account email address is printed.

### Starting a Compute Engine instance

A simple way to access compute on GCP is Google Compute Engine. Compute Engine
instances can be sized to meet computational and storage needs for your project.

Before we get started, [ensure you have adequate quota
provisioned](https://cloud.google.com/compute/quotas) so that you can get all
the CPUs/GPUs that you need. To start with, you might want to request quota for
64 CPUs and 2 GPUs in your zone.

DeepVariant can make use of multiple CPU cores and (currently, a single) GPU
device. For this "quick start" guide, let's allocate an 8-core non-preemptible
instance in your default zone with a single GPU, running Ubuntu 20.04, with a
disk of reasonable size for modest work with genomic data. From our local
command line, we do:

```shell
gcloud beta compute instances create "${USER}-deepvariant-quickstart" \
  --scopes "compute-rw,storage-full,cloud-platform"  \
  --image-family ubuntu-2004-lts --image-project ubuntu-os-cloud \
  --machine-type n1-standard-8  \
  --boot-disk-size=200GB \
  --zone us-west1-b \
  --accelerator type=nvidia-tesla-k80,count=1 --maintenance-policy TERMINATE --restart-on-failure
```

NOTE: To create an instance *without GPU*, simply omit the last line from the
command.

Check that the instance has been created and started:

```shell
gcloud compute instances list
```

which should produce output like:

```
NAME                    ZONE        MACHINE_TYPE    PREEMPTIBLE   INTERNAL_IP  EXTERNAL_IP     STATUS
[USER]-deepvariant-quickstart  us-west1-b  n1-standard-8                 10.138.0.4   35.185.203.59   RUNNING
```

Then connect to your instance via SSH:

```shell
gcloud compute ssh --zone us-west1-b "${USER}-deepvariant-quickstart"
```

You should land at a shell prompt in your new instance!

NOTE: All of these steps can also be completed from the Cloud Console, if you
prefer. Consult [this
guide](https://cloud.google.com/compute/docs/quickstart-linux), but be sure to
choose Ubuntu 20.04 as your image, as DeepVariant has not been tested on other
Linux distributions.

For more information about getting started with Compute Engine, see:

*   [Compute Engine instance creation in `gcloud`
    manual](https://cloud.google.com/sdk/gcloud/reference/compute/instances/create)
*   [Reference to machine
    sizes/types](https://cloud.google.com/compute/docs/machine-types)
    
## More about the dataset
We have chosen to look at the HCC1143 cell line, which is a publicly available illumina whole genome sequencing data. The cell line was generated from a 52 year old caucasian woman with breast cancer tumor. Fastq files of both matched normal and tumor were preprocessed, subjected to GATK best practices. The bam files containing the reads for the cancer cell line and the matched normal consists of these 2 bam files. We chose to only look at reads from chromosome 17 as we wanted to start with a smaller dataset to test our pipeline. The genome sequence reads were aligned to the Human GRCh38 reference genome.

## More about Docker Containers 
To ensure that such tools used in the Nextflow pipeline can run on any machine without running into errors of uninstalled dependencies, Docker containers are used. We are able to encapsulate each process in a Docker container. The configurations needed in the container can be specified in a Docker image, and the image is used like a piece of instruction to build the Docker container. When a process is being run, it would be running in a container, which we can think of as an environment that is specially configured for that process. Therefore, there are no worries about the environment in the deployed machine affecting the execution of each process.

## More about the pipeline 

Some input files ar optional and if not given, they will be automatically created for the user during the preprocessing steps. If these are given, the preprocessing steps are skipped. For more information about preprocessing, please refer to the "INPUT PARAMETERS" section.

The worklow **accepts one reference genome and multiple BAM files as input**. The variant calling for the several input BAM files will be processed completely indipendently and will produce indipendent VCF result files. The advantage of this approach is that the variant calling of the different BAM files can be parallelized internally by Nextflow and take advantage of all the cores of the machine in order to get the results at the fastest.


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









    
    

    
    
