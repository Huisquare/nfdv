# Incorporating Google's DeepVariant variant caller as a Nextflow pipeline on Google Cloud Platform

A ZB4171 Project, incorporating Google DeepVariant with Nextflow pipeline.  

Group Members:  
Vanessa Tan Li Xuan  
Yan Zihao  
Li Huihui  

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
  Support containerisation using Docker
  Portable
  Parallelisation
  Continuous checkpointing
  All the intermediate results produced during the pipeline execution are automatically tracked

### Google Cloud Platform (GCP) 
Firstly, new customers are offered $300 free credits to run, test and deploy programs on Google Cloud. As DeepVariant is the brainchild of Google Life Sciences, its workflows can be run easily on Google Cloud Platform with optimised configurations. There is also support for running Nextflow pipelines with Google Life Sciences API.

## Motivation
DeepVariant specifies that it only takes in bam or cram files and its indexed bam file for aligned reads. It requires fasta files and its index file for reference genomes. There is a need to provide input files in such format to run DeepVariant. In addition, when multiple bam files are passed into DeepVariant, it runs them one by one. Also, DeepVariant requires high computing power, therefore running it on local machine may not be an ideal choice.

#### Preprocessing & Parallelisation 
The Nextflow pipeline aims to automatically handle the creation of some extra needed index files such as the fai and the bai files, which are needed as inputs for DeepVariant. This files were normally produced by users. The support for Docker in Nextflow allows us to bundle each preprocessing step with an individual Docker container. In addition, incorporating DeepVariant into a Nextflow pipeline allows DeepVariant to leverage on the parallelisation that Nextflow offers. This allows variant calling to be performed at the same time on multiple bam files. Finally, the nextflow pipeline can be deployed onto cloud to use the high computing power offered by cloud computing companies.

## Dependencies

[Nextflow](https://www.nextflow.io/) \
[Docker](https://www.docker.com/) \
[Google Cloud](https://cloud.google.com/)

## About the pipeline

Our pipeline only supports whole genome sequencing data. The workflow **accepts one reference genome and one folder containing multiple BAM files as input**. 
The parallelisation supported by Nextflow allows the multiple bam files in the input folder to be processed independently and simultaneously during variant calling. Through parallelisation, Nextflow allows the most optimal usage of computational resources by utilising all available CPU cores.

## Using the pipeline

A run on whole genome sequencing data looks like this: 
```
git clone https://github.com/Huisquare/nfdv
cd nfdv
nextflow run main.nf --fasta path/to/fastaFile --bam_folder path/to/bamFolder
```

## Checking out the pipeline

If you just want to check out the pipeline, you can run the below code. In this run, the fasta files used are from chr20 from hg19 genome assembly and the bam files used are small bam files for fast processing. 

```
git clone https://github.com/Huisquare/nfdv
cd nfdv
nextflow run main.nf --test
```

## Input parameters 

### Reference genome (fasta) input

An input fasta file for the reference genome is required. The path to the fasta file should be specified in the command using:
```
--fasta path/to/fastaFile
```

To allow DeepVariant to run faster, it requires some indexed and compressed versions of the reference genome (in fasta) and the alignment files (in bam). There is a preprocessing step in the pipeline to allow for those files to be produced using samtools and bgzip. If both the files are already at your disposal (.fa.fai, .fa.gz, .fa.gz.fai, .fa.gz.gzi for fasta and .bai for bam), you can specify them when running the command. However, they are optional inputs.

For fasta related optional inputs:
```
--fai   "/path/to/myGenome.fa.fai"
--fastagz "/path/to/myGenome.fa.gz"
--gzfai  "/path/to/myGenome.fa.gz.fai"
--gzi  "/path/to/myGenome.fa"
```

### Alignment file (bam) input 

An input bam folder (containing all the bam files to be processed) is required. The path to the bam folder should be specified in the command using:
```
--bam_folder path/to/bamFolder
```

If the bam_folder specified contain other bam files that are not to be used as input, a file prefix for the bam files to be used can be specified using the below command tag: 
```
--bam_file_prefix prefix_of_bam_files
```

If optional indexed bam inputs for the bam files are present, they must reside in the same input bam folder and have the same prefix as the corresponding bam file that it is indexing. (e.g. file.bam and file.bam.bai):
```
--getBai "true"
```

### Advanced parameters options

#### Number of CPUs 

In the pipeline, the makeExamples process is able to be parallelized and the user can define how many CPUs are to be used in this process. By default, **all the CPUs** of the machine are used in the process.

The **make_example** process can be internally parallelized and it can be defined how many cpus should be assigned to this process.
By default all the cpus of the machine are used.

```
-- numCores int_number_of_cpus_to_use
```

## Google Cloud support
Firstly, a Google Cloud account is required. Go to [Google Cloud Platform](https://cloud.google.com/gcp) and create an account.

To use Gloogle Cloud to run the code, a `nextflow-service-account` is required. Follow the below steps to create one:

(Adapted from [Google Life Sciences Nextflow Guide](https://cloud.google.com/life-sciences/docs/tutorials/nextflow))

### Creating a Service Account
Create a service account using Cloud Console:

In the Cloud Console, go to the Service Accounts page.

[Go to Service Accounts page](https://console.cloud.google.com/projectselector2/iam-admin/serviceaccounts)

Click **Create service account**.

In the **Service account name** field, enter `nextflow-service-account`, and then click **Create**.

In the **Grant this service account access to project** section, add the following roles from the **Select a role** drop-down list:

    - Cloud Life Sciences Workflows Runner
    - Service Account User
    - Service Usage Consumer
    - Storage Object Admin
Click **Continue**, and then click **Done**.

In the [Service Accounts page](https://console.cloud.google.com/projectselector2/iam-admin/serviceaccounts), find the service account you created. In the service account's row, click the **More** (3 dots) button, and then click **Manage keys**.

On the **Keys** page, click **Add key**, and then click **Create new key**.

Select **JSON** for the **Key type** and click **Create**.

A JSON file that contains your key downloads to your computer.

### Providing credentials to your application
You can provide authentication credentials to your application code or commands by setting the environment variable `GOOGLE_APPLICATION_CREDENTIALS` to the path of the JSON file that contains your service account key.

The following steps show how to set the `GOOGLE_APPLICATION_CREDENTIALS` environment variable:

1. Open Cloud Shell (or Terminal of Virtual Machine).

2. From the Cloud Shell More menu, select Upload file, and select the JSON key file you created. The file is uploaded to the home directory of your Cloud Shell instance.

Confirm that the uploaded file is in your present directory and confirm the filename by running the following command:

```
ls
```
3. Set the credentials, replacing KEY_FILENAME.json with the name of your key file.

```
export GOOGLE_APPLICATION_CREDENTIALS=${PWD}/KEY_FILENAME.json
```

### Installing Nextflow in Cloud Shell or Virtual Machine

Go to the terminal and run:

```
export NXF_VER=20.10.0
export NXF_MODE=google
curl https://get.nextflow.io | bash
```

### Installing Docker in Cloud Shell or Virtual Machine

Run the following code in the terminal to install Docker: 

```
sudo apt-get -qq -y install \
  apt-transport-https \
  ca-certificates \
  curl \
  gnupg-agent \
  software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository \
  "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) \
  stable"
sudo apt-get -qq -y update
sudo apt-get -qq -y install docker-ce
```

You have the environment set up to run the pipeline on Google Cloud now!

## More about Docker Containers 
To ensure that such tools used in the Nextflow pipeline can run on any machine without running into errors of uninstalled dependencies, Docker containers are used. We are able to encapsulate each process in a Docker container. The configurations needed in the container can be specified in a Docker image, and the image is used like a piece of instruction to build the Docker container. When a process is being run, it would be running in a container, which we can think of as an environment that is specially configured for that process. Therefore, there are no worries about the environment in the deployed machine affecting the execution of each process.
The Docker files that we wrote for each process can be found in the **Dockerfiles** folder. The images of the Docker files can be found on Dockerhub. \
Links to Docker images on Dockerhub: \
[htslib-and-samtools](https://hub.docker.com/repository/docker/huisquare/htslib-and-samtools) \
[samtools-config](https://hub.docker.com/repository/docker/huisquare/samtools-config) \
[vcftools-config](https://hub.docker.com/repository/docker/huisquare/vcftools-config) 

## More about the dataset
We have chosen to look at the HCC1143 cell line, which is a publicly available illumina whole genome sequencing data. The cell line was generated from a 52 year old caucasian woman with breast cancer tumor. Fastq files of both matched normal and tumor were preprocessed, subjected to GATK best practices. The bam files containing the reads for the cancer cell line and the matched normal consists of these 2 bam files. We chose to only look at reads from chromosome 17 as we wanted to start with a smaller dataset to test our pipeline. The genome sequence reads were aligned to the Human GRCh38 reference genome.

## Acknowledgements
We referenced similar pipelines developed by [lifebit.ai](https://github.com/lifebit-ai/DeepVariant) and [nf-core](https://github.com/nf-core/deepvariant) when building our pipeline.

We were able to run our pipeline on Google Cloud due to the USD$300 free credits that we received as new users in Google Cloud. 

## Improvements from current similar pipelines
The DeepVariant model we used is the 1.2.0 version, which is a huge advancement from v0.6.1 used by [lifebit.ai](https://github.com/lifebit-ai/DeepVariant) and v1.0 used by [nf-core](https://github.com/nf-core/deepvariant).

The model used is enclosed within the DeepVariant docker container, instead of using a model stored on cloud. Our pipeline is therefore more efficient in this aspect as it does not need to download an additional trained DeepVariant model from cloud storage.
