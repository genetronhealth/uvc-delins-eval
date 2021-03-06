This script contains the code used for evaluating UVC-delins (which is available at https://github.com/genetronhealth/uvc-delins-eval). 

The script main-delins-eval-set-vars.sh contains global filenames that should be set up properly before "running main-delins-eval.sh". 
This script should be be sourced from another script instead of being run on its own. 
Please set up these global filenames properly by modifying this script. 
The default locations used with ./install-dependencies.sh can be left unmodified. 

The script ./install-dependencies.sh installs the pre-requisite dependencies to set up for the evaluation. 
Please note that both JAVA8 and GATK4 are not installed by default as the user is highly likely to have already installed them. 
If not, then it is easy to install both tools anyway, and we highly recommend installing them. 
The usage of ./install-dependencies.sh is : ./install-dependencies.sh <keywords> <EVALROOT>
where keywords can be set to skip-ref to skip downloading-and-indexing the human reference genomes
and EVALROOT is the non-default root directory containing evaluation tools. 
The default location is specified in main-delins-eval-set-vars.sh

The script main-delins-eval.sh generates the actual code used for evaluating with different datasets.
Running ./main-delins-eval.sh &#36;{datadir} prints the shell code to evaluate the fastq.gz data at the directory 
&#36;{datadir} to stdout.
This repository does not contain any code for downloading the raw sequencing data because the best way to get the raw sequencing data is geolocation-specific. 
We do recommend, however, to download raw sequencing data by using this website to get the download command: https://sra-explorer.info/

Please note that, for the HCC1395/HCC1395BL dataset, the BAM files pre-aligned by BWA MEM and the BAM files pre-processed by GATK (i.e., prepared BAM files) generated by the script at https://github.com/genetronhealth/uvc-eval/blob/master/seqc2/eval-seqc2.sh (commit 222b8d3) should be used instead of the raw FASTQ files since BaseRecalibration using the entire BAM file is better.
Because it may take some time to generate the prepared BAM files, we uploaded the prepared BAM files at https://doi.org/10.5281/zenodo.6586710 for your convenience. 

