#!/usr/bin/env bash

## This script download the human reference genomes, bioinformatics tools, and other data required for the evaluation of deletion-insertion calling performance.
## JAVA8 and GATK4 are not automatically downloaded and installed because the user may already have installed them
## Also, GRCh38-related files are not downloaded by default because we recommend using the pre-processed files at https://doi.org/10.5281/zenodo.6529871

set -vx

scriptdir=$(dirname $(which "${0}"))
source "${scriptdir}/main-delins-eval-set-vars.sh"

PAT=$1
if [ -n "${2}" ]; then
    EVALROOT="${2}"
fi

function getdata () {
    if [ $(echo "${PAT}" | grep -c skip-download) -eq 0 ]; then
        wget --tries 1000 $@
    fi
}

function prepdata () {
    fastafile=$(echo ${1} | sed 's;.gz$;;g')
    #zcat "${1}" > "${fastafile}"
    #bwa index $fastafile
    #samtools faidx $fastafile
    
    #samtools dict $fastafile > ${fastafile}.dict
    cp ${fastafile}.dict $(echo $fastafile | awk -F"." '{$(NF) = "" ; print $0}' | sed 's; $;;g').dict
}

# START-DOWNLOADING-INSTALLNG-TOOLS

mkdir -p "${EVALROOT}/tools/"
pushd "${EVALROOT}/tools/"
if [ $(echo "${PAT}" | grep -c skip-software) -eq 0 ]; then
    
    if [ $(echo "${PAT}" | grep -c skip-download) -gt 0 ]; then
        skipdown="skip-downloading-samtools,skip-downloading-bcftools"
    else
        skipdown=""
    fi
    
    if [ $(echo "${PAT}" | grep -c skip-clone ) -eq 0 ]; then
        git clone https://github.com/genetronhealth/uvc.git
    fi
    pushd uvc
    git checkout e49bf2112f1591925bac29fa655f4cc300d90f7f
    bash -evx ./install-dependencies.sh "$skipdown" && make && make deploy && cp bin/* ./
    popd

    if [ $(echo "${PAT}" | grep -c skip-clone ) -eq 0 ]; then
        git clone https://github.com/genetronhealth/uvc-delins.git
    fi
    pushd uvc-delins
    git checkout f2987d5
    cp -hr ../uvc/ext/ ./ && make && make deploy && cp bin/* ./
    popd

    getdata https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
    tar jxvf samtools-1.11.tar.bz2
    pushd samtools-1.11
    ./configure --prefix=`pwd -P` --disable-plugins --disable-libcurl --disable-s3 --disable-largefile
    make && make install
    cp samtools ../bin
    cp misc/wgsim ../bin
    popd

    getdata https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
    tar jxvf bcftools-1.11.tar.bz2
    pushd bcftools-1.11
    ./configure --prefix=`pwd -P` --disable-plugins --disable-libcurl --disable-s3 --disable-largefile
    make && make install
    cp bcftools ../bin
    popd

    getdata https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
    tar jxvf bwa-0.7.17.tar.bz2
    pushd bwa-0.7.17
    make
    cp bwa ../bin
    popd
    
    getdata https://github.com/atks/vt/archive/refs/tags/0.57721.tar.gz
    tar -xvf vt-0.57721.tar.gz
    pushd vt-0.57721
    make
    popd
    
    # Get the variant callers freebayes, VarDict, indelseek, pindel, and strelka. 
    # Please note that GATK4 is regularly updated. 
    # Therefore, we recommend the user to download his/her own lateset version of GATK. 
    # If exact reproducibility is desired, then please use the GATK version presented in the sourced script. 

    getdata https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-linux-static-AMD64.gz
    zcat freebayes-1.3.4-linux-static-AMD64.gz > freebayes-1.3.4-linux-static-AMD64
    
    getdata https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.3/VarDict-1.8.3.tar
    tar -xvf VarDict-1.8.3.tar

    if [ $(echo "${PAT}" | grep -c skip-clone ) -eq 0 ]; then
        git clone https://github.com/tommyau/indelseek.git
    fi
    pushd indelseek
    git checkout dba5cf8f84832cbe301ae2d5c4a89846bac94ab4
    popd

    getdata https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
    tar -xvf strelka-2.9.10.centos6_x86_64.tar.bz2
fi
popd

# START-DOWNLOADING-INSTALLNG-REFERENCE-GENOMES

mkdir -p "${EVALROOT}/datafiles/"
pushd "${EVALROOT}/datafiles/"
if [ $(echo "${PAT}" | grep -c skip-ref) -eq 0 ]; then
    getdata ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
    prepdata hs37d5.fa.gz

    # We recommend to use the pre-processed files at https://doi.org/10.5281/zenodo.6529871 instead of downloading this data
    if [ $(echo "${PAT}" | grep -c get-grch38) -gt 1 ]; then
        mkdir -p GRCh38 && pushd GRCh38
        getdata ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
        prepdata GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
        mkdir -p GATK && pushd GATK
        getdata ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
        zcat All_20180418.vcf.gz > All_20180418.vcf
        $gatk4lowmem IndexFeatureFile -I All_20180418.vcf
        popd
        popd
    fi

    getdata http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
    prepdata human_g1k_v37.fasta.gz  
    getdata ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
    zcat dbsnp_138.b37.vcf.gz > dbsnp_138.b37.vcf
fi
$gatk4lowmem IndexFeatureFile -I dbsnp_138.b37.vcf
popd

