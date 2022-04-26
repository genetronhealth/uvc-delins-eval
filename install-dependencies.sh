#!/usr/bin/env bash

## This script download the human reference genomes, bioinformatics tools, and other data required for the evaluation of deletion-insertion calling performance.
## JAVA8 and GATK4 are not automatically downloaded and installed because the user may already have installed them

set -vx

scriptdir=$(dirname $(which "${0}"))
source "${scriptdir}/main-delins-eval-set-vars.sh"

if [ -n "${2}" ]; then
    EVALROOT="${2}"
fi

# START-DOWNLOADING-INSTALLNG-REFERENCE-GENOMES

mkdir -p "${EVALROOT}/datafiles/"
pushd "${EVALROOT}/datafiles/"
if [ $(echo "${1}" | grep -c skip-ref) -eq 0 ]; then
# Download and decompress (extract) GRCh37
wget --tries 1000 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz
# Index GRCh37
bwa index hs37d5.fa
samtools faidx hs37d5.fa
samtools dict hs37d5.fa

# Download and decompress (extract) GRCh38
wget --tries 1000 ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
# Index GRCh38
bwa index GCA_000001405.15_GRCh38_full_analysis_set.fna
samtools faidx GCA_000001405.15_GRCh38_full_analysis_set.fna
samtools dict GCA_000001405.15_GRCh38_full_analysis_set.fna
fi
popd

# START-DOWNLOADING-INSTALLNG-TOOLS

mkdir -p "${EVALROOT}/tools/"
pushd "${EVALROOT}/tools/"

# installing UVC (may take some time)
git clone https://github.com/genetronhealth/uvc.git
pushd uvc
git checkout 8610b54
bash -evx ./install-dependencies.sh && make && make deploy && cp bin/* ./
popd

# installing UVC-delins
git clone https://github.com/genetronhealth/uvc-delins.git
pushd uvc-delins
git checkout f2987d5
cp -rs ../uvc/ext/ ./ && make && make deploy && cp bin/* ./
popd

### install samtools-1.11
wget --tries 1000 https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar jxvf samtools-1.11.tar.bz2
pushd samtools-1.11
./configure --prefix=`pwd -P` --disable-plugins --disable-libcurl --disable-s3 --disable-largefile # which will make some CRAM files produced elsewhere unreadable
make && make install
cp samtools ../bin
cp misc/wgsim ../bin
popd

### install bcftools-1.11
wget --tries 1000 https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
tar jxvf bcftools-1.11.tar.bz2
pushd bcftools-1.11
./configure --prefix=`pwd -P` --disable-plugins --disable-libcurl --disable-s3 --disable-largefile # --disable-bz2 --disable-lzma # which will make some CRAM files produced elsewhere unreadable
make && make install
cp bcftools ../bin
popd

### install bwa-0.7.17
wget --tries 1000 https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar jxvf bwa-0.7.17.tar.bz2
pushd bwa-0.7.17
make
cp bwa ../bin
popd

# Download VT tools
wget --tries 1000 https://github.com/atks/vt/archive/refs/tags/0.57721.tar.gz
tar -xvf 0.57721.tar.gz
pushd vt-0.57721
make
popd

# Get the variant callers freebayes, VarDict, indelseek, pindel, and strelka. Please note that GATK4 is regularly updated. Therefore, we recommend the user to download his/her own lateset version of GATK

wget --tries 1000 https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-linux-static-AMD64.gz
gzip -d freebayes-1.3.4-linux-static-AMD64.gz

wget --tries 1000 https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.3/VarDict-1.8.3.tar
tar -fvx VarDict-1.8.3.tar

wget --tries 1000 https://github.com/tommyau/indelseek/archive/refs/heads/master.zip
unzip master.zip

wget --tries 1000 https://github.com/genome/pindel/archive/refs/tags/v0.2.5b8.zip
unzip v0.2.5b8.zip

wget --tries 1000 https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xvf strelka-2.9.10.centos6_x86_64.tar.bz2

popd

