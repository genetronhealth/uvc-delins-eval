# This code should be sourced. 
# Please modify these variables carefully accordingly to their real path locations. 
# Then, please run ./main-delins-eval.sh ${datadir} where ${datadir} is the directory containing the raw FASTQ files. 
# main-delins-eval.sh will source this script. 

if [ -z "${EVALROOT}" ]; then
    EVALROOT="$I/uvc/eval/"
fi

mkdir -p ${EVALROOT}/systmp/
mkdir -p ${EVALROOT}/datafiles/
mkdir -p ${EVALROOT}/public/software/
mkdir -p ${EVALROOT}/tools/

# human reference genome files
HS37D5=${EVALROOT}/datafiles/hs37d5.fa
GRCH38=${EVALROOT}/datafiles/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna

# common tools
java8=$I/public/software/jdk1.8.0_181/bin/java # The java path has to be manually set

# UVC tools
UVC=${I}/public/software/uvc-20211109/uvc1
UVCTN=${I}/public/software/uvc-20211109/uvcTN.sh
UVCdelins=${I}/public/software/uvc-delins-2022-0207/uvcvcf-raw2delins-all.sh

# other bioinformatics tools
gatk4lowmem="$java8 -Xmx4g -Djava.io.tmpdir=${EVALROOT}/systmp -jar ${EVALROOT}/tools/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar" # The gatk jar package has to be manually set
VT="${EVALROOT}/tools/vt-0.57721/vt"
FREEBAYES=${EVALROOT}/tools/freebayes-1.3.4-linux-static-AMD64
VARDICT_DIR=${EVALROOT}/tools/VarDict-1.8.3/bin/
INDELSEEK=${EVALROOT}/tools/indelseek/indelseek.pl
PINDEL=${EVALROOT}/tools/pindel-0.2.5b8/pindel
PINDEL2VCF=${EVALROOT}/tools/pindel-0.2.5b8/pindel2vcf
STRELKA2="${EVALROOT}/tools/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py"

