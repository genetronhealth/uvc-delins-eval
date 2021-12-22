#!/usr/bin/env bash

EVALROOT="$I/uvc/eval/"
HS37D5=${EVALROOT}/datafiles/hs37d5.fa

VT="${EVALROOT}/tools/vt-0.57721/vt"

java8=$I/software/jdk1.8.0_181/bin/java
gatk4lowmem="$java8 -Xmx4g -Djava.io.tmpdir=${EVALROOT}/systmp -jar ${EVALROOT}/tools/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar"
dbsnp=${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf.gz

UVC=${I}/public/software/uvc-20211109/uvc1
UVCdelins=${I}/public/software/uvc-20211109/uvc-hap-20211129/uvc-rawvcf2allvcfs.sh

FREEBAYES=${EVALROOT}/tools/freebayes-1.3.4-linux-static-AMD64
VARDICT_DIR=${EVALROOT}/tools/VarDict-1.8.3/bin/
INDELSEEK=${EVALROOT}/tools/indelseek/indelseek.pl

EGFR_DEL19='7:55242415-55242513'
EGFR_INS20='7:55248986-55249171'
ERBB2_INS20='17:37880979-37881164 '

TARGETS="${EGFR_DEL19},${EGFR_INS20},${ERBB2_INS20}"

export PATH="${EVALROOT}/tools:${PATH}"

function eval12 {
    truthvcf=$1
    invcf=$2
    field=$3
    for region in $EGFR_DEL19 $EGFR_INS20 $ERBB2_INS20; do
        callvcf=${invcf/.vcf/}.${region/:/-}.vcf
        summary=${invcf/.vcf/}.summary.${region/:/-}.txt
        echo bcftools view $invcf -v indels,mnps,other -i "\"ALT != '*'\"" -r $region \
            '|' bcftools norm - -m-any  \
            '|' $VT normalize -r ${HS37D5} - \
            '|' bcftools norm - -m-any  \
            '|' bcftools view -i "\"ALT != '*'\"" -v indels,mnps,other \
            '>' ${callvcf}
        echo bcftools view $truthvcf -t $region \
            '|' python ${I}/uvc-delins/extdata/base-call-vcfs-to-prec-sens-list.py --base-vcf - --call-vcf ${callvcf} --field ${field} \
            '>' ${summary}
    done
}

datadir=$1

DELINS_BED="${datadir}/delins.bed"
for region in $EGFR_DEL19 $EGFR_INS20 $ERBB2_INS20; do
    printf "${region}\n" | sed 's/:/\t/g' | sed 's/-/\t/g'
done > ${DELINS_BED} || true

for fq1 in $(ls ${datadir}/*_1.fastq.gz); do
    fq2=${fq1/_1.fastq.gz/_2.fastq.gz}
    rawbam=${fq1/_1.fastq.gz/_12.bam}
    rmdupbam=${fq1/_1.fastq.gz/_12.rmdup.bam}
    recalbam=${fq1/_1.fastq.gz/_12.rmdup_recal.bam}
    srr=$(echo $fq1 | awk -F"/" '{print $NF}' | awk -F "_" '{print $1}')
    echo bwa mem -t 24 -R "\"@RG\tID:${srr}.L001\tSM:${srr}\tLB:${srr}\tPL:ILLUMINA\tPM:UNKNOWN\tPU:${srr}.L001\"" "${HS37D5}" $fq1 $fq2 \
        '|' samtools view -bh1 \
        '|' samtools sort -o $rawbam \
        '&&' samtools index -@8 $rawbam
    
    if [ $(echo $fq1 | grep -c SRP268953) -gt 0 ]; then
        QUALTHRES=0
        uvcargs=" --dedup-flag 0x4 "
        mutect2args=" --mitochondria-mode true "
        minFA=0.001
        inbam=${rawbam}
    else
        QUALTHRES=30
        uvcargs=" "
        mutect2args=" "
        minFA=0.01
        echo $gatk4lowmem MarkDuplicates --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true -I ${rawbam} -M ${rmdupbam}.metrics -O ${rmdupbam} '&&' samtools index -@8 ${rmdupbam}
        echo $gatk4lowmem BaseRecalibrator -I ${rmdupbam} --known-sites $dbsnp -O ${recalbam/.bam/.table} -R ${HS37D5}
        echo $gatk4lowmem ApplyBQSR -bqsr ${recalbam/.bam/.table} -I ${rmdupbam} -O ${recalbam} '&&' samtools index -@8 ${recalbam}
        inbam=${recalbam}
    fi
    resdir=${fq1/_1.fastq.gz/_12.resdir}
    mkdir -p ${resdir}
    uvcvcf=${resdir}/uvc.vcf.gz
    uvcdelins=${resdir}/uvc-hap
    #$UVC ${uvcargs} -f ${HS37D5} -o ${uvcvcf} ${inbam} -s ${srr} 2> ${uvcvcf}.stderr  
    #$UVCdelins "${HS37D5}" ${uvcvcf} ${uvcdelins} 2> ${uvcdelins}.stderr
   
    truthvcf=${fq1/_1.fastq.gz/_12.uvc-truth.vcf} # need manual review
    callvcfgz=${uvcdelins}.merged-simple-delins.vcf.gz
    echo bcftools view ${callvcfgz} ${TARGETS} -v indels,mnps,other -i "ALT != '*'"  \
        '|'  bcftools norm - -m-any  \
        '|' $VT normalize -r $HS37D5 -  \
        '|'  bcftools norm - -m-any  \
        '|'  bcftools view -i "\"ALT != '*'  && (STRLEN(REF) >= 4 || STRLEN(ALT) >= 4) && QUAL>=${QUALTHRES}\"" -v indels,mnps,other \
        '>' ${truthvcf}
    eval12 ${truthvcf} ${callvcfgz} QUAL

    callvcf=${resdir}/mutect2_tonly.vcf
    echo $gatk4lowmem Mutect2 -R ${HS37D5} -I ${inbam} -O $callvcf \
        '&&' bcftools view -Oz -o ${callvcf}.gz $callvcf \
        '&&' bcftools index ${callvcf}.gz
    eval12 ${truthvcf} ${callvcf}.gz INFO/TLOD

    callvcfgz=${resdir}/freebayes_tonly.vcf.gz
    echo $FREEBAYES --min-alternate-fraction ${minFA} --pooled-continuous --min-alternate-count 2 -f ${HS37D5} ${inbam} \
        '|' bcftools view - -Oz -o ${callvcfgz} \
        '&&' bcftools index ${callvcfgz}  
    eval12 ${truthvcf} ${callvcfgz} FORMAT/AD
    
    callvcfgz=${resdir}/vardict_tonly.vcf.gz
    echo ${VARDICT_DIR}/VarDict -G ${HS37D5} -f $minFA -N ${srr} -b ${inbam} -c 1 -S 2 -E 3 -g 4 ${DELINS_BED} \
        '|' Rscript ${VARDICT_DIR}/teststrandbias.R \
        '|' ${VARDICT_DIR}/var2vcf_valid.pl -N ${srr} -E -f $minFA \
        '|' bcftools view -Oz -o ${callvcfgz} \
        '&&' bcftools index ${callvcfgz}
    eval12 ${truthvcf} ${callvcfgz} FORMAT/AD
   
    callvcfgz=${resdir}/indelseek_tonly.vcf.gz
    echo samtools view ${inbam} ${TARGETS} \
        '|' ${INDELSEEK} --refseq $HS37D5 \
        '|' bcftools view -fPASS -Oz -o ${callvcfgz} - \
        '&&' bcftools index ${callvcfgz}
    eval12 ${truthvcf} ${callvcfgz} QUAL
done

