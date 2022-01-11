#!/usr/bin/env bash

EVALROOT="$I/uvc/eval/"
HS37D5=${EVALROOT}/datafiles/hs37d5.fa
GRCH38=${EVALROOT}/datafiles/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna

VT="${EVALROOT}/tools/vt-0.57721/vt"

java8=$I/software/jdk1.8.0_181/bin/java
gatk4lowmem="$java8 -Xmx4g -Djava.io.tmpdir=${EVALROOT}/systmp -jar ${EVALROOT}/tools/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar"

UVC=${I}/public/software/uvc-20211109/uvc1
UVCdelins=${I}/public/software/uvc-20211109/uvc-delins-20211223/uvcvcf-raw2delins-all.sh

FREEBAYES=${EVALROOT}/tools/freebayes-1.3.4-linux-static-AMD64
VARDICT_DIR=${EVALROOT}/tools/VarDict-1.8.3/bin/
INDELSEEK=${EVALROOT}/tools/indelseek/indelseek.pl
PINDEL=${EVALROOT}/tools/pindel-0.2.5b8/pindel
PINDEL2VCF=${EVALROOT}/tools/pindel-0.2.5b8/pindel2vcf

NORM_WITH_INDELPOST="${EVALROOT}/tools/norm-with-indelpost.py"

EGFR_DEL19='7:55242415-55242513'
EGFR_INS20='7:55248986-55249171'
ERBB2_INS20='17:37880979-37881164'
KRAS_ALL='12:25357723-25403870'
BRCA1_DEL='17:41196312-41277500'
BRCA2_DEL='13:32889645-32974405'
TP53_DEL=17:7565097-7590856 # https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000141510;r=17:7565097-7590856

datadir=$1
flag=$2

#TARGETS="${EGFR_DEL19},${EGFR_INS20},${ERBB2_INS20},${BRCA1_DEL},${BRCA2_DEL},${TP53_DEL}"

export PATH="${EVALROOT}/tools:${PATH}"

function myecho {
    echo 'for i in 1 2 3; do ' $@ ' ; if [ $? -eq 0 ] ; then break ; fi ; done ; if [ $? -ne 0 ]; then exit $? ; fi ;'
}


function eval12 {
    truthvcf=$1
    invcf=$2
    field=$3
    for region in $TARGETS2 ; do
        callvcf2=${invcf/.vcf/}.${region/:/-}.vcf
        summary=${invcf/.vcf/}.summary.${region/:/-}.txt
        echo "# STEP-TRUTH-EVAL-01"
        myecho bcftools view -v $VARTYPES -i "\"ALT != '*'\"" -r $region $invcf  \
            '|' bcftools norm -m-any -f ${HGREF} - \
            '|' $VT normalize -r ${HGREF} - \
            '|' bcftools norm -m-any -f ${HGREF} - \
            '|' bcftools view -i "\"ALT != '*'\"" -v $VARTYPES - \
            '>' ${callvcf2} 
        echo "# STEP-TRUTH-EVAL-02"
        myecho bcftools view $truthvcf -t $region \
            '|' python ${I}/uvc-delins/extdata/base-call-vcfs-to-prec-sens-list.py --base-vcf - --call-vcf ${callvcf2} --field ${field} \
            '>' ${summary} 
        
    done
}

DELINS_BED="${datadir}/delins.bed"
for region in $TARGETS2; do
    printf "${region}\n" | sed 's/:/\t/g' | sed 's/-/\t/g'
done > ${DELINS_BED} || true

for fq1 in $(ls ${datadir}/*_1.fastq.gz); do
    fq2=${fq1/_1.fastq.gz/_2.fastq.gz}
    rawbam=${fq1/_1.fastq.gz/_12.bam}
    rmdupbam=${fq1/_1.fastq.gz/_12.rmdup.bam}
    recalbam=${fq1/_1.fastq.gz/_12.rmdup_recal.bam}
    srr=$(echo $fq1 | awk -F"/" '{print $NF}' | awk -F "_" '{print $1}')
    
    if [ $(echo $datadir | grep -c SRP162370) -gt 0 ]; then # HCC1395 SEQC2-FDA
        HGREF=${GRCH38}
        dbsnp=${EVALROOT}/datafiles/GRCh38/GATK/All_20180418.vcf.gz
        tpref=$(echo $rawbam | awk -F"/" '{print $NF}' | awk -F"_" '{print $1}')
        TARGETS=chr22:$(echo $tpref | awk -F"-" '{print $3}')-$(echo $tpref | awk -F"-" '{print $4}') # SRR7890897-chr22-10714122-10714124_12.bam 
        VARTYPES=snps,mnps,other
    else
        HGREF=${HS37D5}
        dbsnp=${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf
        TARGETS=${EGFR_DEL19},${EGFR_INS20},${ERBB2_INS20} # ,${KRAS_ALL}
        VARTYPES=indels,mnps,other
        #TARGETS="${EGFR_DEL19},${EGFR_INS20},${KRAS_ALL}"
    fi
    TARGETS2=$(echo $TARGETS | sed 's/,/ /g')
    
    myecho bwa mem -t 24 -R "\"@RG\tID:${srr}.L001\tSM:${srr}\tLB:${srr}\tPL:ILLUMINA\tPM:UNKNOWN\tPU:${srr}.L001\"" "${HGREF}" $fq1 $fq2 \
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
        myecho $gatk4lowmem MarkDuplicates --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true -I ${rawbam} -M ${rmdupbam}.metrics -O ${rmdupbam} '&&' samtools index -@8 ${rmdupbam}
        myecho $gatk4lowmem BaseRecalibrator -I ${rmdupbam} --known-sites $dbsnp -O ${recalbam/.bam/.table} -R ${HGREF}
        myecho $gatk4lowmem ApplyBQSR -bqsr ${recalbam/.bam/.table} -I ${rmdupbam} -O ${recalbam} '&&' samtools index -@8 ${recalbam}
        inbam=${recalbam}
    fi
    
    if [ $(echo $datadir | grep -c SRP162370) -gt 0 ]; then # HCC1395 SEQC2-FDA
        FILTCMD="\" \""
    else
        FILTCMD="\"((STRLEN(REF) >= 4 || STRLEN(ALT) >= 4) || TYPE = 'mnp') && QUAL>=${QUALTHRES}\""
    fi
    
    resdir=${fq1/_1.fastq.gz/_12.resdir}
    mkdir -p -m 777 ${resdir} && chmod uga+s ${resdir}
    echo "mkdir -p -m 777 ${resdir} && chmod uga+s ${resdir}"
    uvcvcf=${resdir}/uvc.vcf.gz
    uvcdelins=${resdir}/uvc-hap
    myecho $UVC ${uvcargs} -f ${HGREF} -o ${uvcvcf} ${rawbam} -s ${srr} '2>' ${uvcvcf}.stderr
    myecho bash -evx $UVCdelins "${HGREF}" ${uvcvcf} ${uvcdelins} '2>' ${uvcdelins}.stderr
   
    truthvcf=${fq1/_1.fastq.gz/_12.uvc-truth.vcf} # need manual review
    callvcfgz=${uvcdelins}.merged-simple-delins.vcf.gz
    
    echo '# STEP-TRUTH-PREP '
    myecho bcftools view -v $VARTYPES -i "\"ALT != '*'\"" ${callvcfgz} ${TARGETS} \
        '|'  bcftools norm -m-any -f ${HGREF} - \
        '|' $VT normalize -r $HGREF -  \
        '|'  bcftools norm -m-any -f ${HGREF} - \
        '|'  bcftools view -v $VARTYPES -i "${FILTCMD}" -  \
        '>' ${truthvcf}
    eval12 ${truthvcf} ${callvcfgz} QUAL
    
    if [ $(echo $flag | grep generate-igv-html -c) -gt 0 ]; then
        bcftools view --no-header ${truthvcf} | while read -r line ; do
            chrom=$(echo $line | awk '{print $1}')
            ref=$(echo $line | awk '{print $4}')
            alt=$(echo $line | awk '{print $5}')
            pos=$(echo $line | awk '{print $2}')
            qual=$(echo $line | awk '{print $6}')
            python $I/scripts/igv-html/GetNTIgv.py -tb ${rawbam} -r ${HGREF} -p $chrom:$(($pos+1))  > $resdir/${qual}_${chrom}_${pos}_${ref}_${alt}_${srr}.html || true
        done
    fi
    
    callvcf=${resdir}/mutect2_tonly.vcf
    myecho $gatk4lowmem Mutect2 -R ${HGREF} -I ${inbam} -O $callvcf \
        '&&' bcftools view -Oz -o ${callvcf}.gz $callvcf \
        '&&' bcftools index ${callvcf}.gz
    eval12 ${truthvcf} ${callvcf}.gz INFO/TLOD
    
    if [ $(echo $flag | grep run-indelpost -c) -gt 0 ]; then
        bcftools view ${callvcf}.gz $TARGETS \
            | ${MY_PYTHON} ${NORM_WITH_INDELPOST} ${HGREF} ${inbam} \
            | bcftools norm -d none \
            | bcftools view -Oz -o ${callvcf}-indelpost.vcf.gz \
            && bcftools index -ft ${callvcf}-indelpost.vcf.gz 
        eval12 ${truthvcf} ${callvcf}-indelpost.vcf.gz INFO/TLOD | bash -vx
    fi

    callvcfgz=${resdir}/freebayes_tonly.vcf.gz
    myecho $FREEBAYES --min-alternate-fraction ${minFA} --pooled-continuous --min-alternate-count 2 -f ${HGREF} ${inbam} \
        '|' bcftools view - -Oz -o ${callvcfgz} \
        '&&' bcftools index ${callvcfgz}  
    eval12 ${truthvcf} ${callvcfgz} FORMAT/AD
    
    callvcfgz=${resdir}/vardict_tonly.vcf.gz
    myecho ${VARDICT_DIR}/VarDict -G ${HGREF} -f $minFA -N ${srr} -b ${inbam} -c 1 -S 2 -E 3 -g 4 ${DELINS_BED} \
        '|' Rscript ${VARDICT_DIR}/teststrandbias.R \
        '|' ${VARDICT_DIR}/var2vcf_valid.pl -N ${srr} -E -f $minFA \
        '|' bcftools view -Oz -o ${callvcfgz} \
        '&&' bcftools index ${callvcfgz}
    eval12 ${truthvcf} ${callvcfgz} FORMAT/AD
    
    callvcfgz=${resdir}/indelseek_tonly.vcf.gz
    myecho samtools view ${inbam} ${TARGETS} \
        '|' ${INDELSEEK} --refseq $HGREF \
        '|' bcftools view -fPASS -Oz -o ${callvcfgz} - \
        '&&' bcftools index ${callvcfgz}
    eval12 ${truthvcf} ${callvcfgz} QUAL
    
    callvcf=${resdir}/pindel_tonly.vcf
    callvcfgz=${resdir}/pindel_tonly.vcf.gz
    myecho printf "\"${inbam}\t200\t${resdir}/pindel_config.txt\n\" > ${resdir}/pindel_config.txt" \
        '&&' ${PINDEL} -f ${HGREF} -i ${resdir}/pindel_config.txt -o ${resdir}/${srr}_pindel
    pindel_files=""
    pindel_all=${resdir}/${srr}_pindel_ALL
    for suffix in BP CloseEndMapped D INT_final INV LI RP SI TD; do
        pindel_files="$pindel_files ${resdir}/${srr}_pindel_${suffix}"
    done
    myecho cat "$pindel_files > $pindel_all"
    myecho ${PINDEL2VCF} -r ${HGREF} -R 'hs37d5.fa' -d 'hs37d5' -p $pindel_all -v ${callvcf}
    myecho bcftools view -Oz -o ${callvcfgz} ${callvcf} '&&' bcftools index ${callvcfgz} 
    eval12 ${truthvcf} ${callvcfgz} FORMAT/AD
done

