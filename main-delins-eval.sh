#!/usr/bin/env bash

scriptdir=$(dirname $(which $0))
source "${scriptdir}/main-delins-eval-set-vars.sh"

EGFR_DEL19='7:55242415-55242513'
EGFR_INS20='7:55248986-55249171'
ERBB2_INS20='17:37880979-37881164'
KRAS_ALL='12:25357723-25403870'
BRCA1_DEL='17:41196312-41277500'
BRCA2_DEL='13:32889645-32974405'
TP53_DEL=17:7565097-7590856 # https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000141510;r=17:7565097-7590856

datadir=$1
flag=$2

if [ "${1}" == '-h' -o "${1}" == '--help' -o "${1}" == '' ]; then
    printf "Usage: ${0} \${datadir} [optional \${flag}]\n"
    printf "  \${datadir} : the directory containing downloaded FASTQ data files. The \${datadir} path string should consist of only alpha-numeric characters. \n"
    printf "  \${flag}    : the flag that is the empty string by default. When properly set, this flag can generate additional results. \n"
    printf "  The stdout output is a script (let us name it as \${stdout}-script.sh) that runs BWA MEM, different variant callers, variant normalization, and variant-call performance evaluation. \n"
    printf "    The user can run (cat \${stdout}-script.sh | grep \${keyword} | bash) to run only certain code containing the \${keyword}. \n"
    printf "    For example, the command (cat \${stdout}-script.sh | grep STEP-TRUTH-EVAL-01 -A1) runs the step associated with TRUTH-EVAL-01. \n"
    printf "    In the end, \${stdout}-script.sh generates \${datadir}/*.resdir/*.summary.*.txt (denoted as summary files) as the final results labeling each call\n"
    printf "    If max_fscore is smaller than 1 in a summary file, then we check the value of n_base_[snv|indel|delins] to determine which type of ground-truth variant is missed as only one type of variant is in the base-line ground truth. \n"
    exit 0
fi

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
    fq2=${fq1/%_1.fastq.gz/_2.fastq.gz}
    rawbam=${fq1/%_1.fastq.gz/_12.bam}
    rmdupbam=${fq1/%_1.fastq.gz/_12.rmdup.bam}
    recalbam=${fq1/%_1.fastq.gz/_12.rmdup_recal.bam}
    srr=$(echo $fq1 | awk -F"/" '{print $NF}' | awk -F "_" '{print $1}')
    
    if [ $(echo $datadir | grep -cP "SRP162370|SRR7890887") -gt 0 ]; then # HCC1395 SEQC2-FDA
        HGREF=${GRCH38}
        dbsnp=${EVALROOT}/datafiles/GRCh38/GATK/All_20180418.vcf.gz
        tpref=$(echo $rawbam | awk -F"/" '{print $NF}' | awk -F"_" '{print $1}')
        TARGETS=chr22:$(echo $tpref | awk -F"-" '{print $3}')-$(echo $tpref | awk -F"-" '{print $4}') # SRR7890897-chr22-10714122-10714124_12.bam 
        VARTYPES=snps,mnps,other
    elif [ $(echo $datadir | grep -c HNF4A) -gt 0 ]; then # HCC1395 SEQC2-FDA
        HGREF=${EVALROOT}/datafiles/Homo_sapiens_assembly19.fasta # https://doi.org/10.1016/j.omtn.2021.07.016 supplementary sequence S1
        dbsnp=${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf
        TARGETS="20" # ,${KRAS_ALL}
        VARTYPES=snps,mnps,other
    elif [ $(echo $datadir | grep -c PRJNA688630) -gt 0 ]; then # HCC1395 SEQC2-FDA
        HGREF=${datadir}/CABE-base-editor.fna # https://doi.org/10.1016/j.omtn.2021.07.016 supplementary sequence S1
        dbsnp=${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf
        TARGETS="" # ,${KRAS_ALL}
        VARTYPES=snps,indels,mnps,other
    else
        HGREF=${HS37D5}
        dbsnp=${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf
        TARGETS=${EGFR_DEL19},${EGFR_INS20},${ERBB2_INS20} # ,${KRAS_ALL}
        VARTYPES=indels,mnps,other
        #TARGETS="${EGFR_DEL19},${EGFR_INS20},${KRAS_ALL}"
    fi
    TARGETS2=$(echo $TARGETS | sed 's/,/ /g')
    
    if [ $(echo $datadir | grep -cP "SRR7890887|HNF4A") -eq 0 ]; then
        myecho bwa mem -t 24 -R "\"@RG\tID:${srr}.L001\tSM:${srr}\tLB:${srr}\tPL:ILLUMINA\tPM:UNKNOWN\tPU:${srr}.L001\"" "${HGREF}" $fq1 $fq2 \
            '|' samtools view -bh1 \
            '|' samtools sort -o $rawbam \
            '&&' samtools index -@8 $rawbam
    fi
    if [ $(echo $fq1 | grep -cP "SRP268953|PRJNA688630") -gt 0 ]; then
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
        if [ $(echo $datadir | grep -cP "SRR7890887|HNF4A") -gt 0 ]; then
            myecho echo NOTE: ${recalbam} is the bam aligned to only chr22 derived from the SRR7890887 reads aligned to all CRCh38 chromosomes. \
                If ${recalbam} does not exist, please generate it manually then apply the MarkDuplicates, ..., ApplyBQSR for all the reads of SRR7890887
            IS_GATK_SKIPPED=TRUE
            if [ $(echo $datadir | grep -cP "HNF4A") -gt 0 ]; then
                recalbam="${rawbam}"
            fi
        else
            myecho $gatk4lowmem MarkDuplicates --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true -I ${rawbam} -M ${rmdupbam}.metrics -O ${rmdupbam} '&&' samtools index -@8 ${rmdupbam}
            myecho $gatk4lowmem BaseRecalibrator -I ${rmdupbam} --known-sites $dbsnp -O ${recalbam/.bam/.table} -R ${HGREF}
            myecho $gatk4lowmem ApplyBQSR -bqsr ${recalbam/.bam/.table} -I ${rmdupbam} -O ${recalbam} '&&' samtools index -@8 ${recalbam}
        fi
        inbam=${recalbam}
    fi
    rawbam_n=$(ls ${rawbam/%.bam/.normal.dir}/*.bam) || true
    inbam_n=$(ls ${inbam/%.bam/.normal.dir}/*.bam) || true
    
    if [ $(echo $datadir | grep -cP "SRP162370|SRR7890887") -gt 0 ]; then # HCC1395 SEQC2-FDA
        FILTCMD="\" \""
    elif [ $(echo $datadir | grep -cP "HNF4A") -gt 0 ]; then
        FILTCMD="\"QUAL>=${QUALTHRES}\""
    else
        FILTCMD="\"((STRLEN(REF) >= 4 || STRLEN(ALT) >= 4) || TYPE = 'mnp') && QUAL>=${QUALTHRES}\""
    fi
    
    resdir=${fq1/_1.fastq.gz/_12.resdir}
    mkdir -p -m 777 ${resdir} && chmod uga+s ${resdir}
    echo "mkdir -p -m 777 ${resdir} && chmod uga+s ${resdir}"
    uvcvcf=${resdir}/uvc.vcf.gz
    uvcdelins=${resdir}/uvc-hap
    if [ $(echo $rawbam_n | grep -c ".bam") -gt 0 ]; then
        nsrr=$(echo $rawbam_n | awk -F"/" '{print $NF}' | awk -F "_" '{print $1}')
        myecho $UVCTN ${HGREF} ${rawbam} ${rawbam_n}  ${uvcvcf} ${srr},${nsrr} ${uvcargs} 
        myecho bash -evx $UVCdelins "${HGREF}" ${uvcvcf} ${uvcdelins} '2>' ${uvcdelins}.stderr
    else
        myecho $UVC ${uvcargs} -f ${HGREF} -o ${uvcvcf} ${rawbam} -s ${srr} '2>' ${uvcvcf}.stderr
        myecho bash -evx $UVCdelins "${HGREF}" ${uvcvcf} ${uvcdelins} '2>' ${uvcdelins}.stderr
    fi
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
            #pos=$(echo $line | awk '{print $2}')
            qual=$(echo $line | awk '{print $6}')
            pos=$(echo $line | awk '{if (length($4) > length($5)) {$2 = $2 + 1} } { print $1, $2, $4, $5} ' | awk '{print $2}')
            if [ $(echo $rawbam_n | grep -c ".bam") -gt 0 ]; then
                echo python $I/scripts/igv-html/GetNTIgv.py -tb ${rawbam} -r ${HGREF} -p $chrom:$pos -nb ${rawbam_n} ' > ' $resdir/${qual}_${chrom}_${pos}_${ref}_${alt}_${srr}.html #|| true
            else
                echo python $I/scripts/igv-html/GetNTIgv.py -tb ${rawbam} -r ${HGREF} -p $chrom:$pos ' > ' $resdir/${qual}_${chrom}_${pos}_${ref}_${alt}_${srr}.html
            fi
        done
    fi
    
    if [ $(echo $rawbam_n | grep -c ".bam") -gt 0 ]; then
        tsname=$(samtools view -H ${inbam}   | grep "^@RG" | sed 's/\t/\n/g' |grep "^SM:" | awk -F":" '{print $2}' )
        nsname=$(samtools view -H ${inbam_n} | grep "^@RG" | sed 's/\t/\n/g' |grep "^SM:" | awk -F":" '{print $2}' )
        callvcf=${resdir}/mutect2_tnpaired.vcf
        myecho $gatk4lowmem Mutect2 -R ${HGREF} -I ${inbam} -I ${inbam_n} --tumor ${tsname} --normal ${nsname} -O $callvcf \
            '&&' bcftools view -Oz -o ${callvcf}.gz $callvcf \
            '&&' bcftools index ${callvcf}.gz
    else 
        callvcf=${resdir}/mutect2_tonly.vcf
        myecho $gatk4lowmem Mutect2 -R ${HGREF} -I ${inbam} -O $callvcf \
            '&&' bcftools view -Oz -o ${callvcf}.gz $callvcf \
            '&&' bcftools index ${callvcf}.gz
    fi
    eval12 ${truthvcf} ${callvcf}.gz INFO/TLOD
    
    callvcf=${resdir}/strelka2_tnpaired.vcf
    if [ $(echo $rawbam_n | grep -c ".bam") -gt 0 ]; then
        tsname=$(samtools view -H ${inbam}   | grep "^@RG" | sed 's/\t/\n/g' |grep "^SM:" | awk -F":" '{print $2}' )
        nsname=$(samtools view -H ${inbam_n} | grep "^@RG" | sed 's/\t/\n/g' |grep "^SM:" | awk -F":" '{print $2}' )
        myecho "rm -r \"${callvcf}.rundir/\" || true"
        myecho $STRELKA2 --referenceFasta="${HGREF}" --runDir="${callvcf}.rundir" --tumorBam="${rawbam}" --normalBam="${rawbam_n}" \
            '1>' "${callvcf}.rundir-config.stdout"  
        myecho "${callvcf}.rundir/runWorkflow.py" -j 8 -m local \
            '1>' "${callvcf}.rundir/runWorkflow.stdout" \
            '2>' "${callvcf}.rundir/runWorkflow.stderr"
        # '|' awk "'OFS=\"\t\" { if ($0 !~ \"^#\") {$9=\"GT:\"$9; $10=\"0/1:\"$10; $11=\"0/1:\"$11; print; } else { print;}}'" \
        myecho bcftools concat -a "${callvcf}.rundir/results/variants/somatic.snvs.vcf.gz" "${callvcf}.rundir/results/variants/somatic.indels.vcf.gz" \
                '|' bcftools view -Oz -o "${callvcf}.gz" \
                '&&' bcftools index -f "${callvcf}.gz"
        eval12 ${truthvcf} ${callvcf}.gz INFO/SomaticEVS
    fi
    
    callvcfgz=${resdir}/freebayes_tonly.vcf.gz
    myecho $FREEBAYES --min-alternate-fraction ${minFA} --pooled-continuous --min-alternate-count 2 -f ${HGREF} ${inbam} \
        '|' bcftools view - -Oz -o ${callvcfgz} \
        '&&' bcftools index ${callvcfgz}  
    eval12 ${truthvcf} ${callvcfgz} FORMAT/AD
    
    callvcfgz=${resdir}/vardict_tonly.vcf.gz
    myecho "export JAVA_HOME=$I/software/jdk1.8.0_181/bin && export PATH=${JAVA_HOME}/bin:${PATH}" \
        '&&' ${VARDICT_DIR}/VarDict -G ${HGREF} -f $minFA -N ${srr} -b ${inbam} -c 1 -S 2 -E 3 -g 4 ${DELINS_BED} \
        '|' Rscript ${VARDICT_DIR}/teststrandbias.R \
        '|' ${VARDICT_DIR}/var2vcf_valid.pl -N ${srr} -E -f $minFA \
        '|' bcftools view -Oz -o ${callvcfgz} \
        '&&' bcftools index ${callvcfgz}
    eval12 ${truthvcf} ${callvcfgz} FORMAT/AD
    
    if [ $(echo $datadir | grep -c SRR7890887) -eq 0 ]; then
    callvcfgz=${resdir}/indelseek_tonly.vcf.gz
    myecho samtools view ${inbam} ${TARGETS} \
        '|' ${INDELSEEK} --refseq $HGREF \
        '|' bcftools view -fPASS -Ou - \
        '|' bcftools sort -Oz -o ${callvcfgz} - \
        '&&' bcftools index ${callvcfgz}
    eval12 ${truthvcf} ${callvcfgz} QUAL
    fi
    
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

