#!/bin/bash
# PRS
#$ -N prepare_UKB_base
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------
# Prepare public GWAS as base datasets for PRS or other sumstats-based downstream analyses
#---------------------------------------------------------------------- 
## following: https://choishingwan.github.io/PRS-Tutorial/base/


#----------------------------------------------------------------------
# Define variables
#----------------------------------------------------------------------
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/sumstats/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/UKB_BIG40/
mkdir -p ${working_dir}
#----------------------------------------------------------------------
# get arguments
base_file=${1}
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/
#---------------------------------------------------------------------- 

cd ${working_dir}
if [ ! -f ${base_pheno}.QC.txt.gz ]
 then
# check structure of files, one file
#chr rsid pos a1 a2 beta se pval(-log10)
#01 rs367896724 10177 A AC -0.014576 0.014648 0.49525
#01 rs201106462 10352 T TA -0.011981 0.015092 0.3693
#01 rs534229142 10511 G A -0.04015 0.20085 0.074913

# extract phenotype name from file
base_pheno=$(echo ${base_file} | sed -e 's/.txt.gz//g' )

# QC base file - on MAF (col4) and INFO, if available: available in different file: variants.txt.gz
##zless variants.txt.gz
#chr rsid pos a1 a2 af info
#01 rs367896724 10177 A AC 0.400522 0.465752
#01 rs201106462 10352 T TA 0.39659 0.44294
#01 rs534229142 10511 G A 0.0013354 0.463
if [ ! -f ${primary_dir}variants.qc.txt ]
then
 zcat ${primary_dir}/../variants.txt.gz | awk 'NR==1 || ($6 > 0.01) && ($6 < 0.99) && ($7>0.7) {print}' > ${primary_dir}variants.qc.txt
fi

# keep only variants after QC
awk 'NR==FNR{a[$1,$2,$3,$4,$5]=$0;next}($1,$2,$3,$4,$5) in a{print $0}' variants.qc.txt <(zcat ${primary_dir}${base_file}) > ${base_pheno}.tmp

##Summary stats downloads
#3,935 files, one file per IDP, containing: chr rsid pos a1 a2 beta se pval(-log10). 
# **The beta coefficient is in the direction of a2.**

# Extract and rename header to match prsice script, and fix chromosome names
## convert -log10p -> p 
#-snp snpid --chr chr --bp pos --A1 a_1 --A2 a_2 --stat beta --se se --pvalue P
awk '(NR==1){print "chr rsid pos a_1 a_2 beta se pval(-log10) P"}(NR>1){P=10^(-$8)}(NR>1){print $1,$2,$3,$4,$5,$6,$7,$8,P}' ${base_pheno}.tmp | \
 sed -e 's/^01/1/g'  -e 's/^02/2/g' -e 's/^03/3/g' -e 's/^04/4/g' \
 -e 's/^05/5/g' -e 's/^06/6/g'  -e 's/^07/7/g' -e 's/^08/8/g' -e 's/^09/9/g'  -e 's/^0X/X/g' > ${base_pheno}.txt

# Ambiguous SNPs
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' ${base_pheno}.txt | \
    gzip > ${base_pheno}.nonamb.txt.gz

# Identify duplicated SNPs
gunzip -c ${base_pheno}.nonamb.txt.gz | awk '{ print $2}' | sort | uniq -d > ${base_pheno}.duplicated.snp
# Remove duplicated SNPs from base file
echo "Remove duplicated SNPs from base file and save clean file"
gunzip -c ${base_pheno}.nonamb.txt.gz | grep -vf ${base_pheno}.duplicated.snp | gzip - > ${base_pheno}.QC.txt.gz

zcat ${base_pheno}.QC.txt.gz | head -n 5
#rsid    chr     pos     a_1     a_0     EAF     beta    se      P
#rs2883059       3       49902160        T       C       0.5746  0.026   0.003   7.149e-25
#rs3796386       3       49899795        A       G       0.4235  -0.026  0.003   8.883e-25
ls -lh ${base_pheno}.QC.txt.gz

# The ${base_pheno}.QC.txt.gz base data are now ready for using in downstream analyses.

# clean intermediate files
rm ${base_pheno}.tmp ${base_pheno}.txt ${base_pheno}.nonamb.txt.gz ${base_pheno}.duplicated.snp


fi

#---------------------------------------------------------------------- 