#!/bin/bash
# PRS
#$ -N prepare_ENIGMA_base
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------
# Prepare public GWAS as base datasets for PRS analyses
#---------------------------------------------------------------------- 
## following: https://choishingwan.github.io/PRS-Tutorial/base/


#----------------------------------------------------------------------
# Define variables
#----------------------------------------------------------------------
#primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/ENIGMA3/enigma.ini.usc.edu/downloads/ENIGMA3_withGlobal/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/ENIGMA3/
mkdir -p ${working_dir}
#----------------------------------------------------------------------
# get arguments
base_file=${1}
base_pheno=${2}
#---------------------------------------------------------------------- 
# From ENIGMA3 study, Grasby et al. 2020
## http://enigma.ini.usc.edu/downloads/
#---------------------------------------------------------------------- 

cd ${working_dir}
if [ ! -f ${working_dir}${base_pheno}.QC.gz ]
 then
# check structure of files, one file
#zless ENIGMA3_mixed_se_wSA_Mean_bankssts_surfavg_20190429.txt.gz
#SNP A1 A2 FREQ1 BETA1 SE P N MARKER CHR BP
#rs28527770 t c 0.8564 1.4212 1.4089 0.3131 24661 1:751756 1 751756
#rs3094315 a g 0.8220 -0.0252 1.2187 0.9835 26492 1:752566 1 752566

# extract phenotype name from file
#base_pheno=$(echo ${base_file} | sed -e 's/ENIGMA3_mixed_se_//g' -e 's/wTHICK_//g' -e 's/wSA_//g' -e 's/wo_//g' -e 's/_20190429//g' -e 's/_20200522//g' -e 's/.txt.gz//g' )
#echo ${base_pheno}

# QC base file - on MAF (col4) and INFO, if available
# convert alleles to uppercase
zcat ${base_file} | awk 'NR==1 || ($4 > 0.01)  {print}' | \
 awk '{ print $1,toupper($2),toupper($3),$4,$5,$6,$7,$8,$9,$10,$11 }' > ${base_pheno}.txt


# Extract and rename header to match prsice script
#-snp snpid --chr chr --bp pos --A1 a_1 --A2 a_0 --stat beta --se se --pvalue P
sed -i -e 's/SNP/rsid/g' -e 's/CHR/chr/g' -e 's/BP/pos/g' -e 's/BETA1/beta/g' -e 's/A1/a_1/g' -e 's/A2/a_0/g'  -e 's/SE/se/g' -e 's/Pval/P/g' \
${base_pheno}.txt

# Ambiguous SNPs
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' ${base_pheno}.txt | \
    gzip > ${base_pheno}.nonamb.txt.gz

# Identify duplicated SNPs
gunzip -c ${base_pheno}.nonamb.txt.gz | awk '{ print $1}' | sort | uniq -d > ${base_pheno}.duplicated.snp
# Remove duplicated SNPs from base file
echo "Remove duplicated SNPs from base file and save clean file"
gunzip -c ${base_pheno}.nonamb.txt.gz | grep -vf ${base_pheno}.duplicated.snp | gzip - > ${working_dir}${base_pheno}.QC.gz

zcat ${base_pheno}.QC.gz | head -n 5
#rsid    chr     pos     a_1     a_0     EAF     beta    se      P
#rs2883059       3       49902160        T       C       0.5746  0.026   0.003   7.149e-25
#rs3796386       3       49899795        A       G       0.4235  -0.026  0.003   8.883e-25

# The ${base_pheno}.QC.gz base data are now ready for using in downstream analyses.

# clean intermediate files
rm ${base_pheno}.txt ${base_pheno}.nonamb.txt.gz ${base_pheno}.duplicated.snp

ls -lh ${working_dir}${base_pheno}.QC.gz

fi

#---------------------------------------------------------------------- 