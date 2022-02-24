#!/bin/bash
#---------------------------------------------------------------------- 
# Prepare public GWAS as base datasets for PRS or other sumstats-based downstream analyses
## since these are study-specific summary statistics, they have been QC'd one by one (given the different formats)
#---------------------------------------------------------------------- 
## following: https://choishingwan.github.io/PRS-Tutorial/base/


#----------------------------------------------------------------------
# Define variables
#----------------------------------------------------------------------
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/
mkdir -p ${working_dir}


#---------------------------------------------------------------------- 
## Educational Attainment 3, Lee et al. 2018
## https://www.thessgac.org/data
#---------------------------------------------------------------------- 
cd ${working_dir}
base_pheno=EA3_Lee2018
base_file=${primary_dir}Lee2018_EA3/GWAS_EA_excl23andMe.txt

# Check that file is not corrupted
md5sum ${base_file} > ${base_pheno}.md5.txt # to compare with original md5sum...
md5=$(awk '{print $1}' ${base_pheno}.md5.txt)
grep ${md5} ${primary_dir}Lee2018_EA3/*

# Check file structure
head ${base_file}
#MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval
#rs13090388      3       49391082        C       T       0.6905  -0.02852        0.00184 4.29e-54
#rs7630869       3       49522543        C       T       0.6922  -0.02848        0.00184 4.61e-54

# QC base file - on MAF (col6) and INFO, if available
awk 'NR==1 || ($6 > 0.01) && ($6 < 0.99) {print}' ${base_file} > ${base_pheno}.txt

# Extract and rename header to match script
#-snp snpid --chr chr --bp pos --A1 a_1 --A2 a_0 --stat beta --se se --pvalue P
sed -i -e 's/MarkerName/rsid/g' -e 's/CHR/chr/g' -e 's/POS/pos/g' -e 's/A1/a_1/g' -e 's/A2/a_0/g' -e 's/Beta/beta/g' -e 's/SE/se/g' -e 's/Pval/P/g' \
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
gunzip -c ${base_pheno}.nonamb.txt.gz | grep -vf ${base_pheno}.duplicated.snp | gzip - > ${base_pheno}.QC.gz

#zless ${base_pheno}.QC.gz
#rsid    chr     pos     a_1     a_0     EAF     beta    se      P
#rs13090388      3       49391082        C       T       0.6905  -0.02852        0.00184 4.29e-54
#rs7630869       3       49522543        C       T       0.6922  -0.02848        0.00184 4.61e-54

# The ${base_pheno}.QC.gz base data are now ready for using in downstream analyses.

# clean intermediate files
rm ${base_pheno}.txt ${base_pheno}.nonamb.txt.gz ${base_pheno}.duplicated.snp

#---------------------------------------------------------------------- 
# Cognitive performance 
## Educational Attainment 3, Lee et al. 2018
## https://www.thessgac.org/data
#---------------------------------------------------------------------- 
cd ${working_dir}
base_pheno=CP_Lee2018
base_file=${primary_dir}Lee2018_EA3/GWAS_CP_all.txt

# Check that file is not corrupted
md5sum ${base_file} > ${base_pheno}.md5.txt # to compare with original md5sum...
md5=$(awk '{print $1}' ${base_pheno}.md5.txt)
grep ${md5} ${primary_dir}Lee2018_EA3/*

# Check file structure
head ${base_file}
#MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval
#rs2352974       3       49890613        T       C       0.5     -0.0319 0.00285 5.19e-29
#rs2271960       3       49878078        C       T       0.5068  -0.03138        0.00285 3.19e-28

# QC base file - on MAF (col6) and INFO, if available
awk 'NR==1 || ($6 > 0.01) && ($6 < 0.99) {print}' ${base_file} > ${base_pheno}.txt

# Extract and rename header to match script
#-snp snpid --chr chr --bp pos --A1 a_1 --A2 a_0 --stat beta --se se --pvalue P
sed -i -e 's/MarkerName/rsid/g' -e 's/CHR/chr/g' -e 's/POS/pos/g' -e 's/A1/a_1/g' -e 's/A2/a_0/g' -e 's/Beta/beta/g' -e 's/SE/se/g' -e 's/Pval/P/g' \
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
gunzip -c ${base_pheno}.nonamb.txt.gz | grep -vf ${base_pheno}.duplicated.snp | gzip - > ${base_pheno}.QC.gz

#zless ${base_pheno}.QC.gz
#rsid    chr     pos     a_1     a_0     EAF     beta    se      P
#rs13090388      3       49391082        C       T       0.6905  -0.02852        0.00184 4.29e-54
#rs7630869       3       49522543        C       T       0.6922  -0.02848        0.00184 4.61e-54

# The ${base_pheno}.QC.gz base data are now ready for using in downstream analyses.

# clean intermediate files
rm ${base_pheno}.txt ${base_pheno}.nonamb.txt.gz ${base_pheno}.duplicated.snp
#---------------------------------------------------------------------- 

#---------------------------------------------------------------------- 
# Developmental dyslexia
# Gialluisi et al. 2020
#---------------------------------------------------------------------- 
cd ${working_dir}
base_pheno=DDyslexia_Gialluisi2020
base_file=${primary_dir}/Gialluisi2020_DD/DD.full.GWASsummaryStats/MetaAnalysis_cc1.tbl

# Check that file is not corrupted
md5sum ${base_file} > ${base_pheno}.md5.txt # to compare with original md5sum...
md5=$(awk '{print $1}' ${base_pheno}.md5.txt)

#head ${base_file}
#MarkerName      Allele1 Allele2 Weight  Zscore  P-value Direction       HetISq  HetChiSq        HetDf   HetPVal
#rs2326918       a       g       8546.00 0.193   0.8468  -+--+++ 0.0     0.617   6       0.9961
#rs977590        a       g       8546.00 0.989   0.3228  ---+-++ 0.0     4.070   6       0.6672

# QC: Keep only variants that are not heterogeneous, i.e. HetPval>0.05
awk 'NR==1 || ($11 > 0.05 ) {print}' ${base_file} > ${base_pheno}.txt

# Extract and rename header to match script
#-snp snpid --chr chr --bp pos --A1 a_1 --A2 a_0 --stat beta --se se --pvalue P
sed -i -e 's/MarkerName/rsid/g' -e 's/CHR/chr/g' -e 's/BP/pos/g' -e 's/Allele1/a_1/g' -e 's/Allele2/a_0/g'  -e 's/Weight/N/g' -e 's/P-value/P/g' \
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
gunzip -c ${base_pheno}.nonamb.txt.gz | grep -vf ${base_pheno}.duplicated.snp | gzip - > ${base_pheno}.QC.gz
# Remove duplicated SNPs from base file
if [ -s ${base_pheno}.duplicated.snp ]
 then
  gunzip -c ${base_pheno}.nonamb.txt.gz | grep -vf ${base_pheno}.duplicated.snp | gzip - > ${base_pheno}.QC.gz
 else
  cp ${base_pheno}.nonamb.txt.gz ${base_pheno}.QC.gz
fi

# Check file structure
#zless ${base_pheno}.QC.gz
#rsid    a_1     a_0     N       Zscore  P       Direction       HetISq  HetChiSq        HetDf   HetPVal
#rs977590        a       g       8546.00 0.989   0.3228  ---+-++ 0.0     4.070   6       0.6672
#rs7929618       c       g       8546.00 1.058   0.2899  ++-+-++ 8.9     6.587   6       0.3607


# clean intermediate files
rm ${base_pheno}.txt ${base_pheno}.nonamb.txt.gz ${base_pheno}.duplicated.snp
