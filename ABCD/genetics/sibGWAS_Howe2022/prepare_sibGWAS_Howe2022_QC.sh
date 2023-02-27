#!/bin/bash
# PRS
#$ -N prepare_sibGWAS_base
#$ -o /export/home/acarrion/acarrion/projects/general_scripts/genotyping/QC_sumstats/sge/
#$ -e /export/home/acarrion/acarrion/projects/general_scripts/genotyping/QC_sumstats/sge/
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
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/ieu-openGWAS/ # mungedSumstats, original files are in vcf format
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/ieu-openGWAS/
mkdir -p ${working_dir}
#----------------------------------------------------------------------
# get arguments
base_file=${1}
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/
#---------------------------------------------------------------------- 

cd ${working_dir}
if [ ! -f ${base_pheno}.QC.txt.gz ] # if output does not exist
 then
# check structure of files, one file
#SNP     CHR     BP      A1      A2      END     FILTER  BETA    SE      LP      N       P       Z
#rs12238997      1       693731  A       G       693731  PASS    -0.008  0.0084  0.470057        61428   0.338799686659859       -0.956539601786072
#rs58276399      1       731718  T       C       731718  PASS    -0.0071 0.0083  0.410497        61428   0.388600182300083       -0.862158719605615
#rs141242758     1       734349  T       C       734349  PASS    -0.0074 0.0083  0.429107        61428   0.372299968951697       -0.892173453065721

# extract phenotype name from file
base_pheno=$(echo ${base_file} | sed -e 's/.tsv.gz//g' )

# QC info
## no MAF or info available in the original vcf file either (AF and SI tags are not available)
zgrep "#" ${primary_dir}../../downloaded_data/ieu-openGWAS/${base_pheno}.vcf.gz > ${base_pheno}_vcf_header.txt
zgrep -v "#" ${primary_dir}../../downloaded_data/ieu-openGWAS/${base_pheno}.vcf.gz | head  >> ${base_pheno}_vcf_header.txt


# Extract and rename header to match prsice script, and fix chromosome names

#-snp snpid --chr chr --bp pos --A1 a_1 --A2 a_2 --stat beta --se se --pvalue P
zcat ${primary_dir}${base_file} | awk '(NR==1){print "chr rsid pos a_1 a_2 beta se pval(-log10) P"}(NR>1){print $2,$1,$3,$4,$5,$8,$9,$10,$12}' > ${base_pheno}.txt

# Ambiguous SNPs
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' ${base_pheno}.txt | \
    gzip > ${base_pheno}.nonamb.txt.gz

# Identify duplicated SNPs
gunzip -c ${base_pheno}.nonamb.txt.gz | awk '{ print $2}' | sort | uniq -d > ${base_pheno}.duplicated.snp

# if duplicated snps file not empty
if [ -s ${base_pheno}.duplicated.snp ]
 then
 # Remove duplicated SNPs from base file
 echo "Remove duplicated SNPs from base file and save clean file"
 gunzip -c ${base_pheno}.nonamb.txt.gz | grep -vf ${base_pheno}.duplicated.snp | gzip - > ${base_pheno}.QC.txt.gz
else
 # copty file to QC.txt.gz
 cp ${base_pheno}.nonamb.txt.gz ${base_pheno}.QC.txt.gz
fi


zcat ${base_pheno}.QC.txt.gz | head -n 5
#chr rsid pos a_1 a_2 beta se pval(-log10) P
#1 rs12238997 693731 A G -0.008 0.0084 0.470057 0.338799686659859
#1 rs58276399 731718 T C -0.0071 0.0083 0.410497 0.388600182300083
#1 rs141242758 734349 T C -0.0074 0.0083 0.429107 0.372299968951697
#1 rs3094315 752566 G A 0.0052 0.0075 0.30883 0.491100074483958


ls -lh ${base_pheno}.QC.txt.gz

# The ${base_pheno}.QC.txt.gz base data are now ready for using in downstream analyses.

# clean intermediate files
rm ${base_pheno}.txt ${base_pheno}.nonamb.txt.gz ${base_pheno}.duplicated.snp


fi

#---------------------------------------------------------------------- 