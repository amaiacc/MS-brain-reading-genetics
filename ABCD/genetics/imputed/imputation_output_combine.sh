#!/bin/bash
# Post-imputation QC: combined QCd files
#$ -N combine_imputed
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------

module load plink/6.15

 
#----------------------------------------------------------------------
# Filter imputed genotype vcf files
#----------------------------------------------------------------------
# Once imputation has been ran, we need to filter the SNPs for r2 > 0.8 and MAF > .01. 

#----------------------------------------------------------------------
# Define variables
#----------------------------------------------------------------------
# required parameters (hard-coded or to be taken as arguments)
# define them as command line arguments:
resource_dir=/export/home/acarrion/acarrion/projects/resources/reference_data/
#project_dir=/export/home/acarrion/acarrion/projects/
working_dir=${1}
data=${2}
#----------------------------------------------------------------------

cd ${working_dir}/tmp/

if [ ! -f ${working_dir}/clean/${data}.nodup.r2.bim ]
then
# merge all chromosomes into one bfile:
ls ${working_dir}/tmp/*.nodup.r2.bim | grep -v "chr1\\." | sed 's/.bim//g'  > chr_files.list
plink --bfile chr1.nodup.r2 --merge-list chr_files.list --make-bed --out ${working_dir}/clean/${data}.nodup.r2

else 
 echo 'Output files ('${data}'.nodup.r2) exist already, so skip.'
fi


if [ ! -f ${working_dir}/clean/${data}.nodup.r2_rs.bim ] 
then
 echo 'Convert SNP names to rs (when available in HRC)'
 # convert SNP names from chr:pos:A1:A2 to rsIDs
 awk 'NR>1 {id=$1":"$2":"$5":"$6; print id,$4}' ${resource_dir}HRC.r1-1.hg38.bed  > HRC_rs_ids.txt # only rs SNPs
 awk '{print $1}' HRC_rs_ids.txt | sort | uniq -c | awk '{ if ($1 >= 2) print $0 }' > HRC_rs_duplicatedpos.txt

 grep -v $(awk '{print $2}' HRC_rs_duplicatedpos.txt | sed 's/^/-e /g' )  HRC_rs_ids.txt > HRC_rs_NOdup_pos.txt

 awk '{print $2}' HRC_rs_NOdup_pos.txt | uniq -c |  awk '{ if ($1 < 2) print $0 }' | awk '{print $2}' > HRC_rs_NOdup_ids.txt # unique rsids, after removing multiallelic SNPs

 # update map, adding rsid's to SNPs in HRC
 plink --bfile ${working_dir}/clean/${data}.nodup.r2 --update-name HRC_rs_NOdup_pos.txt --make-bed --out ${working_dir}/clean/${data}.nodup.r2_rs
 
 ## update map, only for rsids
 plink --bfile ${working_dir}/clean/${data}.nodup.r2_rs --extract HRC_rs_NOdup_ids.txt --make-bed --out ${working_dir}/clean/${data}.nodup.r2_rs_tmp
 
 plink --bfile ${working_dir}/clean/${data}.nodup.r2_rs_tmp --update-map ${resource_dir}HRC_rsid_pos_hg19_unique.txt --make-bed --out ${working_dir}/clean/${data}.nodup.r2_rs_hg19

 # clean all intermediate files
 rm ${working_dir}/clean/${data}.nodup.r2_rs_tmp*
 cd ${working_dir}/tmp/ 
 #rm chr*bim chr*fam chr*bed chr*nosex *duplicated *snplist
else 
 echo 'Output files ('${data}'.nodup.r2_rs) exist already, so skip.'

fi
