#!/bin/sh
#$ -N prsice_cogn_noregress
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------
module load plink/6.15
#module load R/3.5.1
module load R/4.0.3
# required libraries in R (already available)
## "ggplot2","data.table","optparse","methods","tools","grDevices","RColorBrewer"
#----------------------------------------------------------------------
prsice_dir=/export/home/acarrion/acarrion/projects/resources/bin/PRSice_v2.2.12/
#----------------------------------------------------------------------
# get arguments
root=${1}
data=${2}
base_file=${3}
base_pheno=${4}
target_dir=${5}
base_stat=${6}
a1=${7}
a2=${8}
#----------------------------------------------------------------------
#add_dir=/resources/datasets/
cd /export/home/acarrion/acarrion/projects/${root}/data/working_data/genotyping/prsice/

if [ "$base_stat" == beta ] || [ "$base_stat" == Effect ]
 then
 stat="--beta"
 else
 stat=""
fi

#----------------------------------------------------------------------
# QC of target files
if [ ! -f ${data}.qc.bim ]
then
 plink --bfile ${target_dir}${data}.nodup.r2_rs_hg19 \
    --maf 0.05 \
    --mind 0.1 \
    --geno 0.1 \
    --hwe 1e-6 \
    --make-just-bim \
    --make-just-fam \
    --out ${data}.qc
fi
#----------------------------------------------------------------------
if [ ! -f ${base_pheno}.txt ]
then
 zcat ${base_file} > ${base_pheno}.txt
 sed -i -e 's/ /\t/g' ${base_pheno}.txt # to ensure that all cols are separated with tabs
fi
head ${base_pheno}.txt
#----------------------------------------------------------------------
if [ ! -f ${data}_${base_pheno}.valid ]
then
Rscript ${prsice_dir}PRSice.R --dir ${prsice_dir} \
    --prsice ${prsice_dir}PRSice_linux \
    --base ${base_pheno}.txt \
    --target ${target_dir}${data}.nodup.r2_rs_hg19 \
    --keep ${data}.qc.fam \
    --extract ${data}.qc.bim \
    --snp rsid --chr chr --bp pos --A1 ${a1} --A2 ${a2} --stat ${base_stat} --pvalue P ${stat} \
    --binary-target F \
    --bar-levels 0.00000005,0.00001,0.0001,0.001,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --print-snp \
    --no-regress \
    --out ${data}_${base_pheno}
fi


if [ -f ${data}_${base_pheno}.valid ] && [ ! -f ${data}_${base_pheno}.prsice ]
then
Rscript ${prsice_dir}PRSice.R --dir ${prsice_dir} \
    --prsice ${prsice_dir}PRSice_linux \
    --base ${base_pheno}.txt \
    --target ${target_dir}${data}.nodup.r2_rs_hg19 \
    --keep ${data}.qc.fam \
    --extract ${data}.qc.bim \
    --extract ${data}_${base_pheno}.valid \
    --snp rsid --chr chr --bp pos --A1 a_1 --A2 a_0 --stat ${base_stat} --pvalue P ${stat} \
    --binary-target F \
    --bar-levels 0.00000005,0.00001,0.0001,0.001,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --print-snp \
    --no-regress \
    --out ${data}_${base_pheno}
fi
#
#--lower=0.00000005 --upper=1 --interval=0.0005 \

rm ${base_pheno}.txt
