#!/bin/bash
root=ABCD
batch=release_3.0_QCed
#--------------------------------
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/data/primary_data/ABCDgenomicsDataV30/genomics_sample03/imputed/

project_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/
geno_scripts=${project_dir}scripts/genetics/imputed/
working_dir=${project_dir}data/working_data/genotyping/imputed_genotypes/

# create working dir
mkdir -p ${working_dir}
mkdir -p ${geno_scripts}/sge
#----------------------------------------------------------------------
# QC imputed data - ABCD dataset
#----------------------------------------------------------------------
# 1- data has been imputed (TOPmed) by the ABCD team
# imputed data location is defined by the primary_dir

# 2- clean and combine imputed files, and create analysis ready clean files
cd ${geno_scripts}/sge
for i in {1..22} X # {1..15} {21..22} X
 do
 echo chr${i}
 qsub ${geno_scripts}imputation_output_filter.sh ${primary_dir} ${working_dir} ${root}_${batch} chr${i}
done

qsub ${geno_scripts}imputation_output_combine.sh ${working_dir} ${root}_${batch}
#----------------------------------------------------------------------
