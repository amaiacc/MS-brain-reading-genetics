#!/bin/bash
# PRS
#$ -N PRS_run_dataset
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------
# Run PRS analysis on the ABCD dataset, using GenLang GWAS sumstats for reading as base datasets
#----------------------------------------------------------------------
bin_dir=/export/home/acarrion/acarrion/projects/resources/bin/
prsice_dir=${bin_dir}PRSice_v2.2.12/
PATH=${prisce_dir}:${PATH}
base_dir=/export/home/acarrion/acarrion/projects/coeduca/data/working_data/GenLang/projects/GWASMA_freeze_v1/data/QC/
#----------------------------------------------------------------------
# Load modules, define directories, root files, and PATH
#----------------------------------------------------------------------
module load plink/6.15
#module load R/3.5.1
module load R/4.0.3
#----------------------------------------------------------------------
# define variables
root=ABCD
batch=release_3.0_QCed
data=${root}_${batch}
add_dir=/resources/datasets/

#----------------------------------------------------------------------
# build parameters
#----------------------------------------------------------------------
add_dir=/resources/datasets/
analysis_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/scripts/genetics/PRS/
project_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/data/working_data/genotyping/
target_dir=${project_dir}/imputation/${data}/MichiganImputationServer/imputed_genotypes/clean/
if [ ! -d ${target_dir} ]
 then
 target_dir=${project_dir}imputed_genotypes/clean/
fi
working_dir=${project_dir}/prsice/
mkdir -p ${working_dir}
mkdir -p ${working_dir}/readingGenLang/

mkdir -p ${analysis_dir}/sge
#----------------------------------------------------------------------
cd ${analysis_dir}/sge

for base_pheno in WR_RT_EUR_combined_STERR_GCOFF_noABCD
 do
  # define stat beta or OR depending on the base_pheno name (whether it contains a disorder name or not)
  if echo ${base_pheno} | grep -q -e ADHD -e ASD -e SCZ 
  then 
   stat=OR
  elif echo ${base_pheno} | grep -q -e WR
  then 
   stat=Effect
  else
   stat=beta
  fi
  echo ${stat}
  # define base file
  base_file=${base_dir}/${base_pheno}.QC.ChrPos.gz
  if [ ! -f ${working_dir}${data}_${base_pheno}.prsice ]
  then
  if [ ! -f ${working_dir}/readingGenLang/${data}_${base_pheno}.prsice ]
  then
   echo run ${base_pheno}
   qsub ${analysis_dir}PRSice_run.sh ${add_dir}${root} ${data} ${base_file} ${base_pheno} ${target_dir} ${stat} a_1 a_0
  fi; fi
done
#----------------------------------------------------------------------
# reorganize output
mkdir -p ${working_dir}/
cd ${working_dir}
mv *WR* readingGenLang/