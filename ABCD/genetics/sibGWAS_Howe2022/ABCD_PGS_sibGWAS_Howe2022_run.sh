#!/bin/bash
# PRS
#$ -N PRS_run_dataset
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------
# Run PRS analysis on the coeduca dataset, using public GWAS sumstats for cognitive measures as base datasets
#----------------------------------------------------------------------
bin_dir=/export/home/acarrion/acarrion/projects/resources/bin/
prsice_dir=${bin_dir}PRSice_v2.2.12/
PATH=${prisce_dir}:${PATH}
base_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/ieu-openGWAS/
#----------------------------------------------------------------------
# Load modules, define directories, root files, and PATH
#----------------------------------------------------------------------
module load plink/6.15
module load R/3.5.1
#----------------------------------------------------------------------
# define variables
root=ABCD
batch=release_3.0_QCed
data=${root}_${batch}
base_project=sibGWAS_Howe2022
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
mkdir -p ${working_dir}/${base_project}
mkdir -p ${analysis_dir}/sge
#----------------------------------------------------------------------
cd ${base_dir}

for base_file in $(ls *QC.txt.gz )
 do
  # define stat beta or OR depending on the base_pheno name (whether it contains a disorder name or not)
  if echo ${base_pheno} | grep -q -e ADHD -e ASD -e SCZ 
  then 
   stat=OR
  else
   stat=beta
  fi
  echo ${stat}
  # extract phenotype name from file
  base_pheno=$(echo ${base_file} | sed -e 's/.QC.txt.gz//g' )
  if [ ! -f ${working_dir}${data}_${base_pheno}.prsice ] && [ ! -f ${working_dir}/${base_project}/${data}_${base_pheno}.prsice ]
  then
   qsub ${analysis_dir}PRSice_run.sh ${add_dir}${root} ${data} ${base_dir}/${base_file} ${base_pheno} ${target_dir} ${stat} a_1 a_2
  fi
done
#----------------------------------------------------------------------
