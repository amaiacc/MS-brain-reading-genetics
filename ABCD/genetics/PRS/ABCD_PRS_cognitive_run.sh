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
base_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/
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
mkdir -p ${working_dir}/cognitive

mkdir -p ${analysis_dir}/sge
#----------------------------------------------------------------------
cd ${analysis_dir}/sge
for base_pheno in CP_Lee2018 EA3_Lee2018 
 do
  # define stat beta or OR depending on the base_pheno name (whether it contains a disorder name or not)
  if echo ${base_pheno} | grep -q -e ADHD -e ASD -e SCZ 
  then 
   stat=OR
  else
   stat=beta
  fi
  echo ${stat}
  # define base file
  base_file=${base_dir}/${base_pheno}.QC.gz
  if [ ! -f ${working_dir}${data}_${base_pheno}.prsice ] && [ ! -f ${working_dir}/enigma/${data}_${base_pheno}.prsice ]
  then
   qsub ${analysis_dir}PRSice_run.sh ${add_dir}${root} ${data} ${base_file} ${base_pheno} ${target_dir} ${stat} a_1 a_0
  fi
done
#----------------------------------------------------------------------
# reorgnaize output
mkdir -p cognitive
mv *CP* *EA* cognitive/