#!/bin/bash
# Prepare public GWAS sumstats from ENIGMA3 to be used as base datasets for PRS analyses
#---------------------------------------------------------------------- 
# From ENIGMA3 study, Grasby et al. 2020
## http://enigma.ini.usc.edu/downloads/
#---------------------------------------------------------------------- 

working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/ENIGMA3/
analysis_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/QC_sumstats/
mkdir -p ${working_dir}

# Select regions to format:
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/ENIGMA3/enigma.ini.usc.edu/downloads/ENIGMA3_withGlobal/

# enigma global measures, to have as reference
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/ENIGMA3/enigma.ini.usc.edu/downloads/ENIGMA3_Global/
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/ENIGMA3/enigma.ini.usc.edu/downloads/ENIGMA3_Global_noGC/
cd ${primary_dir}
files2format=$(ls *SurfArea* *Thickness*)

cd ${analysis_dir}
for base_file in ${files2format}
 do
  # extract phenotype name from file
  base_pheno=$(echo ${base_file} | sed -e 's/ENIGMA3_mixed_se_//g' -e 's/wTHICK_//g' -e 's/wSA_//g' -e 's/wo_//g' -e 's/_20190429//g' -e 's/_20200522//g' -e 's/.txt.gz//g' )
  if [ ! -f ${working_dir}${base_pheno}.QC.gz ]
  then
  echo $base_pheno
   qsub prepare_publicGWAS_QC_ENIGMA3.sh ${primary_dir}${base_file} ${base_pheno}
  fi
 done
