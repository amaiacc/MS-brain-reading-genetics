#!/bin/bash
# Heritability and genetic correlation analyses of selected imaging derived phenotypes ( UKB (BIG40) )
## for regions selected by the ABCD brain-behaviour analysis, i.e. associated with reading
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/

#----------------------------------------------------------------------
base_dir=/export/home/acarrion/acarrion/projects/resources/datasets/
ukb_dir=${base_dir}GWAS_sumstats/downloaded_data/UKB/BIG40/
primary_dir=${base_dir}GWAS_sumstats/QC/UKB_BIG40/
working_dir=${base_dir}GWAS_sumstats/ldsc/UKB_BIG40/
#
gwas_dir=${base_dir}GWAS_sumstats/ldsc/
ieu_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/ieu-openGWAS/

#analysis_dir=${base_dir}ABCD/scripts/genetics/ldsc/
analysis_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/

scripts_ldsc_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
out_dir=${base_dir}ABCD/data/working_data/genotyping/ldsc/
#----------------------------------------------------------------------
mkdir -p ${working_dir} ${out_dir} ${analysis_dir}
mkdir -p ${out_dir}/logs/
## requirements
# Select and download regions to format :  ${ukb_dir}/scripts/select_phenos2download.R
# Run QC: $scripts_ldsc_dir/QC_sumstats/submit_publicGWAS_QC_UKB_BIG40.sh

#---------------------------------------------------------------------- 
# select brain measures (UKB BIG40) to include in analysis
pheno_file=${ukb_dir}IDPs.csv

global=$(grep Total ${ukb_dir}IDPs_summary.csv  | grep -e IntraCranial -e aparc-Desikan_ | \
 awk -F "," '{print $1}' | sed 's/"//g')
global=$(echo 0$global | sed 's/ / 0/g') # quick and dirty to get to file names with a 0 (i.e. 4 digit pheno codes)
rois=$(awk '{print $1}' ${ukb_dir}IDPs_ABCD_BrainBehaviour.txt | grep -v pheno)

roisCSAlh=$(grep area ${ukb_dir}IDPs_ABCD_BrainBehaviour.txt | grep lh | awk '{print $1}' | grep -v pheno)
roisCTlh=$(grep thick ${ukb_dir}IDPs_ABCD_BrainBehaviour.txt | grep lh | awk '{print $1}' | grep -v pheno)

# total left CSA
globalCSA=0648
# total left CT
globalCT=1020

#------------------------------------
# run
mkdir -p ${scripts_ldsc_dir}sge
cd ${scripts_ldsc_dir}sge
#-----------------------------------
# LDSC format and h2 already run: submit_ldsc_UKB_BIG40_readingROIs.sh
#------------------------------------
# Genetic correlations with cognitive traits - published GWASes from Howe et al. 2022 (sibGWAS paper)
suffix=Howe2022
# create variable containing a concatenation of all the phenotypes to be tested
vars2test=$(ls ${ieu_dir}*tsv.gz | sed 's|.tsv.gz|.tsv.gz-AND-|g')
vars2test=$(echo $vars2test | sed 's| ||g' | sed 's|-AND-$||g')
 

for base_pheno in ${global} ${roisCSAlh} ${roisCTlh}
 do
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
  if [ ! -f ${working_dir}/${base_pheno}_${suffix}_rg.log ] && [ ! -f ${working_dir}/logs/${base_pheno}_${suffix}_rg.log ]
  then
   echo Compute rg between ${base_pheno} and traits of interest using LDSC
   qsub -v phenos2test=${vars2test},p1=${base_pheno},suffix=${suffix} ${scripts_ldsc_dir}ldsc_rg_cognitive.sh
  fi
  fi
 done

 
#------------------------------------
# copy all relevant files to output directory
# rm ${out_dir}/logs/*
for r in ${global} ${rois} ${phensLR}
 do
 cp ${working_dir}/*${r}*log ${out_dir}/logs/
 done

#------------------------------------
# create summary tables
#------------------------------------
cd ${out_dir}/logs/
## h2 estimates
# columns to get:
##Total Observed scale h2: 0.0103 (0.0027)
##Lambda GC: 1.0466
##Mean Chi^2: 1.0469
##Intercept: 1.0093 (0.006)
##Ratio: 0.1984 (0.1285)
h2_file=$(ls -1 *h2.log )
paste <(echo 'File') <(echo 'h2 (se)') <(echo 'Lambda GC') <(echo 'Mean Chi^2') <(echo 'Intercept (se)') <(echo 'Ratio (se)')> ${out_dir}summary_h2_ldsc.table
paste <(ls -1 ${h2_file}) \
      <(grep "h2:" ${h2_file} | awk '{print $5,$6}' ) \
      <(grep "Lambda GC:" ${h2_file} | awk '{print $3}') \
      <(grep "Mean Chi^2" ${h2_file} | awk '{print $3}') \
      <(grep "Intercept" ${h2_file} | awk '{print $2,$3}') \
      <(grep "Ratio" ${h2_file} | awk '{print $2,$3}') >> ${out_dir}summary_h2_ldsc.table

## genetic correlations between left and right
grep -A11 p2 *${suffix}_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > ${out_dir}summary_${suffix}_rg_ldsc.table
#------------------------------------