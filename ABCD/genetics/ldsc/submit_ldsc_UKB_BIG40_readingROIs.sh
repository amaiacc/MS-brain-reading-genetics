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
genlang_dir=/export/home/acarrion/acarrion/projects/coeduca/data/working_data/GenLang/projects/GWASMA_freeze_v1/data/ldsc/

scripts_dir=${base_dir}ABCD/scripts/genetics/ldsc/
scripts_ldsc_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
out_dir=${base_dir}ABCD/data/working_data/genotyping/ldsc/
#----------------------------------------------------------------------
mkdir -p ${working_dir} ${out_dir} ${scripts_dir}
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

#------------------------------------
# run
mkdir -p ${scripts_ldsc_dir}sge
cd ${scripts_ldsc_dir}sge
#-----------------------------------
# LDSC format: munge sumstats 
for base_pheno in ${global} ${rois}
 do
  base_file=${base_pheno}.QC.txt.gz
  if [ ! -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
   echo Submit LDSC formatting ${base_pheno}
   qsub ${scripts_ldsc_dir}ldsc_format_UKB_BIG40.sh ${base_pheno} ${pheno_file}
  fi
 done
# heritability
for base_pheno in ${global} ${rois}
 do
  base_file=${base_pheno}.QC.txt.gz
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
  if [ ! -f ${working_dir}/${base_pheno}_h2.log ]
  then
   echo Compute h2 using LDSC ${base_file}
   qsub ${scripts_ldsc_dir}ldsc_h2.sh ${base_pheno}
  fi
  fi
 done

# Genetic correlation between global and regional measures
## CSA and CSA measures
roisCSAlh=$(grep area ${ukb_dir}IDPs_ABCD_BrainBehaviour.txt | grep lh | awk '{print $1}' | grep -v pheno)
roisCTlh=$(grep thick ${ukb_dir}IDPs_ABCD_BrainBehaviour.txt | grep lh | awk '{print $1}' | grep -v pheno)

# total left CSA
globalCSA=0648

for base_pheno in ${roisCSAlh}
 do
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
  if [ ${base_pheno} != ${globalCSA} ]
  then
  if [ ! -f ${working_dir}/${globalCSA}_${base_pheno}_rg.log ]
  then
   echo Compute rg between ${globalCSA} and ${base_pheno} using LDSC 
   qsub ${scripts_ldsc_dir}ldsc_rg.sh ${globalCSA} ${base_pheno}
  fi
  fi; fi
 done

# total left CT
globalCT=1020

for base_pheno in ${roisCTlh}
 do
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
  if [ ${base_pheno} != ${globalCT} ]
  then
  if [ ! -f ${working_dir}/${globalCT}_${base_pheno}_rg.log ]
  then
   echo Compute rg between ${globalCT} and ${base_pheno} using LDSC 
   qsub ${scripts_ldsc_dir}ldsc_rg.sh ${globalCT} ${base_pheno}
  fi
  fi; fi
 done


#------------------------------------
# Genetic correlations with cognitive traits - published GWASes
for base_pheno in ${global} ${roisCSAlh} ${roisCTlh}
 do
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
  if [ ! -f ${working_dir}/${base_pheno}_cognitive_rg.log ]
  then
   echo Compute rg between ${base_pheno} and cognitive traits using LDSC 
   qsub ${scripts_ldsc_dir}ldsc_rg_cognitive.sh ${base_pheno} ${gwas_dir}
  fi
  fi
 done

#------------------------------------
# Run genetic correlations for L and R
lat_file=${ukb_dir}/IDPs_lat_summary.txt
dos2unix ${lat_file}

for r in ${global} ${rois}
 do
 grep -e ${r} ${lat_file} >> tmp
 done

sort tmp | uniq > ${out_dir}IDPs_lat_summary_ABCD_BrainBehaviour.txt
rm tmp

phensLR=$(awk '{print $3}' ${out_dir}IDPs_lat_summary_ABCD_BrainBehaviour.txt | uniq | grep -v -e "NA$")
for name in ${phensLR}
do
if [ ! -f ${working_dir}${name}_LR_rg.log ]
 then
 
 name2=$(echo ${name}"\t" | sed 's/+/and/g')
 p1=$(sed 's/+/and/g' ${lat_file} | grep -P ${name2} | awk '{print $5}') # use -P to grep including tab, otherwise it may give multiple matches (e.g. Crus_I and Crus_II)
 p2=$(sed 's/+/and/g' ${lat_file} | grep -P ${name2} | awk '{print $6}')
 if [ -f ${working_dir}/${p1}.sumstats.gz ] && [ -f ${working_dir}/${p2}.sumstats.gz ] 
 then
  echo 'Run LDSC to calculate genetic correlation between left and right for ' ${name}': '${p1}' and '${p2}
  qsub ${scripts_ldsc_dir}ldsc_rg_LR.sh ${name} ${p1} ${p2}
 fi
fi
done

#------------------------------------
# copy all relevant files to output directory
rm ${out_dir}/logs/*
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
grep -A11 p2 *LR_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > ${out_dir}summary_LRrg_ldsc.table
grep -A11 p2 *rg.log | grep -v -e LR -e cognitive -e readingGenLang | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > ${out_dir}summary_rg_ldsc.table
grep -A11 p2 *cognitive_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > ${out_dir}summary_cognitive_rg_ldsc.table
#------------------------------------