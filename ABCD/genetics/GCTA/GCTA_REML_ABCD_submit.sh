module load R/4.0.3

# submit reml runs

# run
pheno_v=ABCDv3_DEAP
root=ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned_rm05_adj


# working directory
project_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/
scripts_dir=${project_dir}/scripts/genetics/GCTA/
working_dir=${project_dir}/data/working_data/genotyping/GCTA/
pheno_dir=${project_dir}/data/working_data/phenotypes/${pheno_v}/all/reading/destrieux/summary/tables/

cd ${scripts_dir}
#dos2unix *

mkdir -p sge_jobs
cd sge_jobs

# define phenotypes to check
cogn=$(echo nihtbx_reading_uncorrected nihtbx_picvocab_uncorrected nihtbx_cryst_uncorrected nihtbx_fluidcomp_uncorrected pea_wiscv_tss)
dos2unix ${pheno_dir}/ROImeasures2follow.txt
roisL=$(sed ':a;N;$!ba;s/\n/ /g' ${pheno_dir}/ROImeasures2follow.txt)
roisR=$(echo ${roisL} | sed 's/.lh/.rh/g')
rois2check=$(echo $roisL $roisR)

global=$(echo smri_vol_cort.desikan_total \
smri_area_cort.destrieux_total.rh smri_area_cort.destrieux_total.lh smri_area_cort.destrieux_total \
smri_thick_cort.destrieux_mean.lh smri_thick_cort.destrieux_mean.rh )

# all
phenos2check=$(echo $cogn \
 $global $rois2check
)
#-------------------------------------------------#
# univariate analysis (heritability)              #
#-------------------------------------------------#
# run regressed phenotypes (residuals from reml run in R)
for phenotype in ${phenos2check}
 do
 if [ ! -f ${working_dir}/pheno/residuals_${phenotype}.txt ]; then 
  Rscript ${scripts_dir}pheno_lm4reml.R ${root} ${phenotype} ${pheno_v}
 fi
done

for phenotype in ${phenos2check}
 do
 if [ ! -f ${working_dir}/reml/${phenotype}_residuals.hsq ]
 then
  qsub ${scripts_dir}GCTA_REML_ABCD_template.sh ${root} ${phenotype}_residuals ${working_dir}/pheno/residuals_${phenotype}.txt
 fi
done

for phenotype in ${rois2check}
 do
 if [ ! -f ${working_dir}/reml/${phenotype}_residuals_adjGlobal.hsq ]
 then
  f=$(ls ${working_dir}/pheno/residuals_adjGlobal*.*_${phenotype}.txt )
  qsub ${scripts_dir}GCTA_REML_ABCD_template.sh ${root} ${phenotype}_residuals_adjGlobal ${f}
 fi
done

#-------------------------------------------------

#-------------------------------------------------#
# bivariate analyses (genetic correlation)        #
#-------------------------------------------------#
## across cognitive measures
for phenotype1 in ${cogn}
 do
 for phenotype2 in ${cogn}
 do
  if [ ${phenotype1} != ${phenotype2} ]
  then
  if [ ! -f ${working_dir}/reml/${phenotype1}_residuals_${phenotype2}_residuals.log ]
   then
   if [ ! -f ${working_dir}/reml/${phenotype2}_residuals_${phenotype1}_residuals.log ] 
   then
   echo  ${phenotype1}_residuals ${phenotype2}_residuals
     qsub ${scripts_dir}GCTA_REML_bivar_ABCD_template.sh ${root} \
         ${phenotype1}_residuals ${phenotype2}_residuals \
         ${working_dir}/pheno/residuals_${phenotype1}.txt ${working_dir}/pheno/residuals_${phenotype2}.txt
  fi
  fi; fi
 done
done

# between global and regional measures
for phenotype1 in smri_area_cort.destrieux_total.lh smri_thick_cort.destrieux_mean.lh
 do
 m1=$(echo ${phenotype1} | awk -F "_" '{print $2}')
 for phenotype2 in ${roisL}
 do
 m2=$(echo ${phenotype2} | awk -F "_" '{print $2}')
  if [ "${m1}" == "${m2}" ]
  then
  if [ ${phenotype1} != ${phenotype2} ]
  then
  if [ ! -f ${working_dir}/reml/${phenotype1}_residuals_${phenotype2}_residuals.log ]
   then
   if [ ! -f ${working_dir}/reml/${phenotype2}_residuals_${phenotype1}_residuals.log ] 
   then
   echo  ${phenotype1}_residuals ${phenotype2}_residuals
     qsub ${scripts_dir}GCTA_REML_bivar_ABCD_template.sh ${root} \
         ${phenotype1}_residuals ${phenotype2}_residuals \
         ${working_dir}/pheno/residuals_${phenotype1}.txt ${working_dir}/pheno/residuals_${phenotype2}.txt
  fi; fi
  fi; fi
 done
done


## between left and right rois
for phenotype1 in ${roisL}
 do
  phenotype2=$(echo $phenotype1 | sed 's/lh/rh/')
  if [ ! -f ${working_dir}/reml/${phenotype1}_residuals_${phenotype2}_residuals.log ]
   then
   if [ ! -f ${working_dir}/reml/${phenotype2}_residuals_${phenotype1}_residuals.log ] 
   then
   echo  ${phenotype1}_residuals ${phenotype2}_residuals
     qsub ${scripts_dir}GCTA_REML_bivar_ABCD_template.sh ${root} \
         ${phenotype1}_residuals ${phenotype2}_residuals \
         ${working_dir}/pheno/residuals_${phenotype1}.txt ${working_dir}/pheno/residuals_${phenotype2}.txt
  fi; fi
 done

## between left and right rois, adjusting for global measures
for phenotype1 in ${roisL}
 do
  phenotype2=$(echo $phenotype1 | sed 's/lh/rh/')
  if [ ! -f ${working_dir}/reml/${phenotype1}_residuals_adjGlobal_${phenotype2}_residuals_adjGlobal.log ]
   then
   if [ ! -f ${working_dir}/reml/${phenotype2}_residuals_adjGlobal_${phenotype1}_residuals_adjGlobal.log ] 
   then
   echo  ${phenotype1}_residuals_adjGlobal ${phenotype2}_residuals_adjGlobal
    f1=$(ls ${working_dir}/pheno/residuals_adjGlobal*.*_${phenotype1}.txt )
    f2=$(ls ${working_dir}/pheno/residuals_adjGlobal*.*_${phenotype2}.txt )
     qsub ${scripts_dir}GCTA_REML_bivar_ABCD_template.sh ${root} \
         ${phenotype1}_residuals_adjGlobal ${phenotype2}_residuals_adjGlobal \
         ${f1} ${f2}
  fi; fi
 done


#--------------------------------------------
# organize output
#--------------------------------------------
# extract summary tables for h2 and rho
cd ${working_dir}/reml/
#mv *log logs/

########################################
## h2 ##						########
########################################

hsq_file=$(ls -1 *hsq | grep -e residuals | \
             grep -v -e _residuals_nihtbx_ -e _residuals_smri -e _residuals_adjGlobal_nihtbx_ -e _residuals_adjGlobal_smri_ -e _residuals_pea_wiscv_tss_ | sort | uniq ) # exclude rg runs (i.e. two phenos)
if [ $(echo ${hsq_file} | wc | awk '{print $2}') -ne 0 ]
 then
  paste <(grep "^n" -H ${hsq_file} ) <(grep "/" ${hsq_file} | awk '{print $2,$3}') <(grep "Pval" ${hsq_file} | awk '{print $2}') > hsq_h2_summary_${root}.table
fi

########################################
## rho ##						########
########################################
# extract genetic correlation values
# and with pvalue (diff to 0)
hsq_file=$(ls -1 *hsq | grep -e residuals | \
             grep -e _residuals_nihtbx_ -e _residuals_smri -e _residuals_adjGlobal_nihtbx_ -e _residuals_adjGlobal_smri_ -e _residuals_pea_wiscv_tss | sort | uniq)
if [ $(echo ${hsq_file} | wc | awk '{print $2}') -ne 0 ]
 then
 paste <(grep "^rG" -H ${hsq_file} ) <(grep "Pval" ${hsq_file} | awk '{print $2}')> hsq_rg_summary_${root}.table
fi
