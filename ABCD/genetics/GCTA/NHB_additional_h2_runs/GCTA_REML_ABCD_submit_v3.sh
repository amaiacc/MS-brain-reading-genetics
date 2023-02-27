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
rois2check=$(echo $roisL) # $roisR

global=$(echo smri_vol_cort.desikan_total \
smri_area_cort.destrieux_total.rh smri_area_cort.destrieux_total.lh smri_area_cort.destrieux_total \
smri_thick_cort.destrieux_mean.lh smri_thick_cort.destrieux_mean.rh )

# all
phenos2check=$(echo $cogn \
 $global $roisL
)
#-------------------------------------------------#
# univariate analysis (heritability)              #
#-------------------------------------------------#
# run regressed phenotypes (residuals from reml run in R)
for phenotype in ${phenos2check}
 do
 if [ ! -f ${working_dir}/pheno_v3/residuals_${phenotype}.txt ]; then 
  Rscript ${scripts_dir}pheno_lm4reml_v3.R ${root} ${phenotype} ${pheno_v}
 fi
done

mkdir -p ${working_dir}/reml

for phenotype in ${phenos2check}
 do
 if [ ! -f ${working_dir}/reml/${phenotype}_residuals.hsq ]
 then
  qsub ${scripts_dir}GCTA_REML_ABCD_template.sh ${root} ${phenotype}_residuals ${working_dir}/pheno_v3/residuals_${phenotype}.txt
 fi
done

for phenotype in ${rois2check}
 do
 if [ ! -f ${working_dir}/reml/${phenotype}_residuals_adjGlobal.hsq ]
 then
  f=$(ls ${working_dir}/pheno_v3/residuals_adjGlobal*.*_${phenotype}.txt )
  qsub ${scripts_dir}GCTA_REML_ABCD_template.sh ${root} ${phenotype}_residuals_adjGlobal ${f}
 fi
done

#-------------------------------------------------


#--------------------------------------------
# organize output
#--------------------------------------------
mv ${working_dir}/reml ${working_dir}/reml_v3
# extract summary tables for h2 and rho
cd ${working_dir}/reml_v3/
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
