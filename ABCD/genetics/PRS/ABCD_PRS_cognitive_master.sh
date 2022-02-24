#----------------------------------------------------------------------
# Run PRS analysis on the ABCD dataset, using public GWAS sumstats for cognitive measures as base datasets
#----------------------------------------------------------------------
bin_dir=/export/home/acarrion/acarrion/projects/resources/bin/
prsice_dir=${bin_dir}PRSice_v2.2.12/
PATH=${prisce_dir}:${PATH}
#----------------------------------------------------------------------
# Load modules, define directories, root files, and PATH
#----------------------------------------------------------------------
# define variables
root=ABCD
batch=release_3.0_QCed
data=${root}_${batch}
base_project=cognitive
#----------------------------------------------------------------------
# build parameters
#----------------------------------------------------------------------
add_dir=/resources/datasets/
analysis_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/scripts/genetics/PRS/
#----------------------------------------------------------------------
cd ${analysis_dir}/sge
#bash ../ABCD_PRS_cognitive_run.sh

qsub ../ABCD_PRS_cognitive_parse.sh