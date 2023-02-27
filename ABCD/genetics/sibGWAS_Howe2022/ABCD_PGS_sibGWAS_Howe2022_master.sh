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
base_project=sibGWAS_Howe2022
#----------------------------------------------------------------------
# build parameters
#----------------------------------------------------------------------
add_dir=/resources/datasets/
analysis_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/scripts/genetics/${base_project}/
working_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/data/working_data/genotyping/prsice/
#----------------------------------------------------------------------
# RUN PGS calculation
bash ${analysis_dir}/ABCD_PGS_${base_project}_run.sh
# move log files to sge folders
cd ${working_dir}
mv *.e* *.o* ${analysis_dir}../PRS/sge
# reorgnaize output
cd ${working_dir}
mkdir -p ${base_project}
mv *ieu* ${base_project}/
#----------------------------------------------------------------------
# PARSE output
qsub ${analysis_dir}/ABCD_PGS_${base_project}_parse.sh
#----------------------------------------------------------------------
