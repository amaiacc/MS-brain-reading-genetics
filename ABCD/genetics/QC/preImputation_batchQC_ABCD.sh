#!/bin/bash
# Pre-imputation QC - ABCD dataset
#$ -N preimputationQC_batch
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Define variables
#----------------------------------------------------------------------
# required parameters (hard-coded or to be taken as arguments)
# define them as command line arguments:
project_dir=${1}
root=${2}
batch=${3}
input_dir=${4}
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Load modules, define directories, root files, and PATH
#----------------------------------------------------------------------
module load plink/6.15
module load R/4.0.3
module load gnuplot/5.8.2
module load gcc/7.3.0
#module load eigensoft/7.2.1
bin_dir=/export/home/acarrion/acarrion/projects/resources/bin/
PATH=${bin_dir}GenGen-1.0.1:${PATH}
PATH=${bin_dir}EIG-6.1.4/bin:${PATH}

# Reference samples (Population Stratification)
# - Hapmap3 (following ENIGMA2 protocol) (plink 1)
# - 1KG (following Coleman 2014) (plink v1.9)
# donwload as in:
resource_dir=/export/home/acarrion/acarrion/projects/resources/reference_data/
ref_HM=${resource_dir}HapMap/phase2/
ref_1KG=${resource_dir}1KG/phase1/
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Published protocols (/acarrion/refs/) and associated github repositories (/acarrion/projects/resources/github/)
#----------------------------------------------------------------------
# Quality control, imputation and analysis of genome-wide genotyping data from the Illumina HumanCoreExome microarray (Coleman et al. 2015)
gwas_scripts=/export/home/acarrion/acarrion/projects/resources/github/gwas_scripts/
# ENIGMA2 protocol for imputation
## /export/home/acarrion/acarrion/github/ENIGMA/Genetics/ENIGMA2/Imputation/
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# build parameters
geno_scripts=${project_dir}scripts/genetics/QC/
working_dir=${project_dir}data/working_data/genotyping/preImputation_QC/
#input_dir=${project_dir}data/working_data/genotyping/
#
data=${root}_${batch}
input_file=${input_dir}/${data}
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Preprocessing
#---------------------------------------------------------------------- 
# Create working directory if it doesn't exist, and go there
mkdir -p ${working_dir}/clean/
cd ${working_dir}

# only run if final files do not exist

if [ ! -f ${data}.IBD_PopStr_cleaned.bed ] && [ ! -f ./clean/${data}.IBD_PopStr_cleaned.bed ] 
then
echo Run

#----------------------------------------------------------------------
# Initial format conversions
#----------------------------------------------------------------------
if [ ! -f ${data}.updated_genders.bed ]
then
echo "----------------------------------------------------------------------"
echo "Convert files to plink format, update allele coding and gender information"
echo "----------------------------------------------------------------------"

# Update gender
echo "Update genders into ped files"

# famID is missing, so get it from the fam file (after double checking that iID is unique)
awk '{print $2}' ${input_file}.fam | sort | uniq -c | grep -v "1 " # individual IDs are unique

dos2unix ${project_dir}/data/working_data/phenotypes/*.txt
awk 'FNR==NR{a[$1]=$0;next}{print $1,a[$2]}' <(sort -k 1 ${project_dir}/data/working_data/phenotypes/ABCD_release2.0.1_sex.txt ) <(sort -k 2 ${input_file}.fam) | grep -e "F$" -e "M$" | \
 sed 's/F$/2/g' |  sed 's/M$/1/g' > ${data}.genders

plink --bfile ${input_file} --update-sex ${data}.genders --make-bed --out ${data}.updated_genders

echo "----------------------------------------------------------------------"
fi
#---------------------------------------------------------------------- 

#----------------------------------------------------------------------
# SNP QC
#----------------------------------------------------------------------
if [ ! -f ${data}.hw_dropped.bed ]
then
echo "----------------------------------------------------------------------"
echo "SNP QC: filter based on genotyping rate, MAF, HWE and missingness"
echo "----------------------------------------------------------------------"

echo "Filter SNPs"
plink --bfile ${data}'.updated_genders' --freq --out ${data}'_allSNPs'
awk '$5==0 { print $0 }' ${data}'_allSNPs.frq' > ${data}'_maf0.frq'

# missingness rates for SNPs and samples
plink --bfile ${data}'.updated_genders' --missing -out ${data}'_raw'

# remove samples with genotyping rate < 0.90 (initial QC to get SNP missingness)
plink --bfile ${data}'.updated_genders' --mind 0.1 -make-bed --out ${data}'_filtered1'
wc ${data}'_filtered1'.fam # 

# remove snps
### genotyping rate < 0.99
plink --bfile ${data}'_filtered1' --geno 0.01 --set-hh-missing --make-bed --out ${data}'_commonNrareSNPs'
plink --bfile  ${data}'_commonNrareSNPs' --missing --out ${data}'_commonNrareSNPs' # get missingness
### maf < 0.01
plink --bfile ${data}'_commonNrareSNPs' --geno 0.01 --maf 0.01 --set-hh-missing --make-bed --out ${data}'_commonSNPs'
plink --bfile  ${data}'_commonSNPs' --missing --out ${data}'_commonSNPs'

# remove samples that have a genotyping rate < 95% 
# check number of removed samples in logfile
echo "Filter samples"
plink --bfile ${data}'_commonSNPs' --mind 0.05 -make-bed --out ${data}'_filtered' # no samples removed
# Generate files for individual call rates and variant call rates.
echo "Check missingness after filtering"
plink --bfile ${data}'_filtered' --missing --out ${data}'_filtered'
#Examine the lowest call rates for variants:
## CHR                  SNP   N_MISS   N_GENO   F_MISS
sort -k 5 -gr ${data}'_filtered.lmiss' | grep -v '^  23' | head
#Examine the lowest call rates for individuals:
## FID     IID MISS_PHENO   N_MISS   N_GENO   F_MISS
sort -k 6 -gr ${data}'_filtered.imiss' | head

# keep copy of genotypes including rare, for imputation using HRC as a reference
plink --bfile ${data}'_commonNrareSNPs' --mind 0.05 -make-bed --out ${data}'_filtered2' 

wc -l ${data}*SNPs.bim
wc -l ${data}*filtered*.bim

# HWE: Assess SNPs for deviation from Hardy-Weinberg Equilibrium
echo "Assess SNPs for deviation from Hardy-Weinberg Equilibrium"
# --hardy calculates HWE test p-values:
plink --bfile ${data}'_filtered' --hardy --out ${data}'_filtered'
# --hwe removes deviant SNPs past a given threshold, 1x10^-6 below:
plink --bfile ${data}'_filtered' --hwe 0.000001 --make-bed --out  ${data}'.hw_dropped'

# filter on --hwe for hrc copy as well
plink --bfile ${data}'_filtered2' --hwe 0.000001 --make-bed --out  ${data}'.hw_dropped2'

# number of SNPs after applying HWE filterwc *.hw_dropped*.bim
wc -l ${data}*.hw_dropped*.bim

echo "----------------------------------------------------------------------"
fi
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Sample QC
#----------------------------------------------------------------------
echo "----------------------------------------------------------------------"
echo "Sample QC - using common SNPs"
echo "----------------------------------------------------------------------"

if [ ! -f ${data}_pruned_LD2.bed ]
then
echo "LD pruning - select independent SNPs for sex check, IBD and PCAs"
echo "----------------------------------------------------------------------"
# select independent SNPs
# prune LD using a window of 1500 variants and a shift of 150 variants between windows, with an r2 cut-off of 0.2:
plink --bfile ${data}'.hw_dropped' --indep-pairwise 1500 150 0.2 --make-bed --out ${data}'_pruned_LD' 
# exclude high-LD and non-autosomal regions from the pruned file 
# justification: https://sites.google.com/site/mikeweale/
# get file from https://github.com/JoniColeman/gwas_scripts/blob/master/highLDregions4bim_b37.awk
awk -f ${gwas_scripts}highLDregions4bim_b37.awk ${data}'_pruned_LD.bim' > highLDexcludes
awk '($1 < 1) || ($1 > 22) {print $2}' ${data}'_pruned_LD.bim' > autosomeexcludes
cat highLDexcludes autosomeexcludes > highLD_and_autosomal_excludes 
plink --bfile ${data}'_pruned_LD' --exclude highLD_and_autosomal_excludes --make-bed --out ${data}'_pruned_LD2' 
fi
#----------------------------------------------------------------------

if [ ! -f ${data}_sex_check_chrX.sexcheck ]
then
# Sex check
echo "Sex check"
echo "----------------------------------------------------------------------"
#Ensure there is a separate XY region for the pseudoautosomal region on X
# plink --bfile ${data}'_pruned_LD' --split-x b37 --make-bed --out ${data}'_pruned_split'
# Error: --split-x cannot be used when the dataset already contains an XY region.
plink --bfile ${data}'_pruned_LD' --check-sex ycount 0.2 0.8 0 1 --out ${data}'_sex_check' # considering both chrX and chrY
Rscript ${geno_scripts}sexcheck_plot.R ${data}'_sex_check'.sexcheck
# ignoring chrY, because there are many nonmale chrY calls - from GenomeStudio
plink --bfile ${data}'_pruned_LD' --check-sex 0.25 0.8 --out ${data}'_sex_check_chrX' # considering chrX, all OK.
Rscript ${geno_scripts}sexcheck_plot.R ${data}'_sex_check_chrX'.sexcheck

# check output
grep PROBLEM *chrX.sexcheck
# save samples to exclude given sex mismatches
grep PROBLEM *chrX.sexcheck | awk '($3 != "0") {print $1,$2}' > ${data}.sexcheck_mismatch_samples.txt

# update sex for samples where it was unkown
grep PROBLEM *chrX.sexcheck | awk '($3 == "0") {print $1,$2,$4}' *chrX.sexcheck > ${data}.sex2update

# remove samples with sex mismatches and update sex
plink --bfile ${data}'_pruned_LD2' --remove ${data}'.sexcheck_mismatch_samples.txt' --update-sex ${data}.sex2update --make-bed --out ${data}'_sexchecked'

fi

#----------------------------------------------------------------------

if [ ! -f ${data}'_pruned.no_close_relatives'.bed ]
then
echo "IBD test"
echo "----------------------------------------------------------------------"
plink --bfile ${data}'_sexchecked' --genome --out ${data}'_pruned.sex_ibd' # 

# Remove one sample from each pair with pi-hat (% IBD) above threshold (0.1 below, 0.1875 by JoniColeman pipeline):
# flag pair of outliers if they have pi_hat> 0.05 (?) - too stringent
awk '$10 >= 0.05 {print $1, $2, $3, $4, $10}' ${data}'_pruned.sex_ibd.genome' > ${data}'_pruned_ibd_pairs_pihat.txt' 
## flag twins/sibs
awk '$5 >= 0.9 {print $0}' ${data}'_pruned_ibd_pairs_pihat.txt' > ${data}'_pruned_ibd_MZpairs.txt'
awk '$5 >= 0.4 && $5 < 0.8 {print $0}' ${data}'_pruned_ibd_pairs_pihat.txt' > ${data}'_pruned_ibd_sibDZpairs.txt'
# flag outliers
awk '$5 >= 0.1875 {print $1, $2, $3, $4}' ${data}'_pruned_ibd_pairs_pihat.txt'  > ${data}'_pruned_ibd_outlier_pairs.txt'
out_n=$(wc ${data}'_pruned_ibd_outlier_pairs.txt' -l | awk '{print $1}')
if [ ${out_n} -gt 1 ]
then 
 Rscript ${geno_scripts}/selectRandom_ibdOutliers.R ${data}'_pruned_ibd_outlier_pairs.txt' # will generate ${data}'_pruned_ibd_outliers.txt'
 wc *'_pruned_ibd_outliers.txt'
fi

if [ -f ${data}'_pruned_ibd_outliers.txt' ]
 then
  plink --bfile ${data}'_sexchecked' --remove ${data}'_pruned_ibd_outliers.txt' --make-bed --out ${data}'_pruned.sex.no_close_relatives'
  else
  mv ${data}'_sexchecked.bim' ${data}'_pruned.sex.no_close_relatives.bim'
  mv ${data}'_sexchecked.bed' ${data}'_pruned.sex.no_close_relatives.bed'
  mv ${data}'_sexchecked.fam' ${data}'_pruned.sex.no_close_relatives.fam'
fi

# Check there are no relatives left!
plink --bfile ${data}'_pruned.sex.no_close_relatives' --genome --out ${data}'_pruned.sex.no_close_relatives_ibd'
awk '$10 >= 0.1875 {print $1, $2}' ${data}'_pruned.sex.no_close_relatives_ibd.genome'

# Calculate average IBD per individual using R, output outliers (defined as more than sigma standard deviations above the mean, as provided by the user):
R --file=${gwas_scripts}/IndividualIBD.R --args ${data}'_pruned.sex_ibd.genome' 5
R --file=${gwas_scripts}/IndividualIBD.R --args ${data}'_pruned.sex.no_close_relatives_ibd.genome' 5

rm  ${data}'_pruned.sex.no_close_relatives_ibd.genome' # remove genome file, which is huge

# Exclude outliers from both LD-stripped and all SNP binary files
#plink --bfile ${data}'_pruned.sex.no_close_relatives' --remove ${data}'_pruned.sex.no_close_relatives_ibd.genome.IBD_INDIV_outliers.txt' --make-bed --out ${data}'.LD_IBD' # 


fi
#----------------------------------------------------------------------

if [ ! -f ${data}'_het.ibc' ]
then
echo "Test for unusual patterns of genome-wide heterogeneity in LD-pruned data"
echo "----------------------------------------------------------------------"
# Test for unusual patterns of genome-wide heterogeneity in LD-pruned data
plink --bfile ${data}'_sexchecked' --ibc --out ${data}'_het'

echo "Exclude samples identified as outliers"
R --file=${gwas_scripts}IdHets.R --args ${data}'_het'

#----------------------------------------------------------------------
fi

if [ ! -f ${data}.pop_strat.eigenstratgeno ]
then
echo "Population stratification by principal component analysis: all samples included"
echo "----------------------------------------------------------------------"
# Run EIGENSOFT using LD-pruned binary
# Convert files to EIGENSOFT format using CONVERTF
## Requires par file to convert from packedped format to eigenstrat format
convertf -p <(printf "genotypename: "${working_dir}/${data}"_sexchecked.bed
snpname:  "${working_dir}/${data}"_sexchecked.bim
indivname:  "${working_dir}/${data}"_sexchecked.fam
outputformat: EIGENSTRAT
genotypeoutname:  "${data}".pop_strat.eigenstratgeno
snpoutname: "${data}".pop_strat.snp
indivoutname: "${data}".pop_strat.ind")

fi

echo "SmartPCA"
# Run SmartPCA, removing no outliers
## Produces 100 PCs ## too many for such a small sample size?
# m: outlier removal interactions 
# t: number of components from which outliers should be removed 
# k: number of components to output 
# s: minimum number of standard dev from the mean of each component an individual must be to be counted as an outlier.
if [ ! -f ${data}.pop_strat.pca.evec ]
then
smartpca.perl \
-i ${data}.pop_strat.eigenstratgeno \
-a ${data}.pop_strat.snp \
-b ${data}.pop_strat.ind \
-o ${data}.pop_strat.pca \
-p ${data}.pop_strat.plot \
-e ${data}.pop_strat.eval \
-l ${data}.pop_strat_smartpca.log \
-m 0 -t 100 -k 100 \
-s 6
# Minor edit to allow import into R
##Remove leading tab and split ID into two columns.
sed -i -e 's/^[ \t]*//' -e 's/:/ /g' ${data}.pop_strat.pca.evec
# Plot PCs, 1 vs 2, 2 vs 3, 3 vs 4
R --file=${gwas_scripts}/PlotPCs.R --args ${data}.pop_strat 1 2
R --file=${gwas_scripts}/PlotPCs.R --args ${data}.pop_strat 2 3
R --file=${gwas_scripts}/PlotPCs.R --args ${data}.pop_strat 3 4
fi

#Run SmartPCA again to remove outliers: skipped this step, because there is clearly a lot of structure (not a homogeneous population) so cannot exclude outliers as such... will do so for different ancestry groups

echo "----------------------------------------------------------------------"
echo "Visualize ancestry grouping - 1000G "
echo "----------------------------------------------------------------------"
# Use 1000 Genomes Phase 1 data to achieve the assess ancestry grouping, rerun without pruned samples
cat ${data}'.sexcheck_mismatch_samples.txt' \
 ${data}'_pruned.sex.no_close_relatives_ibd.genome.IBD_INDIV_outliers.txt' \
 ${data}'_het.LD_het_outliers_sample_exclude' | grep -v ID | \
 sort | uniq  > ${data}_sex_ibd_het_excluded.outliers
# do not remove related individuals:  ${data}'_pruned_ibd_outliers.txt', so that also related indviduals get the ancestry grouping (e.g. for GWAS)

plink --bfile ${data}'.hw_dropped' --remove ${data}_sex_ibd_het_excluded.outliers --make-bed --out ${data}_cleaned

# exclude SNPs within highLD regions or not autosomes
plink --bfile  ${data}_cleaned --exclude  ${working_dir}/highLD_and_autosomal_excludes --make-bed  -out ${data}_cleaned_excHighLD_autosomes
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")|| \
($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' ${data}_cleaned_excHighLD_autosomes.bim | grep -v ambig > ${data}_cleaned_excHighLD_autosomes.snplist.txt

echo "Use 1000 Genomes Phase 1 data to achieve the assess ancestry grouping, rerun without pruned samples"
# Extract rs IDs from 1KG (and add phenotypes)
plink --bfile ${ref_1KG}/1kg_phase1_all --pheno ${ref_1KG}/1KG_Phenos.txt --extract ${data}_cleaned_excHighLD_autosomes.snplist.txt --make-bed --out 1kg_phase1_all.rsids.autosomal #

# Obtain SNPs present in both files
[ -f 1kg_phase1_all.rsids_names.txt ] && rm 1kg_phase1_all.rsids_names.txt
awk '{print $2}' 1kg_phase1_all.rsids.autosomal.bim > 1kg_phase1_all.rsids_names.txt

# Extract 1KG SNPs from data, and recode to 12 alleles, in order to match it to the external 1KG allele definitions
plink --bfile  ${data}_cleaned_excHighLD_autosomes --extract 1kg_phase1_all.rsids_names.txt --make-bed --out ${data}_cleaned_1KGsnps
# Dry run bmerge to identify SNPs PLINK will fail on
plink --bfile ${data}_cleaned_1KGsnps --bmerge 1kg_phase1_all.rsids.autosomal --merge-mode 6 --out ${data}_cleaned_1KGsnps_failures
# flip failed items
plink --bfile ${data}_cleaned_1KGsnps --flip ${data}_cleaned_1KGsnps_failures.missnp --make-bed --out ${data}_cleaned_1KGsnps_flipped
# and merge
plink --bfile ${data}_cleaned_1KGsnps_flipped --bmerge 1kg_phase1_all.rsids.autosomal --merge-mode 6 --out ${data}_cleaned_1KGsnps_failures2

# remove if exists
[ -f ${data}_cleaned_1KGsnps_failures.multiple.positions.txt ] && rm ${data}_cleaned_1KGsnps_failures.multiple.positions.txt
#Add variants with multiple positions to missnp
fgrep \'rs ${data}_cleaned_1KGsnps_failures2.log |
grep 'Multiple' |\
awk '{print $7}' |\
sed -e "s/the//g" -e "s/\\.//g" -e "s/'//g"  > ${data}_cleaned_1KGsnps_failures.multiple.positions.txt

# skip if multiple.possitions.txt or .missnp is empty
cat ${data}_cleaned_1KGsnps_failures2.missnp ${data}_cleaned_1KGsnps_failures.multiple.positions.txt > ${data}_cleaned_1KGsnps_failures.multiple.positions.missnp
# Exclude mismatched SNPs and variants with multiple positions, if multiple.positions is not empty
plink --bfile ${data}_cleaned_1KGsnps_flipped --exclude ${data}_cleaned_1KGsnps_failures2.missnp --make-bed --out ${data}_cleaned_1KGsnps_4merge

# Merge
plink --bfile ${data}_cleaned_1KGsnps_4merge --bmerge 1kg_phase1_all.rsids.autosomal --make-bed --out ${data}'_1KG'

# Filter missing variants, rare variants and HWE
plink --bfile ${data}'_1KG' --geno 0.01 --maf 0.05 --hwe 0.0001 --make-bed --out ${data}'_1KG.for_prune'
# LD Pruning - merged file
plink --bfile ${data}'_1KG.for_prune' --indep-pairwise 2000 250 0.1 --out ${data}'_1KG.prune'
plink --bfile ${data}'_1KG.for_prune' --extract ${data}'_1KG.prune.prune.in' --maf 0.1  --make-bed --out ${data}'_1KG.pruned'
#----------------------------------------------------------------------

# run MDS on merged file, following ENIGMA
plink --bfile ${data}'_1KG.pruned' --cluster --mind .05 --mds-plot 4 --out ${data}'_pruned_1KG_mds'
[ -f ${data}'_pruned_1KG_mds2R.mds' ] && rm ${data}'_pruned_1KG_mds2R.mds'
awk '{print $1, $2, $3, $4, $5, $6, $7}' >> ${data}'_pruned_1KG_mds2R.mds' ${data}'_pruned_1KG_mds.mds'
# specify input file and output prefix
if [ ! -f PopulationCodes.csv ]
then
 cp ${ref_1KG}PopulationCodes.csv ./
fi
Rscript ${geno_scripts}/PlotMDS_1KG_ENIMGA2.R ${data}'_pruned_1KG_mds2R.mds' ${data}'_1KG.pruned' ${data}'_pruned_1KG_mds' ${data}'.pop_strat_outliers.outliers' ${data}

#c("ABCD_release_3.0_QCed_pruned_1KG_mds2R.mds","ABCD_release_3.0_QCed_1KG.pruned","ABCD_release_3.0_QCed_pruned_1KG_mds","ABCD_release_3.0_QCed.pop_strat_outliers.outliers","ABCD_release_3.0_QCed")


#----------------------------------------------------------------------
# for Eigensoft, following Coleman et al.
if [ ! -f "${data}"_1KG.pruned_pop_strat.pca.evec_RENAMED ]
then
# edit fam file so that the in-house dataset has it's own label (instead of -9):
sed -i 's/\-9$/100$/g' "${data}"_1KG.pruned.fam
sed -i 's/100\$/100/g' "${data}"_1KG.pruned.fam
# Run convertf to make EIGENSTRAT file
convertf -p <(printf "genotypename: "${data}"_1KG.pruned.bed
             snpname: "${data}"_1KG.pruned.bim
             indivname: "${data}"_1KG.pruned.fam
             outputformat: EIGENSTRAT
             genotypeoutname: "${data}"_1KG.pruned.eigenstratgeno
             snpoutname: "${data}"_1KG.pruned.snp
             indivoutname: "${data}"_1KG.pruned.ind")

# Generate poplist for projection
awk '{print $3}' ${ref_1KG}/1KG_Phenos.txt | sort | uniq > ${data}'_1KG.pruned.poplist.txt'

# Run Smartpca, projecting on 1KG samples only
smartpca.perl \
-i "${data}"_1KG.pruned.eigenstratgeno \
-a "${data}"_1KG.pruned.snp \
-b "${data}"_1KG.pruned.ind \
-o "${data}"_1KG.pruned_pop_strat.pca \
-p "${data}"_1KG.pruned_pop_strat.plot \
-e "${data}"_1KG.pruned_pop_strat.eigenvalues \
-l "${data}"_1KG.pruned_pop_strat.log \
-w "${data}"_1KG.pruned.poplist.txt \
-m 0

# Note that the command below relabels the phenotype column as xCHANGE, where x is the phenotype, and then relabels the 1KG populations with their names for graphing. Modify the sed command to allow your samples to be labelled usefully!
# Modify $root.1kg.LD_pop_strat.pca.evec for R
awk 'NR > 1  {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "CHANGE"$12"CHANGE"}' "${data}"_1KG.pruned_pop_strat.pca.evec  >  "${data}"_1KG.pruned_pop_strat.pca.evec_RENAMED
sed -i -e 's/CHANGE3CHANGE/ASW/g' -e 's/CHANGE4CHANGE/CEU/g' -e 's/CHANGE5CHANGE/CHB/g' -e 's/CHANGE6CHANGE/CHS/g' -e 's/CHANGE7CHANGE/CLM/g' \
 -e 's/CHANGE8CHANGE/FIN/g' -e 's/CHANGE10CHANGE/GBR/g' -e 's/CHANGE11CHANGE/IBS/g' -e 's/CHANGE12CHANGE/JPT/g' -e 's/CHANGE13CHANGE/LWK/g' \
 -e 's/CHANGE14CHANGE/MXL/g' -e 's/CHANGE15CHANGE/PUR/g' -e 's/CHANGE16CHANGE/TSI/g' -e 's/CHANGE17CHANGE/YRI/g' -e "s/CHANGE100CHANGE/${data}/g" \
 "${data}"_1KG.pruned_pop_strat.pca.evec_RENAMED
# _Plot PCs_
Rscript ${gwas_scripts}/PC_Plot_1KG.R "${data}"_1KG.pruned_pop_strat.pca.evec_RENAMED
fi

if [ ! -f "${data}".pruned_pop_strat3.phy_RENAMED ]
then

# edit par file and rerun to get the weight of each SNP into the PCs 
# generate new lines to add into par file
echo snpweightoutname: ${data}.pruned_pop_strat.weighted.snps > weight
echo phylipoutname:  ${data}.pruned_pop_strat.phy > phy
# add lines to par file:
cat "${data}"_1KG.pruned_pop_strat.pca.par weight > "${data}"_1KG.pruned_pop_strat.pca.par2 
sed -i 's/_pop_strat./_pop_strat2./g' "${data}"_1KG.pruned_pop_strat.pca.par2
#
cat "${data}"_1KG.pruned_pop_strat.pca.par phy > "${data}"_1KG.pruned_pop_strat.pca.par3 
sed -i "s/poplistname: ${data}_1KG.pruned.poplist.txt//g" "${data}"_1KG.pruned_pop_strat.pca.par3 
sed -i 's/_pop_strat./_pop_strat3./g' "${data}"_1KG.pruned_pop_strat.pca.par3 
# run smartpca using new par files
smartpca -p ${data}_1KG.pruned_pop_strat.pca.par2 > "${data}"_1KG.pruned_pop_strat_2.log
smartpca -p ${data}_1KG.pruned_pop_strat.pca.par3 > "${data}"_1KG.pruned_pop_strat_3.log
# rename phy results 
awk 'NR > 1 {print "CHANGE"$1"CHANGE", $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' "${data}".pruned_pop_strat3.phy > "${data}".pruned_pop_strat3.phy_RENAMED
sed -i -e 's/CHANGE3CHANGE/ASW/g' -e 's/CHANGE4CHANGE/CEU/g' -e 's/CHANGE5CHANGE/CHB/g' -e 's/CHANGE6CHANGE/CHS/g' -e 's/CHANGE7CHANGE/CLM/g' -e 's/CHANGE8CHANGE/FIN/g' -e 's/CHANGE10CHANGE/GBR/g' -e 's/CHANGE11CHANGE/IBS/g' -e 's/CHANGE12CHANGE/JPT/g' \
 -e 's/CHANGE13CHANGE/LWK/g' -e 's/CHANGE14CHANGE/MXL/g' -e 's/CHANGE15CHANGE/PUR/g' -e 's/CHANGE16CHANGE/TSI/g' -e 's/CHANGE17CHANGE/YRI/g' -e "s/CHANGE100CHANGE/${data}/g" "${data}".pruned_pop_strat3.phy_RENAMED
fi 
#----------------------------------------------------------------------
fi

# run R script to define ancestry subject groups: sample_ancestry_check.R
Rscript ${geno_scripts}sample_ancestry_groups.R /export/home/acarrion/acarrion/projects/resources/datasets/ABCD/data/primary_data/ABCDStudyNDAv3/acspsw03.txt ${data}.updated_genders.fam
dos2unix *samples

## Run smartpca again, including single ancestry subsamples
for ancestry in EUR #EAS AFR AMR 
 do
 if  [ $(wc ancestry_0.9${ancestry}.samples -l | awk '{print $1}') -gt 100 ]
 then
 echo run

 if [ ! -f ${data}_${ancestry}'.LD.fam' ]
 then
 plink --bfile ${data}'_sexchecked' --keep ancestry_0.9${ancestry}.samples --make-bed --out ${data}_${ancestry}'.LD' # file for PCA
 fi
 
 if [ ! -f ${data}_${ancestry}.pop_strat.eigenstratgeno ]
 then
 # Convert files to EIGENSOFT format using CONVERTF
 convertf -p <(printf "genotypename: "${working_dir}/${data}_${ancestry}".LD.bed
 snpname:  "${working_dir}/${data}_${ancestry}".LD.bim
 indivname:  "${working_dir}/${data}_${ancestry}".LD.fam
 outputformat: EIGENSTRAT
 genotypeoutname:  "${data}"_"${ancestry}".pop_strat.eigenstratgeno
 snpoutname: "${data}"_"${ancestry}".pop_strat.snp
 indivoutname: "${data}"_"${ancestry}".pop_strat.ind")
 fi
 
 for sd in 3 4 5 6
 do
 if [ ! -f  ${data}_${ancestry}_${sd}sd.pop_strat_outliers_smartpca.log ]
 then
 # Run smartpca
 smartpca.perl \
-i ${data}_${ancestry}.pop_strat.eigenstratgeno \
-a ${data}_${ancestry}.pop_strat.snp \
-b  ${data}_${ancestry}.pop_strat.ind \
-o  ${data}_${ancestry}_${sd}sd.pop_strat_outliers.pca \
-p  ${data}_${ancestry}_${sd}sd.pop_strat_outliers.plot \
-e  ${data}_${ancestry}_${sd}sd.pop_strat_outliers.eval \
-l  ${data}_${ancestry}_${sd}sd.pop_strat_outliers_smartpca.log \
-m 5 \
-t 25 \
-k 100 \
-s ${sd}
 
sed -i -e 's/^[ \t]*//' -e 's/:/ /g' ${data}_${ancestry}_${sd}sd.pop_strat_outliers.pca.evec
# Plot PCs after extracting outliers
R --file=${gwas_scripts}/PlotPCs.R --args ${data}_${ancestry}_${sd}sd.pop_strat_outliers 1 2
R --file=${gwas_scripts}/PlotPCs.R --args ${data}_${ancestry}_${sd}sd.pop_strat_outliers  2 3
R --file=${gwas_scripts}/PlotPCs.R --args ${data}_${ancestry}_${sd}sd.pop_strat_outliers  3 4

awk '/REMOVED/ {print $3}'  ${data}_${ancestry}_${sd}sd.pop_strat_outliers_smartpca.log | sed 's/:/ /g' > ${data}_${ancestry}_${sd}sd.pop_strat_outliers.outliers
fi
done

 fi
done



#----------------------------------------------------------------------
fi

#----------------------------------------------------------------------
# Save clean files
#----------------------------------------------------------------------
# Remove outlier samples from files for analyses
# all outliers, to exclude
cat ${data}'.sexcheck_mismatch_samples.txt' \
 ${data}'_pruned.sex.no_close_relatives_ibd.genome.IBD_INDIV_outliers.txt' \
 ${data}'_het.LD_het_outliers_sample_exclude' | \
 grep -v ID | \
 sort | uniq  > ${data}_sex_ibd_het_excluded.outliers
 
# ancestry specific, but all
plink --bfile ${data}'.hw_dropped2' --remove ${data}_sex_ibd_het_excluded.outliers --make-bed --out ${data}_all.cleaned
# exclude related individuals
plink --bfile ${data}_all.cleaned --remove ${data}'_pruned_ibd_outliers.txt' --make-bed --out ${data}_all_unrelated.cleaned

 
#----------------------------------------------------------------------
# to remove related individuals: ${data}'_pruned_ibd_outliers.txt' 
# to remove ancestry specific outliers: 
ancestry=EUR
for sd in 3 6
do

file2keep=ancestry_0.9${ancestry}.samples
file2rm=${data}_sex_ibd_het_${ancestry}${sd}sd_excluded.outliers

# heterozygosity outliers for subset
if [ ! -f ${data}'_het_'${ancestry}${sd}'sd.ibc' ]
then
echo "Test for unusual patterns of genome-wide heterogeneity in LD-pruned data"
# Test for unusual patterns of genome-wide heterogeneity in LD-pruned data
plink --bfile ${data}'_sexchecked' --keep ${file2keep} --remove ${data}'_'${ancestry}'_'${sd}'sd.pop_strat_outliers.outliers' --ibc --out ${data}'_het_'${ancestry}${sd}'sd'



echo "Exclude samples identified as outliers"
R --file=${gwas_scripts}IdHets.R --args ${data}'_het_'${ancestry}${sd}'sd'
fi

# all outliers, to exclude
cat ${data}'.sexcheck_mismatch_samples.txt' \
 ${data}'_pruned.sex.no_close_relatives_ibd.genome.IBD_INDIV_outliers.txt' \
 ${data}'_het.LD_het_outliers_sample_exclude' \
 ${data}'_'${ancestry}'_'${sd}'sd.pop_strat_outliers.outliers' \
 ${data}'_het_'${ancestry}${sd}'sd.LD_het_outliers_sample_exclude' | \
 grep -v ID | sort | uniq  > ${file2rm}

 #${data}'_'${ancestry}'_'${sd}'sd.pop_strat_outliers.outliers'

# ancestry specific, but all (including related)
plink --bfile ${data}'.hw_dropped2' --keep ${file2keep} --remove ${file2rm} --make-bed --out ${data}_${ancestry}.${ancestry}${sd}sd.cleaned
# exclude related individuals
plink --bfile ${data}_${ancestry}.${ancestry}${sd}sd.cleaned --remove ${data}'_pruned_ibd_outliers.txt' --make-bed --out ${data}_${ancestry}.${ancestry}${sd}sd_unrelated.cleaned

done

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# organize output files
#----------------------------------------------------------------------
# only if final files exist
if [ -f ${data}_all.cleaned.bed ]
then
cd ${working_dir}

# clean files
mkdir -p clean
mv *.cleaned.* clean/

# intermediate files
mkdir -p intermediate_files/
mv ${data}* intermediate_files/
rm -f weight phy
cd intermediate_files
rm -f *nosex *hh
rm -f $(ls *ps *eigenstratgeno *snp) # gnuplot not installed
## sampleQC
mkdir -p sampleQC/PopStr/1KG sampleQC/relatedness
mv ancestry*samples ./sampleQC/PopStr/
mv *EUR* ./sampleQC/PopStr
mv *1KG* *1kg* ./sampleQC/PopStr/1KG
mv *pop_str* ./sampleQC/PopStr
mv *imiss ./sampleQC
mv *genome* *IBD* *ibd* *relatives* ./sampleQC/relatedness
mv *sexcheck* ./sampleQC
mv *het* ./sampleQC
mv *pruned_LD* ./sampleQC
## snpQC
mkdir -p snpQC/
mv *lmiss ./snpQC
mv *frq *hwe *hw* ./snpQC
mv *filtered* *raw* *SNP* ./snpQC

## clean some intermediate files, keep log files for reference
##
rm -f $(ls sampleQC* | grep -e bim -e bed -e fam -e nosex -e hh )
rm -f $(ls ./sampleQC/PopStr/1KG/* | grep -e bim -e bed -e fam -e nosex -e hh )
rm -f $(ls ./sampleQC/relatedness/* | grep -e bim -e bed -e fam -e nosex -e hh -e genome | grep -v -e IBD_INDIV )
## 
rm -f $(ls ./snpQC/* | grep -e bim -e bed -e fam -e nosex -e hh )
# keep only:
## updated genders
rm -f $(ls * | grep -e bim -e bed -e fam -e nosex -e hh | grep -v updated_genders)

#----------------------------------------------------------------------
fi

else
echo "Final files already exist, skip."
fi
