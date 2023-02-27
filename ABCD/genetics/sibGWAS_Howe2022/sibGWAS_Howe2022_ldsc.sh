source activate ldsc

PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
ldscores_dir=${resource_dir}/LDscores/


# testing data from Howe2022 sibGWAS
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/ieu-openGWAS/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/ieu-openGWAS/
cd ${working_dir}
ids2check=$(awk '{print $1}' ${primary_dir}Howe2022_GWASlist2download.txt | grep -v id | sed 's/"//g')


# get h2 for each phenotype
for id in ${ids2check}
 do 
 # get estimates for h2
 # id=ieu-b-${i}
 if [ ! -f ${id}_h2.log ]
 then
 ldsc.py \
    --h2 ${id}.tsv.gz \
    --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
    --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
    --out ${id}_h2
 fi
 done
