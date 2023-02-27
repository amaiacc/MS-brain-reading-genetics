library(MungeSumstats)
# If your GWAS summary statistics file of interest relates to GRCh37, you will need to install 
# SNPlocs.Hsapiens.dbSNP144.GRCh37 and BSgenome.Hsapiens.1000genomes.hs37d5 from Bioconductor as follows:
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
# library(BSgenome)

library(dplyr)

# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
input_dir=paste0(dir,"/resources/datasets/GWAS_sumstats/downloaded_data/ieu-openGWAS/")
output_dir=paste0(dir,"/resources/datasets/GWAS_sumstats/ldsc/ieu-openGWAS/") # directory for ieu-open GWAS project
if(!dir.exists(output_dir)){dir.create(output_dir)}

path2files <- paste0(input_dir,list.files(input_dir,pattern="vcf.gz$"))

for (ifile in path2files){
  ofile=ifile %>% gsub(input_dir,output_dir,.) %>% gsub("vcf.gz","tsv.gz", .)
  if(file.exists(ofile)==FALSE){
    cat('Format file:',ofile,"\n")
    reformatted_vcf <- MungeSumstats::format_sumstats(path=ifile,
                                                      save_path = ofile,
                                                      return_data = FALSE,
                                                      ldsc_format = TRUE,
                                                      log_mungesumstats_msgs=TRUE,
                                                      log_folder = output_dir,
                                                      ref_genome="GRCh37")
  } else {
    cat('Output file already exists, skip.\n')
  }
  gc()
}
