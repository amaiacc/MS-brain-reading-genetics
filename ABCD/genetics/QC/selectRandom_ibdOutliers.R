# random selection of individual to exclude from the IBD pairs, pi-hat > 0.1

#ibd_file<-list.files(pattern="_pruned_ibd_outlier_pairs.txt")
args <- commandArgs(TRUE)
ibd_file<-args[1]

root_name<-gsub("_pruned_ibd_outlier_pairs.txt","",ibd_file)
ibd_table<-read.table(ibd_file,header=TRUE,stringsAsFactors=FALSE)


# create list to exclude with first indv from each pair
ibd_table_exclude<-ibd_table[,1:2]

# get random indices (1,2) for samples to exclude
index<-sample(c(1,2),replace=TRUE,size=NROW(ibd_table))
print(index)
# define indices, to avoid using them if their length is 0
w1<-which(index==1)
w2<-which(index==2)

# replace individual to second of pair, for every index that is 2
if (length(w2)>0){
  ibd_table_exclude[which(index==2),]
  # when index=2, replace 1st ind in pair for the 2nd indv
  ibd_table_exclude[ which(index==2),]<-ibd_table[ which(index==2),3:4] 
}
# check
if (length(w1)>0){ ibd_table_exclude[w1,]==ibd_table[w1,1:2] }
if (length(w2)>0){ ibd_table_exclude[w2,]==ibd_table[w2,3:4] }

# save file with IDs to exclude
write.table(ibd_table_exclude,file=paste(root_name,"_pruned_ibd_outliers.txt",sep=""),sep=" ",col.names=TRUE,quote=FALSE,row.names=FALSE)
