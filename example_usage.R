#Console usage
#rm(list=ls())
#devtools::load_all()
#library("model3D",lib="..")
#histon_data<-read.table("../60histone_24hr.txt", sep="\t", header=T, row.names=1)
#phospho_data<-read.table("../peptide_24hr_final.txt", sep="\t", header=T, row.names=1)
#write.csv(histon_data, histon_path)
#write.csv(phospho_data, phospho_path)

install.packages(".",type="source",repos=NULL)
library("gen3DNet")

histon_path <- "histon_data.csv"
phospho_path <- "phospho_data.csv"

result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,
    k_picker=max_cophenetic
)
#result <- gen3DNet(
#    histon_path, 
#    phospho_path,
#    nmf_nrun=10,
#    k=4,#2:10,
#    k_picker=nbclust_pick_k,
#    seed=123456,
#    out_folder="model3d_4",
#    p_threshold=.001
#)
