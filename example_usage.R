#install.packages(c("BiocManager", "plsdof","MASS","NbClust","cli","progress"), repos="http://cran.us.r-project.org")
#BiocManager::install("NMF")
install.packages(".",type="source",repos=NULL)
library("gen3DNet")

histon_path <- "data/histon_data.csv"
phospho_path <- "data/phospho_data.csv"

result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,
    p_val_threshold=0.01, 
    k_picker=max_ward_kl
)
