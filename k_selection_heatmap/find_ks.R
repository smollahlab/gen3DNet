load("datasets")
library("NbClust")
source("../generate_cluster_data.R")
source("../generate_nmf_modules_functions.R")

measures <- lapply(
    1:600,
    function(x) {
        load(paste("measures",x,sep=""))
        measures
    }
)

print("Cophenetic")

kneedle_cophenetic <- lapply(
    1:length(measures),
    function(i) {
        print(i)
        measure <- measures[[i]]
        (2:10)[kneedle(measure$cophenetic, 1)]
    }
)

save(kneedle_cophenetic, file="kneedle_cophenetic")

print("Silhouette")

kneedle_silhouette <- lapply(
    1:length(measures),
    function(i) {
        print(i)
        measure <- measures[[i]]
        (2:10)[kneedle(measure$silhouette.consensus, 1)]
    }
)

save(kneedle_silhouette, file="kneedle_silhouette")

print("NbClust")

nbclust <- lapply(1:length(datasets),
   function(i) {
        print(i)
        dataset <- datasets[[i]]
        NbClust(t(dataset), diss=NULL, distance = "euclidean", min.nc=2, max.nc=10, 
            method = "ward.D2", index = "kl")$Best.nc["Number_clusters"]
   }
)
save(nbclust, file="nbclust")
