load("kneedle_cophenetic")
load("parameters")
load("kneedle_silhouette")
load("nbclust")

#Needs refactoring

library(dplyr)
library(tidyr)
library(maditr)

parameters$kneedle_cophenetic <- kneedle_cophenetic
parameters$kneedle_silhouette <- kneedle_silhouette
parameters$nbclust <- nbclust

names(parameters)

print(parameters)

print(cor(as.numeric(parameters$nbclust), parameters$k))

#                     NBCLUST

pivot <- parameters %>% 
    select(k, nbclust, variances, nrows) %>%
    group_by(variances, nrows) %>%
    summarise(correlation=cor(k, as.numeric(nbclust)))

pivot$variances <- as.character(pivot$variances)
pivot <- data.frame(pivot)
nbclust_pivot <- reshape(pivot, idvar="nrows", timevar="variances", direction="wide")
print(pivot)
save(nbclust_pivot, file="nbclust_pivot")

#                   KNEEDLE SILHOUETTE

pivot <- reshape(pivot, idvar="nrows", timevar="variances", direction="wide")
print(pivot)

pivot <- parameters %>%
    select(k, kneedle_silhouette, variances, nrows) %>%
    group_by(variances, nrows) %>%
    summarise(correlation=cor(k, as.numeric(kneedle_silhouette)))

pivot$variances <- as.character(pivot$variances)
pivot <- data.frame(pivot)
print(pivot)
kneedle_silhouette_pivot <- reshape(pivot, idvar="nrows", timevar="variances", direction="wide")
save(kneedle_silhouette_pivot,file="kneedle_silhouette_pivot")
print(pivot)

#                   KNEEDLE COPHENETIC

pivot <- parameters %>%
    select(k, kneedle_cophenetic, variances, nrows) %>%
    group_by(variances, nrows) %>%
    summarise(correlation=cor(k, as.numeric(kneedle_cophenetic)))

pivot$variances <- as.character(pivot$variances)
pivot <- data.frame(pivot)
print(pivot)
kneedle_cophenetic_pivot <- reshape(pivot, idvar="nrows", timevar="variances", direction="wide")
print(pivot)
save(kneedle_cophenetic_pivot, file="kneedle_cophenetic_pivot")
