#How good is pick_k?

#source("../generate_nmf_modules_functions.R")
source("../generate_cluster_data.R")
lib="/home/paul.morrison/bids/iPhDNet/"
#library(NMF,lib=lib)
#library(MASS,lib=lib)
library(NMF)
library(MASS)
k_range <- 2:10 #Low for testing
n_nmf_run <- 10 #Low for testing
n_gen_run <- 10 #Low for testing
seed <- 123
#k_finder <- kneedle
metric <- "cophenetic"
n_rows <- 10
n_cols <- 1000 #More cols than rows because these are clustered
f_score <- 100 #super distinct clusters
data_generator <- function(...)  2^generate_n_clusters(...)

compare_k <- function(k_range, 
                      n_nmf_run, 
                      seed,
                      k_finder, 
                      metric,
                      data_generator,
                      n_gen_run, 
                      n_rows,
                      n_cols,
                      f_score) {
  true_k <- sample(k_range, size=n_gen_run, replace=TRUE)
  sim_datasets <- lapply(true_k,
    function(k) {
      data_generator(k, n_rows, n_cols, f_score)
    }
  )
  print(true_k)
  fit <- lapply(sim_datasets,
    function(dataset) {
      pick_k(dataset, k_range, n_nmf_run, seed, k_finder, metric)
    }
  )
  pred_k <- lapply(fit,
    function(fit_i) {
      fit_i$k
    }
  )
  list(true_k=true_k,
       pred_k=pred_k,
       sim_datasets=sim_datasets,
       fit=fit
  )
}

maximize = 1
minimize = -1

metrics <- c(
    "cophenetic",
    "rss",
    "residuals",
    "evar",
    "silhouette.coef",
    "silhouette.basis",
    "silhouette.consensus",
    "dispersion"
)
signs <- c(
  maximize,
  minimize,
  minimize,
  maximize,
  maximize,
  maximize,
  maximize,
  minimize
)

k_compared <- list()

#for (i in 1:length(metrics)) {
#  k_compared[length(list)+1] <- compare_k(
#    k_range,
#    n_nmf_run, 
#    seed,
#    function(data) kneedle(data, signs[i]),
#    metrics[i],
#    data_generator,
#    n_gen_run, 
#    n_rows,
#    n_cols, 
#    f_score
#  )
#}



#k_compared_cophenetic_1000_cols <- compare_k(
#  k_range, 
#  n_nmf_run, 
#  seed,
#  which.max,
#  "sparseness.basis", 
#  data_generator, 
#  n_gen_run, 
#  n_rows, 
#  n_cols, 
#  f_score
#)

#random_k <- sample(k_range, 100, replace=TRUE)
#save(random_k, file="random_k_2")
#random_data <- dapply(
#  random_k,
#  function(k) {
#    generate_n_clusters(k, n_rows, n_cols, f_score)
#  }
#)

#random_pos <- dapply(
#  random_data,
#  function(data_func) {
#    2 ^ data_func()
#  }
#)

#nmf_results <- dapply(
#  random_pos,
#  function(data_func) {
#    nmf(data_func(), k_range, "brunet",nrun=n_nmf_run, seed=seed)
#  }
#)

#best_k <- dapply(
#  nmf_results,
#  function(nmf_result) {
#    k_range[which.max(nmf_result()$measures$sparseness.coef)]
#  }
#)

#lapply(nmf_results[1:30], function(x) x())

#for (i in 1:length(random_k)) {
#  print(i)
#  plot(random_k[1:i],
#       lapply(best_k[1:i], function(x) x()),
#       col=rgb(red=0, green=0, blue=0, alpha=.2))
#}
#save(random_data,file="random_data_2")
#lapply(1:length(random_data),
#        function(i){
#            print(i)
#            result <- nmf(
#                2^random_data[[i]](), 
#                k_range, 
#                "brunet", 
#                nrun=n_nmf_run, 
#                seed=seed, 
#                .options=list(
#                    parallel=TRUE, 
#                    verbose=TRUE
#
#                )
#            )
#            save(result, file=paste("nmf_round_two",i,sep=""))
#        })

nrows <- 2^seq(2:6)
ncols <- 1000
variances <- seq(1,11,2)
k_range <- 2:10
print(nrows)
print(ncols)
print(variances)
print(k_range)

create_parameters <- function(nrows, variances, k_range, ncol, n_k) {
    original_df <- data.frame(expand.grid(nrows=nrows, variances=variances))
    #df$k <- sample(k_range, dim(df)[1], replace=TRUE)
    #Repeat dataframe and add different k's
    print(original_df)
    print(dim(original_df))
    df <- original_df[rep(1:(dim(original_df)[1]),each=n_k),]
    df$k <- do.call(c, lapply(
        original_df$nrows,
        function(nrows) {
            sample(k_range, n_k, replace=TRUE)
        }
    ))
    #df <- do.call(rbind, rep(list(original_df), n_k))
    print(dim(df))
    print(df)
    df$ncols <- ncols
    df
}

tryCatch({    
    load("parameters")
},
error=function(cond) {
    print("creating parameters")
    parameters <<- create_parameters(nrows, variances, k_range, ncols, 20)
    #datasets <- create_datasets(parameters)
    save(parameters, file="parameters")
})

#save(parameters,file="parameters")i
#print(parameters)
#print("Parameter dimension")
#print(dim(parameters))

create_datasets <- function(parameters) {
    ncols <- parameters$ncols
    variances <- parameters$variances
    ks <- parameters$k
    nrows <- parameters$nrows
    lapply(
        1:length(ncols),
        function(i) {
#            print(ks[i])
#            print(nrows
            generate_n_clusters(ks[i], nrows[i], ncols[i], variances[i])
        }
    )
}

tryCatch({
    load("datasets")
},
error=function(cond) {
    print("creating datasets")
    datasets <<- create_datasets(parameters)
    save(datasets, file="datasets")
})

print(datasets[[1]])

for (i in 1:length(datasets))  {
    print(i)
    filename <- paste("nmf_result",i,sep="")
    tryCatch({
        load(filename)
    },
    error=function(cond) {
        print("running nmf")
        dataset <- datasets[[i]]
        nmf_result <<- nmf(
            as.matrix(2^dataset),
            k_range,
            #"brunet",
            nrun=n_nmf_run,
            #seed="nndsvd",
            .options=list(
                parallel=FALSE, 
                verbose=TRUE
            )
        )
        save(nmf_result,file=filename)
    })
}


#ncols <- c(2, 4, 8, 16, 32, 64)

#variances <- c(1, 3, 5, 7, 8, 9, 10, 11) 

#create_many_datasets(ncols, variances, "datasets")

#diskapply <- function(x, FUN, folder) {
#    in_memory <- list()
#    lapply(
#        1:length(x),
#        function(n) {
#            if (list_has(in_memory, n)) {
#                return(in_memory[[n]])
#            }
#            else {
#            }
#        }
#    )
#}

