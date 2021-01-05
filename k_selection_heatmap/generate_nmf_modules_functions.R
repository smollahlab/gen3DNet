# ----------------------------------------------------------------------------
# Author: Shamim Mollah  
# Created: 12-10-2016
# 
# Generate functional modules using NMF
#-----------------------------------------------------------------------------

#rm(list = ls())

library(NMF)

normalize_histon <- function(histon_data) {
  YNormZ<-scale(as.matrix(histon_data[,]))
  #Normalized Data to [0,1]
  fun <- function(x){
    a <- min(x) 
    b <- max(x)
    (x - a)/(b - a) 
  }
  # centered histone data
  hi=as.data.frame(YNormZ)
  #2^ log fc value
  YNormZ_anti=as.data.frame(sapply(hi, function(x) 2^x))
  #[0-1] range
  mat_hist <- apply(YNormZ_anti, 2, fun)
  mat_hist
}

a <- 0

#from https://stackoverflow.com/questions/35194048/using-r-how-to-calculate-the-distance-from-one-point-to-a-line
dist2d <- function(a,b,c) {
 v1 <- b - c
 v2 <- a - b
 m <- cbind(v1,v2)
 d <- det(m)/sqrt(sum(v1*v1)) #take out abs
 d
}

up_down <- function(metrics) {
    last <- metrics[1]
    for (i in 2:length(metrics)) {
        current <- metrics[i]
        if (current < last) {
          return(i-1)
        }
        last <- current
    }
}

kneedle <- function(cophenetic, sign) {
    start = c(1, cophenetic[1])
    end = c(length(cophenetic), cophenetic[length(cophenetic)])
    k <- which.max(lapply(1:length(cophenetic),
                     function(idx) {
                         sign * -1 * dist2d(c(idx, cophenetic[idx]),
                                start,
                                end
                         )
                     })
              )
    plot(1:length(cophenetic),cophenetic,type="l")
    lines(c(1,length(cophenetic)),
          c(cophenetic[1],cophenetic[length(cophenetic)]),
          col="red")
    points(k, cophenetic[k],
           col="blue")
    k
}

avg = 0
pick_k <- function(mat_hist, k_range, nrun, seed, k_finder = which.max, metric) {
  #browser()
  ##to identify optimal factorization rank
  #V.random <- randomize(mat_hist)
  #estim.r.random <- nmf(V.random, 2:10, nrun = nrun, seed = 123456)
  #estim.r <- nmf(mat_hist, 2:10, nrun = nrun, seed = 123456)
  #V.random <- randomize(mat_hist)
  estim.r.random <- nmf(mat_hist, 
                        k_range,
                        nrun = nrun,
                        seed = seed)
  #! Add estim.r and use heuristic to combine
  #  Must consider both.
  #estim.r <- nmf(mat_hist, 
  #               k_range,
  #               nrun = nrun,
  #               seed = seed)
  cophenetic1 <- estim.r.random$measures[[metric]]
  #cophenetic2 <- estim.r$measures[[metric]]
  print(cophenetic1)
  #print(cophenetic2)
  k <- k_range[k_finder(cophenetic1)]
  #l <- k_range[k_finder(cophenetic2)]
  list(k=k,#c(k,l))
       estim.r.random=estim.r.random #un-commented
  )
}

generate_nmf_modules <- function(histon_data,
                                 k_options,
                                 nrun,
                                 seed,
                                 k_finder=which.max,
                                 metric="cophenetic"
                                 ) {
  browser()
  
  print("Normalizing...")

  mat_hist <- normalize_histon(histon_data)

  print("Picking k...")
  
  #switch(length(k_options))
  #case '1'
  #Switch case more elegant
  #https://stackoverflow.com/questions/10393508/how-to-use-the-switch-statement-in-r-functions
  
  
  if (length(k_options) == 1) {
    k <- k_options[1]
  }
  else {
    k_results <- pick_k(mat_hist, k_options, nrun, seed, k_finder, metric)
    k <- k_results[1]
    l <- k_results[2]
  }
  estim.r.random = k_results$estim.r.random
  
  print(paste("chosen k is", k))
  #! Separate into different function, from this point on
  print("finding final W, H")
  # with nndsvd seed
  res_hist <- nmf(mat_hist, k, "brunet", nrun=nrun, seed = "nndsvd")

  browser()
  list(res_hist=res_hist,
       estim.r.random=estim.r.random,
       k=k)
}
