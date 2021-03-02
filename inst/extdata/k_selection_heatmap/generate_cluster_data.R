library(MASS)

get_center <- function(n_rows, size) {
  center <- mvrnorm(n_rows, mu=0, Sigma=1)
  center_as_size <- do.call(cbind, rep(list(center),size))
  center_as_size
}

get_error <- function(n_rows, size, f_stat) {
  error_as_size <- sapply(1:size,
                          function(i) {
                            mvrnorm(n_rows, mu=0, Sigma=1) / f_stat
                          })
  error_as_size
}

get_sizes <- function(n_clusters, n_cols) {
  dividend <- n_cols %/% n_clusters
  remainder <- n_cols %% n_clusters
  c(rep(dividend, n_clusters - remainder),
    rep(dividend + 1, remainder))
}

generate_n_clusters <- function(n_clusters, n_rows, n_cols, f_stat) {
  do.call(cbind,
          lapply(get_sizes(n_clusters, n_cols),
                 function(size) {
                   get_center(n_rows, size) + get_error(n_rows, size,f_stat)
                 }))
}

generate_clustering_datasets <- function(ns_clusters,
                                         n_of_each,
                                         n_rows,
                                         n_cols,
                                         f_score) {
  do.call(c,
          lapply(ns_clusters,
                 function(n_clusters){
                   lapply(1:n_of_each,
                          function(n) {
                            list(n=n_clusters, clusters=generate_n_clusters(n_clusters, n_rows, n_cols, f_score))
                          }
                   )
                 }))
}

make_list <- function(x) {
  if (is.function(x)) {
    x <- c(x)
  }
  as.list(x)
}

dapply <- function(param_sets, func, name=NULL) {
  results_length <- length(param_sets)
  all_results <- lapply(1:results_length,
                        function(i) NULL)
  completed <- lapply(1:results_length,
                      function(i) FALSE)
  results <- lapply(
    1:results_length,
    function(i) {
      parameterized <- function() {
        if (!completed[[i]]) {
          result <- do.call(func, 
                            make_list(param_sets[[i]]))
          if (!is.null(name)) {
            result <- list(result)
            names(result) <- name
          }
          all_results[[i]] <<- result
          completed[[i]] <<- TRUE
        }
        all_results[[i]]
      }
      return(parameterized)
    }
  )
  attr(results,"all_results") <- function() {
    all_results
  }
  attr(results,"completed") <- function() {
    completed
  }
  results
}



generic_genetic <- function(initializer,
                            update,
                            scorer,
                            generate,
                            n_create,
                            n_keep,
                            iterations) {
  initials <- lapply(n_create, function(i) initializer())
  scores <- lapply(initials, scorer)
  best <- order(scores)[1:n_keep]
  for (i in i:iterations) {
    next <- generate(best)
    scores <- lapply(initials, scorer)
    best <- order(scores)[1:n_keep]
  }
  best
}