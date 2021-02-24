# ----------------------------------------------------------------------------
# Author: Shamim Mollah  
# Created: 12-10-2016
# 
# Generate functional modules using NMF
#-----------------------------------------------------------------------------


#' Scale 0 to 1
#'
#' This function scales data from 0 to 1.
#' 
#' @param x the data to be scaled from 0 to 1.
scale_0_1 <- function(x) {
     a <- min(x) 
     b <- max(x) 
     (x - a)/(b - a) 
}

#' User picks k
#'
#' This is a k_picker callback that simply plots the performance
#' of NMF on every number of clusters in k_range, and asks the
#' user for input.
#'
#' @param nmf_data The data whose number of clusters is in question
#' @param k_range The numbers of clusters (k values) being considered.
#' 
#' @export
user_pick_k <- function(nmf_data, k_range) {
    V.random <- NMF::randomize(nmf_data)
    estim.r.random <- NMF::nmf(V.random, k_range, nrun = 10, seed = "nndsvd")
    estim.r <- NMF::nmf(nmf_data, k_range, nrun = 10, seed = "nndsvd")
    plot(estim.r, estim.r.random)
    print("Please pick k")
    as.numeric(readline())
}

#' Normalize NMF
#'
#' This function applies an antilog transformation to the data,
#' ensuring that all values are positive, and then scales the
#' result from 0 to 1.
#'
#' @param left_data The data to be normalized
normalize_nmf <- function(left_data) {
    YNormZ <- scale(as.matrix(left_data[,]))
    # centered data
    hi = as.data.frame(YNormZ)
    #2^ log fc value
    YNormZ_anti = as.data.frame(sapply(hi, function(x) 2^x))
    #[0-1] range
    apply(YNormZ_anti, 2, scale_0_1)
}

#how to pass seed to NMF if seed = "nndsvd"

#' Generate NMF Modules
#'
#' This function uses NMF to generate functional modules representing interactions
#' between the entities in the rows and columns of the dataframe. This has two
#' components:
#' Each entity is assigned to a cluster
#' "Column signatures" are generated for each row object. For example, if rows are histones and columns are drugs, this function will generate "histone signatures" that express
#' Because NMF relies on knowing the number of clusters (k) in advance, this function also determines this automatically.
#' Currently, the default option uses the number of clusters that produces the lowest KL index with Ward clustering. The package depends on NbClust to implement this feature. However, the user can override this option with a custom function.
#'
#' @param left_data This is the data that will be clustered using NMF. Signatures will be generated for the object type that is represented by rows. 
#' @param nmf_nrun The number of r
#' @param k_range Either a) A consecutive range of possible k values, b) a single k value, or c) nothing, meaning that generate_nmf_modules will assume the range is from 2 to the number of rows.
#' @param k_picker Any function that takes a dataframe and a consecutive range of potential k values
#' @param seed The seed to use with NMF
#' @param verbose Whether to print output.
generate_nmf_modules <- function(left_data, nmf_nrun, k_range, k_picker=gen3DNet::max_kl_ward, seed, verbose=FALSE) {

    cli::cli_alert_success("Running NMF")

    left_data <- t(left_data)

    mat_hist <- normalize_nmf(left_data)

    if (verbose) {
        if (length(k_range) != 1) {
            cli::cli_alert_info(paste(" Picking k from [", min(k_range), ",", max(k_range), "]."))
        }
    }

    k <- k_picker(left_data, k_range)

    if (verbose) {
        cli::cli_alert_info(paste(" Using k =", k))
    }

    #with nndsvd seed
    res_hist <- NMF::nmf(mat_hist, k, "brunet", nrun=nmf_nrun, seed = "nndsvd")

    #Add NMF plot of basis and loading heatmaps
    pdf(file="basis_loading_heatmap")

    layout(cbind(1, 2))
    basismap(res_hist, labRow= row.names(left_data))
    coefmap(res_hist)

    dev.off()

    res_hist
}
