
show <- function(a) {
    print(a)
    return(a)
}



dist2d <- function(a,b,c) {
 v1 <- b - c
 print(v1)
 v2 <- a - b
 print(v2)
 m <- cbind(v1,v2)
 d <- show(det(m))/show(sqrt(sum(v1*v1))) #take out abs
 show(d)
}
kneedle <- function(cophenetic, sign) {
    start = c(1, cophenetic[1])
    end = c(length(cophenetic), cophenetic[length(cophenetic)])


    k <- which.max(show(lapply(1:length(cophenetic),
                     function(idx) {
                         show(sign * -1 * dist2d(c(idx, cophenetic[idx]),
                                start,
                                end
                         ))
                     })
              ))
    #plot(1:length(cophenetic),cophenetic,type="l")
    #lines(c(1,length(cophenetic)),
    #      c(cophenetic[1],cophenetic[length(cophenetic)]),
    #      col="red")
    #points(k, cophenetic[k],
    #       col="blue")
    print(k)
    k
}

remove.na <- function(a) {

    a[!is.na(a)]
}



#' Pick k based on the maximum KL index of Ward clusterings
#' 
#' This function uses the NbClust package to 
#' 
#' @param data The data with an unknown number of clusters.
#' @param k_range The range of possible k values
#' 
#' @export
max_ward_kl <- function(data, k_range) {
  NbClust::NbClust(t(data), diss=NULL, distance = "euclidean", min.nc=min(k_range), max.nc=max(k_range),
          method = "ward.D2", index = "kl")$Best.nc["Number_clusters"]
}

#' @export
max_silhouette_consensus <- function(data, k_range) {
  k_range[which.max(nmf(data, k_range, nrun=10)$measures$silhouette.consensus)]
}

#' @export
max_cophenetic <- function(data, k_range) {
  k_range[which.max(nmf(data, k_range, nrun=10)$measures$cophenetic)]
}

#' @export
kneedle_silhouette_consensus <- function(data, k_range) {
  k_range[kneedle(show(remove.na(nmf(data, k_range, nrun=10)$measures$silhouette.consensus)), 1)]
}

#' @export
kneedle_cophenetic <- function(data, k_range) {
  k_range[kneedle(remove.na(nmf(data, k_range, nrun=10)$measures$cophenetic), 1)]
}

#' Create a general 3D network.
#' 
#' This function performs calculations to generate a 3D network.
#' It is not meant to be called directly. You should use gen3DNet
#' instead. See help(gen3DNet) for further documentation.
create_gen3DNet <- function(
    left,
    right,
    nmf_nrun=10,
    k=c(),
    k_picker=max_ward_kl, 
    seed=123456,
    p_val_threshold=.0001,
    out_folder="gen3DNet",
    verbose=TRUE
) {
    if (verbose) {
        cli::cli_alert_success("Generating 3D network")
    }

    #Simplify 'k' argument
    if (length(k) == 1) {
        #The user has provided a single k
        #No k selection is necessary
        #k_final <- k
        min_k <- k
        max_k <- k
    }
    else {
       #The user has not provided a single k
       #We must choose it automatically
       if (length(k) == 0) {
           #The user has provided no min or max k.
           #Choose min and max k automatically.
           min_k <- 1            #Min k is 1
           max_k <- min(dim(left) - 1) #Max k is num. cols or 30
       }
       else {
           #The user has provided a range of k values.
           #So first, ensure the k values are consecutive.
           if (!identical(min_k:max_k, k)) {
               stop("model3d must consider every choice from min_k to max_k")
           }
           #Select min and max
           min_k <- min(k)
           max_k <- max(k)
       }
    }

    #1. Create and enter new folder for network results
    prev_wd = normalizePath(getwd())
    dir.create(out_folder)
    setwd(out_folder)

    #2. Create NMF folder and run NMF
    dir.create("nmf")
    setwd("nmf")
    nmf_result <- generate_nmf_modules(
        left,
        nmf_nrun=nmf_nrun,
        k_range=min_k:max_k,
        k_picker=k_picker,
        seed=0,
        verbose=verbose
    )
    setwd("..")

    #3. Create PLSR folder and run PLSR
    dir.create("plsr")
    setwd("plsr")
    plsr_result <- generate_plsr(left, right, p_val_threshold, verbose=verbose)
    setwd("..")

    #4. Merge and re-group data

    p_list <- plsr_result$p_list
    non_p_list <- plsr_result$non_p_list

    if (verbose) {
        cli::cli_alert_success("Merging data...")
    }

    loading_original <- coef(nmf_result)
    loading <- t(loading_original)

    left_names <- colnames(left)

    basis <- basis(nmf_result)
    rownames(basis) <- left_names

    #4. Restore original working directory
    setwd(prev_wd)

    #5. Create 3D Network from NMF and PLSR results

    # Create dataframe storing the max cluster of each
    # basis ("left") element (and the corresponding 
    # coefficient)
    left_cluster_raw <- apply(basis, 1, which.max) #in NMF, max coef determines cluster
    left_max <- apply(basis, 1, max)
    left_cluster <- data.frame(
        left_names = names(left_cluster_raw), 
        cluster = left_cluster_raw,
        left_weight = left_max
    )

    # Create dataframe storing the max cluster of each
    # loading ("common") element (and the corresponding
    # coefficient)
    common_cluster_raw <- apply(loading, 1, which.max) #in NMF, max coef determines cluster
    common_max <- apply(loading, 1, max)
    common_cluster <- data.frame(
        common_names = names(common_cluster_raw), 
        cluster = common_cluster_raw,
        common_weight = common_max
    )

    # Merge loading coefficients and basis coefficients by clusters
    merged_left_common <- merge(common_cluster, left_cluster, by="cluster")

    # Merge left and common with p_list ("right")
    merged_left_common_right <- merge(merged_left_common, p_list, by="left_names")

    # 5. Return 3D Network
    list(
        plsr=plsr_result,
        nmf=list(
            nmf_result=nmf_result,
            basis=data.frame(basis),#data.frame(basis(nmf_result)),
            loading=data.frame(loading_original)#data.frame(coef(nmf_result))
        ),
        merged_final=merged_left_common_right
    )
}

write_all_to_disk <- function(name, object, folder) {
    if (class(object) == "data.frame") {
        final_path = file.path(folder, paste(name,".csv",sep=""))
        write.csv(
            object,
            final_path
        )
    }
    else if (class(object) == "list") {
        new_path <- file.path(folder, name)
        dir.create(new_path, showWarnings=FALSE)
        lapply(
            1:length(object),
            function(idx) {
                object_name <- names(object)[idx]
                if (length(object_name) > 0) {
                    write_all_to_disk(
                         object_name,
                         object[[idx]],
                         new_path
                    )
                }
            }
        )
    }
}

# need to reference iPhDNet to explain
#     "generalized into model3D"
# how to formulate problems
# e.g. metabolites, peptides
# important to understand in order to relate matrix files to package inputs

# left = output
# right = input
# "common"

#' Create a 3D network model.
#'
#' This function creates a 3D network model relating two matrices.
#' The columns of each matrix should represent two sets of objects,
#' with the rows of both matrices corresponding to observations under
#' the same circumstances. For example, the columns might be 
#' phosphoproteins and histones, and the rows might correspond to drugs.
#' In this case, the data values would represent the transcription
#' or GCP activation level in cells treated with different drugs.
#' 
#' The final 3D network model consists of connections between the rows
#' and both sets of columns, as follows:
#'
#'   - Connections between the rows and the left columns. These are found using 
#'     NMF (nonnegative matrix factorization). The data can contain negative
#'     values, as an antilog transformation is used internally. The number of
#'     NMF clusters can be specified, but otherwise it will be chosen
#'     automatically based on the value that gives the highest KL-index
#'     when using Ward hierarchical clustering.
#'     The connection weight between a given left column and row is the 
#'     maximum basis value for the left column object.
#'     See the function generate_nmf_modules for further information.
#' 
#'   - Connections between the rows of the right and left columns.
#'     These are found using PLSR. The connection strength is the absolute
#'     value of the PLSR regression coefficient. See the function
#'     generate_plsr for more information.
#'
#' @param left Matrix of left objects
#' @param right Matrix of right objects
#' @param nmf_nrun Number of iterations to use for NMF
#' @param k Number of clusters to use for NMF. This can be specified as a 
#'   - single number, or
#'   - consecutive range.
#' If unspecified, k is picked from 1 to min(num_cols - 1, num_rows - 1)
#' @param k_picker Method for picking k. If unspecified, k values are compared
#' using the KL-index of Ward clusterings based on euclidean distance.
#' @param seed Seed to use for NMF.
#' @param p_val_threshold Threshold for significant p-values in PLSR. 
#' @param out_folder Folder used for outputting results.
#' @param verbose Whether to print output (default TRUE).
#' 
#' @export
gen3DNet <- function(
    left,
    right,
    nmf_nrun=10,
    k=NULL, 
    k_picker=max_ward_kl,
    seed=123456,
    p_val_threshold=.0001,
    out_folder=NULL,
    verbose=TRUE
) {
    #Read input files
    if (is.character(left)) {
        left=read.csv(left, row.names=1)
    }
    if (is.character(right)) {
        right=read.csv(right, row.names=1)
    }
    if (is.null(out_folder)) {
        out_folder <- paste("gen3DNet", strptime(Sys.time(),"%Y-%m-%d %H:%M:%S"))
    }
    #Generate 3D model
    result <- create_gen3DNet(
        left=left,
        right=right,
        nmf_nrun=nmf_nrun,
        k=k,
        k_picker=k_picker,
        seed=seed,
        p_val_threshold=p_val_threshold,
        out_folder=out_folder
    )
    #Write output files
    if (verbose) {
        cli::cli_alert_success(paste("Writing to", out_folder))
    }
    write_all_to_disk(
        out_folder,
        result,
        "."
    )
    #Return same output data
    result
}
