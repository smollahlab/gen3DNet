#' @importFrom grDevices dev.off pdf
#' @importFrom graphics layout lines par plot
#' @importFrom stats dt vcov
#' @importFrom utils read.csv write.csv
#' @importFrom NMF coef nmf basismap coefmap

dist2d <- function(a,b,c) {
 v1 <- b - c
 v2 <- a - b
 m <- cbind(v1,v2)
 d <- det(m)/sqrt(sum(v1*v1))
 d
}


#' Find the knee/elbow point in a vector using the Kneedle algorithm.
#'
#' This function uses the Kneedle algorithm (Satopaa 2011)
#' to find the index of the knee point in the provided vector.
#' If the values are mostly increasing, use sign = 1. If they are
#' mostly decreasing, use sign = -1.
#'
#' @param values The values to find a knee/elbow in.
#' @param sign -1 if the values are mostly decreasing, 1 if the 
#' values are mostly decreasing.
#'
#' @return The index of the knee/elbow.
#'
kneedle <- function(values, sign) {
    start = c(1, values[1])
    end = c(length(values), values[length(values)])
    k <- which.max(lapply(1:length(values),
                     function(idx) {
                         sign * -1 * dist2d(c(idx, values[idx]),
                                start,
                                end
                         )
                     })
              )
    #plot(1:length(values),values,type="l")
    #lines(c(1,length(values)),
    #      c(values[1],values[length(values)]),
    #      col="red")
    #points(k, values[k],
    #       col="blue")
    k
}

remove.na <- function(a) {
    a[!is.na(a)]
}




#' Pick k based on the maximum KL index of Ward clusterings
#'
#' This function uses the NbClust package (Charrad 2014)
#' to estimate the number of clusters present.
#'
#' First, the data is repeatedly clustered with 1 cluster, then 
#' 2, and so on, using the Ward hierarchical clustering algorithm
#' (Ward 1963). Afterwards, each set of clusters is scored
#' using the Krzanowski-Lai index (Krzanowski 1988). The
#' number of clusters with the highest index is then chosen. 
#'
#' @param data The data with an unknown number of clusters.
#' @param k_range The range of possible k values
#'
#' @return The chosen value of k.
#' 
#' @export
max_ward_kl <- function(data, k_range) {
  NbClust::NbClust(t(data), diss=NULL, distance = "euclidean", min.nc=min(k_range), max.nc=max(k_range),
          method = "ward.D2", index = "kl")$Best.nc["Number_clusters"]
}

#' Pick k based on the maximum silhouette score
#'
#' This function uses the average silhouette width (Rousseeuw 1987), 
#' as implemented in the NMF library (Gaujoux 2010). 
#'
#' @param data The data with an unknown number of clusters.
#' @param k_range The range of possible k values
#'
#' @return The chosen value of k.
#' 
#' @export
max_silhouette_consensus <- function(data, k_range) {
  k_range[which.max(nmf(normalize_nmf(data), k_range, nrun=10)$measures$silhouette.consensus)]
}

#' Pick k based on the maximum cophenetic score 
#' 
#' This function uses the cophenetic correlation (Lessig 1972),
#' which is implemented in the NMF library (Gaujoux 2010).
#'
#' @param data The data with an unknown number of clusters.
#' @param k_range The range of possible k values
#'
#' @return The chosen value of k.
#' 
#' @export
max_cophenetic <- function(data, k_range) {
  k_range[which.max(nmf(normalize_nmf(data), k_range, nrun=10)$measures$cophenetic)]
}

#' Pick k based on the knee point in the silhouette score
#'
#' This function uses the knee point (Satopaa 2011) of the average silhouette
#' width (Lessig 1972), which is implemented in the NMF library 
#' (Gaujoux 2010).
#'
#' @param data The data with an unknown number of clusters.
#' @param k_range The range of possible k values
#'
#' @return The chosen value of k.
#' 
#' @export
kneedle_silhouette_consensus <- function(data, k_range) {
  k_range[kneedle(remove.na(nmf(normalize_nmf(data), k_range, nrun=10)$measures$silhouette.consensus), 1)]
}

#' Pick k based on the knee point in the cophenetic correlation
#'
#' This function uses the knee point (Satopaa 2011) of the cophenetic 
#' correlation (Lessig 1972), which is implemented in the NMF library 
#' (Gaujoux 2010).
#'
#' @param data The data with an unknown number of clusters.
#' @param k_range The range of possible k values
#'
#' @return The chosen value of k.
#' 
#' @export
kneedle_cophenetic <- function(data, k_range) {
  k_range[kneedle(remove.na(nmf(normalize_nmf(data), k_range, nrun=10)$measures$cophenetic), 1)]
}

#' Create a general 3D network.
#' 
#' This function performs calculations to generate a 3D network.
#' It is not meant to be called directly. You should use the gen3DNet
#' function, which sanitizes arguments and reads in data from file paths.
#' See help(gen3DNet) for further documentation.
#' 
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
#' Possible values:
#'   - max_cophenetic
#'   - kneedle_silhouette_consensus
#'   - kneedle_cophenetic 
#'   - max_silhouette_consensus
#'   - max_cophenetic
#'   - max_ward_kl
#' To learn more, use help() for each of these functions.
#'
#' @param seed Seed to use for NMF.
#' @param p_val_threshold Threshold for significant p-values in PLSR.
#' @param out_folder Folder used for outputting results.
#' @param verbose Whether to print output (default TRUE).
#' 
#' @examples
#' 
#' library("gen3DNet")
#' histon_path <- system.file("extdata", "histon_data.csv", package="gen3DNet")
#' phospho_path <- system.file("extdata", "phospho_data.csv", package="gen3DNet")
#' result <- gen3DNet(
#'    histon_path,
#'    phospho_path,
#'    nmf_nrun = 10,
#'    p_val_threshold = 0.01, 
#     # Use one of the following k_selection functions 
#'    # k_picker = max_cophenetic
#'    # k_picker = kneedle_silhouette_consensus
#'    # k_picker = kneedle_cophenetic 
#'    # k_picker = max_silhouette_consensus
#'    # k_picker = max_cophenetic
#'    # k_picker = max_ward_kl
#'    k_picker = max_ward_kl
#' )
#'
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
           min_k <- 2            #Min k is 2  
           max_k <- min(dim(left) - 1) #Max k is num. cols or 30
       }
       else {
           #The user has provided a range of k values.
           #So first, ensure the k values are consecutive.

           #Select min and max
           min_k <- min(k)
           max_k <- max(k)
           if (!identical(min_k:max_k, k)) {
               stop("gen3DNet must consider every choice from min_k to max_k")
           }
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

    loading_original <- as.matrix(coef(nmf_result))
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
    
    # Create column storing the direction ("pos" or "neg") 
    # of left and common correlation based on values in left matrix 
    left_common_directions <- data.frame()
    for (i in 1:nrow(merged_left_common)) {
      cname = merged_left_common$common_names[i]
      lname = merged_left_common$left_names[i]
      if (left[cname,lname] < 0){
        direction = "neg"
      } 
      else {
        direction = "pos" 
      } 
      left_common_directions = append(left_common_directions, direction)
    }
    
    #Add left_common direction to results dataframe 
    merged_left_common$left_common_directions <- left_common_directions
    merged_left_common <- apply(merged_left_common,2,as.character)

    # Merge left and common with p_list ("right")
    merged_left_common_right <- merge(merged_left_common, p_list, by="left_names")

  # Create column storing the direction ("pos" or "neg") 
    # of right and common correlation based on values in right matrix 
    right_common_directions <- data.frame()
    for (i in 1:nrow(merged_left_common_right)) {
      cname = merged_left_common_right$common_names[i]
      rname = merged_left_common_right$right_names[i]
      if (right[cname,rname] < 0){
        direction = "neg"
      } 
      else {
        direction = "pos" 
      } 
      right_common_directions = append(right_common_directions, direction)
    }

    #Add right_common direction to results dataframe 
    merged_left_common_right$right_common_directions <- right_common_directions
    merged_left_common_right <- apply(merged_left_common_right,2,as.character)
    merged_left_common_right <- as.data.frame(merged_left_common_right)
    
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

#' Write a nested list as nested folders, with data.frame objects as csvs.
#'
#' @param name The name of the top-level folder to be written
#' @param object The nested list
#' @param folder The parent folder where the top-level folder should be placed.
#' 
#' @examples
#' # For example, if the object is
#' object = list(
#'     item1 = list(
#'         df1 = data.frame(a=c(1,2),b=c(3,4)),
#'         df2 = data.frame(a=c(1,2),b=c(3,4))
#'     ),
#'     item2 = list(
#'         df3=data.frame(a=c(1,2),b=c(3,4))
#'     ),
#'     df4 = data.frame(a=c(1,2),b=c(3,4))
#' )
#' # and the call is
#' write_all_to_disk("toplevelfolder", object, ".")
#' # The following files and folders would be written:
#' # toplevelfolder/
#' #     item1/
#' #         df1.csv
#' #         df2.csv
#' #     item2/
#' #         df3.csv
#' #     df4.csv
#' #
#'
#' @export
write_all_to_disk <- function(name, object, folder) {
    if (class(object) == "data.frame") {
        final_path = file.path(folder, paste(name,".csv",sep=""), fsep = .Platform$file.sep)
        write.csv(
            object,
            final_path
        )
    }
    else if (class(object) == "list") {
        new_path <- file.path(folder, name, fsep = .Platform$file.sep)
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

MIN_SIZE = 5 

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
#' Possible values:
#'   - max_cophenetic
#'   - kneedle_silhouette_consensus
#'   - kneedle_cophenetic 
#'   - max_silhouette_consensus
#'   - max_cophenetic
#'   - max_ward_kl
#' 
#' To learn more, use help() for each of these functions.
#'
#' @param seed Seed to use for NMF.
#' @param p_val_threshold Threshold for significant p-values in PLSR.
#' @param out_folder Folder used for outputting results.
#' @param verbose Whether to print output (default TRUE).
#' 
#' @examples
#' 
#' library("gen3DNet")
#' histon_path <- system.file("extdata", "histon_data.csv", package="gen3DNet")
#' phospho_path <- system.file("extdata", "phospho_data.csv", package="gen3DNet")
#' result <- gen3DNet(
#'    histon_path,
#'    phospho_path,
#'    nmf_nrun = 10,
#'    p_val_threshold = 0.01, 
#     # Use one of the following k_selection functions 
#'    # k_picker = max_cophenetic
#'    # k_picker = kneedle_silhouette_consensus
#'    # k_picker = kneedle_cophenetic 
#'    # k_picker = max_silhouette_consensus
#'    # k_picker = max_cophenetic
#'    # k_picker = max_ward_kl
#'    k_picker = max_ward_kl
#' )
#'
#' @export
#' 

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
    if (any(c(dim(left), dim(right)) < MIN_SIZE)) {
        stop("Please ensure that each provided dataframe has at least 3 rows and columns.")
    }
    if (is.null(out_folder)) {
        a = tempdir()
        out_folder <- paste(
            a, "gen3DNet", chartr(old=" ", new="_", chartr(old=":",new="-",strptime(Sys.time(),"%Y-%m-%d %H:%M:%S"))),
            sep="_"
        )
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
        ""
    )
    #Return same output data
    result
}
