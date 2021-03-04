library(gen3DNet)

histon_path <- system.file("extdata", "histon_data.csv", package="gen3DNet")
phospho_path <- system.file("extdata", "phospho_data.csv", package="gen3DNet")

if (dir.exists("test_output")) {
    unlink("test_output",recursive=TRUE)
}

verify_folders_same <- function(folder_a, folder_b, parent_folder = ".") {
   subfolders_a = list.files(folder_a)
   subfolders_b = list.files(folder_b)
   if (length(subfolders_a) != length(subfolders_b)) {
       result = FALSE
   } 
   else {
       result = all(sapply(
           1:length(subfolders_a),
           function(idx) {
               item_a = subfolders_a[idx]
               item_b = subfolders_b[idx]
               n_folders = sum(dir.exists(item_a), dir.exists(item_b))
               #If both are files, not folders, check whether the files are equal.  
               if (n_folders == 0) {
                   sub_result <- verify_files_same(file.path(folder_a, item_a), file.path(folder_b, item_b))
               }
               #If one is a file, but the other is a folder, they can't be the same. 
               else if (n_folders == 1) {
                   sub_result <- FALSE
               }
               #If both are folders, recurse.
               else if (n_folders == 2) {
                   sub_result <- verify_folders_same(file.path(folder_a, item_a[idx]), file.path(folder_b, item_b[idx]))
               }
               sub_result
           }
       ))
   }
   result
}

verify_files_same <- function(file_a, file_b) {
    file_a_lines <- readLines(file_a)
    file_b_lines <- readLines(file_b)
    if (length(file_a_lines) == length(file_b_lines)) {
        result <- all(file_a_lines == file_b_lines)
    }
    else {
        result <- FALSE
    }
    result
}




test_that("Results match iPhDNet", {

result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,
    p_val_threshold=0.01, 
    #k_picker=max_cophenetic
    #k_picker=kneedle_silhouette_consensus
    #k_picker=kneedle_cophenetic 
    #k_picker=max_silhouette_consensus
    #k_picker=max_cophenetic
    k_picker=max_ward_kl,
    out_folder = "test_output"
)

expect(verify_folders_same("test_output","correct_output"))

})

test_that("k_selection_functions work", {

k_selection_functions = list(max_cophenetic, kneedle_silhouette_consensus, kneedle_cophenetic, max_silhouette_consensus, max_cophenetic, max_ward_kl)

for (k_selection_function in k_selection_functions) {

result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=3,
    p_val_threshold=0.01,
    #k_picker=max_cophenetic
    #k_picker=kneedle_silhouette_consensus
    #k_picker=kneedle_cophenetic
    #k_picker=max_silhouette_consensus
    #k_picker=max_cophenetic
    k_picker=k_selection_function
)

}

})


