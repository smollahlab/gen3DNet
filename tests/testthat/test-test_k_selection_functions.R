histon_path <- "../../histon_data.csv"
phospho_path <- "../../phospho_data.csv"

test_that("k_selection_functions work", {
  result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,
    k_picker=max_ward_kl
  )


result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,
    k_picker=max_silhouette_consensus
)
result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,

    k_picker=max_cophenetic
)
result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,

    k_picker=kneedle_silhouette_consensus
)

result <- gen3DNet(
    histon_path,
    phospho_path,
    nmf_nrun=10,
    k_picker=kneedle_cophenetic
)


})
