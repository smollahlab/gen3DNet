# gen3DNet: An R Package for Generating 3D Network Models
## Abstract
### Motivation: 
Networks are ubiquitous to organize and represent data in the real world such as the biomedical networks, social networks, etc. Due to the attractive property of networks which can explicitly capture the rich relationships among data samples, there is an increasing demand to transform linkage-free data to graph-structure data for efficient downstream analyses. Existing works typically focus on representing a single data object like gene expression profile in a homogeneous 2D network, which fail to deal with situations where two different data objects are involved to create a 3D heterogeneous network.
### Results: 
We introduce an R package, gen3DNet (a generic version of the iPhDNet [1]), for generating 3D network models from two correlated objects with shared common factors. Specifically, gen3DNet builds the relationships between samples and shared factors using the non-negative matrix factorization, where three clustering techniques are evaluated to determine the number of functional modules. In addition, it builds the relationships between samples from the two distinct data objects based on a partial least squares regression model. Usage of the package is illustrated through a real-world application (see gen3DNet paper below). 
## How to use
You can install gen3DNet R package from CRAN using: 
```
install.packages('gen3DNet')
```

Or directly from Github using:
```
install.packages("remotes")
remotes::install_github("MollahLab/gen3DNet")
```

### Installation of other dependencies
* Install NMF (>= 0.23.0) using `install.packages('NMF')`

If you have any problems running our code, please feel free to contact us (smollah@wustl.edu)

### Example of usage
```
library("gen3DNet")
left <- system.file("extdata", "histon_data.csv", package="gen3DNet")
right <- system.file("extdata", "phospho_data.csv", package="gen3DNet")
result <- gen3DNet(
   left,
   right,
   nmf_nrun = 10,
   p_val_threshold = 0.01, 
   # k_picker = max_cophenetic
   # k_picker = kneedle_silhouette_consensus
   # k_picker = kneedle_cophenetic 
   # k_picker = max_silhouette_consensus
   # k_picker = max_cophenetic
   # k_picker = max_ward_kl
   k_picker = max_ward_kl
)
```
## Citation
Please cite the following paper if you are using gen3DNet for your research:

1. S. A. Mollah and S. Subramaniam, “Histone Signatures Predict Therapeutic Efficacy in Breast Cancer”. 2020, IEEE Open Journal of Engineering in Medicine and Biology, vol. 1, pp. 74-82.  (doi.org/10.1109/OJEMB.2020.2967105)

2. gen3DNet: An R Package for Generating 3D Network Models. Paul Morrison, Tina Tang, Charles Lu, Shamim Mollah bioRxiv 

