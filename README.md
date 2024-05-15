# gen3DNet: An R Package for Generating 3D Network Models
## Abstract
### Motivation: 
Networks are ubiquitous to organize and represent data in the real world such as the biomedical networks, social networks, etc. Due to the attractive property of networks which can explicitly capture the rich relationships among data samples, there is an increasing demand to transform linkage-free data to graph-structure data for efficient downstream analyses. Existing works typically focus on representing a single data object like gene expression profile in a homogeneous 2D network, which fail to deal with situations where two different data objects are involved to create a 3D heterogeneous network.
### Results: 
In this paper, we introduce an R package, gen3DNet (a generic version of the iPhDNet), for generating 3D network models from two correlated objects with shared common factors. Specifically, gen3DNet builds the relationships between samples and shared factors using the non-negative matrix factorization, where three clustering techniques are evaluated to determine the number of functional modules. In addition, it builds the relationships between samples from the two distinct data objects based on a partial least squares regression model. Usage of the package is illustrated through a real-world application.
## How to use
You can install gen3DNet R package from CRAN using: 

`install.packages(c('generate_nmf_modules', 'generate_plsr', 'gen3DNet'))`

### Installation of other dependencies
* Install NMF (>= 0.23.0) using `install.packages('NMF')`

If you have any problems running our code, please feel free to contact us (smollah@wustl.edu)
## Citation
Please cite the following paper if you are using gen3DNet for your research:

gen3DNet: An R Package for Generating 3D Network Models. Paul Morrison, Charles Lu, Tina Tang, Shamim Mollah bioRxiv 

