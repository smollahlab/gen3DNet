% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_nmf_modules.R
\name{generate_nmf_modules}
\alias{generate_nmf_modules}
\title{Generate NMF Modules}
\usage{
generate_nmf_modules(
  left_data,
  nmf_nrun,
  k_range,
  k_picker = gen3DNet::max_ward_kl,
  seed,
  verbose = FALSE
)
}
\arguments{
\item{left_data}{This is the data that will be clustered using NMF. Signatures will be generated for the object type that is represented by rows.}

\item{nmf_nrun}{The number of r}

\item{k_range}{Either a) A consecutive range of possible k values, b) a single k value, or c) nothing, meaning that generate_nmf_modules will assume the range is from 2 to the number of rows.}

\item{k_picker}{Any function that takes a dataframe and a consecutive range of potential k values}

\item{seed}{The seed to use with NMF}

\item{verbose}{Whether to print output.}
}
\description{
This function uses NMF to generate functional modules representing interactions
between the entities in the rows and columns of the dataframe. This has two
components:
Each entity is assigned to a cluster
"Column signatures" are generated for each row object. For example, if rows are histones and columns are drugs, this function will generate "histone signatures" that express
Because NMF relies on knowing the number of clusters (k) in advance, this function also determines this automatically.
Currently, the default option uses the number of clusters that produces the lowest KL index with Ward clustering. The package depends on NbClust to implement this feature. However, the user can override this option with a custom function.
}
