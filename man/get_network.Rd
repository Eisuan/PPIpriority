% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions.R
\name{get_network}
\alias{get_network}
\title{Get PPI network from a PPIpriority_obj}
\usage{
get_network(object)
}
\arguments{
\item{object}{a PPIpriority object.}
}
\value{
An igraph object
}
\description{
Accessor to the igraph object stored in a PPIpriority obj.
}
\examples{
#Set extended timeout, not necessary while using the package
options(timeout=600)

#Build human PPI
network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
score_threshold=800, input_directory=getwd())

# Access network igraph
get_network(network_ppi_obj)
}
