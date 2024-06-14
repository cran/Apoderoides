\name{getAllCentroids}
\alias{getAllCentroids}
\title{getAllCentroids}
\usage{
getAllCentroids(tree,OTUrankData=NULL,show_progress=FALSE,num_threads=1)
}
\description{Calculate all the centroids of ranks in the tree. The centroid of a rank is equivalent to S-centroid by Slater (1978).}
\arguments{
\item{tree}{A phylogenetic tree to be checked. This is loaded by ape::read.tree() from a file.}
\item{OTUrankData}{A list composed of two character vectors. The first vector is tips of tree. The second vector is the upper rank of the tips. When this is NULL, the function assumes that all the tree tips are expressed as Genus_species like Homo_sapience, and calculate the centroids for genera. When this is not NULL, the function returns centroids based on the upper rank in this list.}
\item{show_progress}{If TRUE, calculation progress is shown on the R console.}
\item{num_threads}{A positive integer to specify the number of threads to calculate the scores.}
}
\value{
A list containing vectors of integers of centroid node number(s).
}
\examples{
data(testTree)
#calculate centroids for genus
getAllCentroids(testTree)
data(testRankList)
#calculate centroids for the rank in the list
getAllCentroids(testTree,testRankList)
}
\references{
Slater P. J. 1978. Centers to centroids in graphs. Journal of Graph Theory 2: 209–222.
}
