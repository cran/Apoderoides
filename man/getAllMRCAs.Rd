\name{getAllMRCAs}
\alias{getAllMRCAs}
\title{getAllMRCAs}
\usage{
getAllMRCAs(tree,OTUrankData=NULL)
}
\description{Calculate all the most recent common ancestors (MRCAs) of ranks in the tree. Unlike getMRCA() in ape package, this function returns a tip node number when the rank is monotypic.}
\arguments{
\item{tree}{A phylogenetic tree to be checked. This is loaded by ape::read.tree() from a file.}
\item{OTUrankData}{A list composed of two character vectors. The first vector is tips of tree. The second vector is the upper rank of the tips. When this is NULL, the function assumes that all the tree tips are expressed as Genus_species like Homo_sapience, and calculate the MRCAs for genera. When this is not NULL, the function returns MRCAs based on the upper rank in this list.}
}
\value{
A list containing vectors of an MRCA node number.
}
\examples{
data(testTree)
#calculate MRCAs for genus
getAllMRCAs(testTree)
data(testRankList)
#calculate MRCAs for the rank in the list
getAllMRCAs(testTree,testRankList)
}
