\name{get.upperRank}
\alias{get.upperRank}
\title{get.upperRank}
\usage{
get.upperRank(data,OTUrankData=NULL)
}
\description{
Obtain upper rank of scientific names in data. When OTUrankData is not provided, this function returns genus names assuming the elements in data are scientific names connected by underlines like "Homo_sapiens". When OTUrankData is provided, this function searches data in OTUrankData[[1]] and returns OTUrankData[[2]] of the corresponding index.
}
\arguments{
\item{data}{A vector of characters.}
\item{OTUrankData}{A list composed of two character vectors. The first vector is tips of tree. The second vector is the upper rank of the tips. When this is NULL, the function assumes that all the tree tips are expressed as Genus_species like Homo_sapience, and calculate the score for genera. When this is not NULL, the function returns scores based on the upper rank in this list.}
}
\value{
A vector of characters of upper rank.
}
\examples{
#obtain genus name
get.upperRank(c("Oxalis_nipponica","Homo_sapiens"))
data(testTree)
data(testRankList)
#obtain higher rank names
get.upperRank(testTree$tip[1:3],testRankList)
}
