\name{deleteAnomaly}
\alias{deleteAnomaly}
\title{deleteAnomaly}
\usage{
deleteAnomaly(tree,scores,OTUrankData=NULL,drop=FALSE,prior="MRCA")
}
\description{
Delete tip(s) with the highest score from a tree.
}
\arguments{
\item{tree}{A phylogenetic tree to be checked. This is loaded by ape::read.tree() from a file.}
\item{scores}{A list of scores calculated by calc.Score function.}
\item{OTUrankData}{A list composed of two character vectors. The first vector is tips of tree. The second vector is the upper rank of the tips. When this is NULL, the function assumes that all the tree tips are expressed as Genus_species like Homo_sapience and that the score is calculated based on genera. When this is not NULL, the function assumes the score is calculated based on the upper rank in this list.}
\item{drop}{Whether the dropped OTU(s) is included in the returned tree.}
\item{prior}{Used only when the length of "scores" is two. "MRCA" or "centroid". This argument defines the prioritized score when scores based on MRCA and centroid are equal.}
}
\value{
A list of the length two. The first element is a vector of characters of deleted tip label(s). The second is a list of a phylogenetic tree without the deleted tip(s).
}
\examples{
data(testTree)
data(testRankList)
#calculate scores for the rank in the list
score<-calc.Score(testTree,testRankList)
#delete tip with the highest score from tree
deleteAnomaly(testTree,score,testRankList)
}
