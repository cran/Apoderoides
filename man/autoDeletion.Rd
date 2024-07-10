\name{autoDeletion}
\alias{autoDeletion}
\title{autoDeletion}
\usage{
autoDeletion(
tree,OTUrankData=NULL,
show_progress=TRUE,num_threads=1,
prior="MRCA",criteria="composite"
)
}
\description{
Iterate calc.Score() and deleteAnomaly() until all the tree tips have 0 score or the number of the tips becomes three or lower.
}
\arguments{
\item{tree}{A phylogenetic tree to be checked. This is loaded by ape::read.tree() from a file.}
\item{OTUrankData}{A list composed of two character vectors. The first vector is tips of tree. The second vector is the upper rank of the tips. When this is NULL, the function assumes that all the tree tips are expressed as Genus_species like Homo_sapience, and calculates for genera. When this is not NULL, the function calculates based on the upper rank in this list.}
\item{show_progress}{If TRUE, calculation progress is shown on the R console.}
\item{num_threads}{A positive integer to specify the number of threads to calculate.}
\item{prior}{Used only when "criteria" is "both". "MRCA" or "centroid". This argument defines the prioritized score when scores based on MRCA and centroid are equal.}
\item{criteria}{Criteria nodes to calculate the scores."composite", "both", "MRCA" or "centroid". "MRCA" and "centroid" use their corresponding node to calculate both intruder and outlier scores. "composite" calculates intruder scores using MRCA and outlier ones using centroid, which is empirically known to be most effective. "both" calculates both of MRCA-based and centroid-based scores and uses the highest one to select taxa to be deleted.}
}
\value{
A list of the length three or four. The first element is a list of phylogenetic tree from which erroneous taxa are deleted. The second is a character vector of deleted taxa. The third and fourth are a list of lists showing the transition of the scores. When criteria is "both", third and fourth elements correspond to scores based on MRCA and centroid, respectively. See \link{calc.Score} about the contents of the third and fourth elements.
}
\examples{
data(testTree)
data(testRankList)
#calculate scores for the rank in the list, and delete all the erroneous tips
#this takes tens of seconds for calculation
result<-autoDeletion(testTree,testRankList)
#tree without erroneos tips
result[[1]]
#deleted tips
result[[2]]
#scores during iteration of score calculation and tip deletion
result[[3]]
}
