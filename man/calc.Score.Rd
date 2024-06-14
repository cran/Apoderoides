\name{calc.Score}
\alias{calc.Score}
\title{calc.Score}
\usage{
calc.Score(tree,OTUrankData=NULL,
allRankNames=NULL,allCores=NULL,dropIndex=NULL,
sort=TRUE,show_progress=TRUE,num_threads=1)
}
\description{
Calculate scores of a phylogenetic tree to find and prioritize erroneous taxa to delete.
}
\arguments{
\item{tree}{A phylogenetic tree to be checked. This is loaded by ape::read.tree() from a file.}
\item{OTUrankData}{A list composed of two character vectors. The first vector is tips of tree. The second vector is the upper rank of the tips. When this is NULL, the function assumes that all the tree tips are expressed as Genus_species like Homo_sapience, and calculate the score for genera. When this is not NULL, the function returns scores based on the upper rank in this list.}
\item{allRankNames}{This can be omitted.
This is a unique character vector of the upper ranks of the tree tips. If given, the calculation will be a little faster.}
\item{allCores}{This can be omitted.
This is a list of numeric vectors of the centroids or MRCAs of ranks. If given, the calculation will be a little faster.}
\item{dropIndex}{This can be omitted.
A numeric vector of indices of tree tips. The tree tips indicated by this dropIndex will be removed from the score calculation.}
\item{sort}{If TRUE, the calculation result is sorted by descending order of the total score.}
\item{show_progress}{If TRUE, calculation progress is shown on the R console.}
\item{num_threads}{A positive integer to specify the number of threads to calculate the scores.}
}
\value{
A list containing two matrices of characters. The first one is the score based on the centroids, and the second is that based on the MRCAs. The following explains the columns in the matrix.
\item{OTU}{The name of tree tip.}
\item{perCladeOTUScore}{The final score calculated by "sum" divided by the number of OTUs with the same "#clade".}
\item{sum}{The sum of "intruder" and "outlier" for the OTU.}
\item{intruder}{The intruder score showing how many ranks the OTU intruding into.}
\item{outlier}{The outlier score showing how the OTU is far away from the core clade of the belonging rank.}
\item{#clade}{The clade number. Monophyletic OTUs with the same rank has the same #clade.}
}
\examples{
data(testTree)
#calculate scores for genus
calc.Score(testTree)
data(testRankList)
#calculate scores for the rank in the list
calc.Score(testTree,testRankList)
}
