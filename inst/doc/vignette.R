## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("Apoderoides")

## ----results="hide", warning=FALSE, message=FALSE-----------------------------
library(Apoderoides)

## ----eval=FALSE---------------------------------------------------------------
#  data(testTree)

## ----eval=FALSE---------------------------------------------------------------
#  testTree <- read.tree(file="/directory/yourTree.tre")

## ----eval=FALSE---------------------------------------------------------------
#  #for the test tree
#  calc.Score(testTree)
#  #for the user imported tree
#  calc.Score(tree)

## -----------------------------------------------------------------------------
calc.Score(testTree,show_progress=FALSE)[[1]][1:10,]

## ----eval=FALSE---------------------------------------------------------------
#  data("testRankList")

## -----------------------------------------------------------------------------
data("testRankList")
testRankList[[1]][1:10]
testRankList[[2]][1:10]

## ----eval=FALSE---------------------------------------------------------------
#  calc.Score(testTree,testRankList)

## ----eval=FALSE---------------------------------------------------------------
#  #for genus level
#  autoDeletion(testTree)
#  #for family level
#  autoDeletion(testTree,testRankList)

