#get upperRank from OTU name. List is necessary if getting rank upper than genus.
#this function is vectorized for the lowerRank.
get.upperRank<-function(data,OTUrankData=NULL)
{
  #when getting genus from species name
  if(is.null(OTUrankData))
  {
    match<-regexpr("_",data)
    result<-substr(data,1,match-1)
  }
  else
    #when getting upper rank from lower rank
  {
    match<-match(data,OTUrankData[[1]])
    result<-OTUrankData[[2]][match]
  }
  return(result)
}

deleteAnomaly<-function(tree,score,OTUrankData=NULL,drop=FALSE)
{
  droppingIndex<-list()
  score<-score[order(as.numeric(score[,2]),decreasing=TRUE),]
  topOTUScore<-score[1,2][[1]]
  #when the top OTU score is 0, return
  if(topOTUScore==0)
  {
    return("No anomaly OTU in the tree.")
  }
  #find index whose clade score is identical to the top clade score
  topScoreOTU<-score[topOTUScore==score[,2],1]
  if(length(topScoreOTU)==1)
  {
    #check if the topScoreOTU is monotypic
    topRank<-get.upperRank(topScoreOTU,OTUrankData)
    if(is.null(OTUrankData))
    {
      OTUrankData<-vector("list",2)
      OTUrankData[[2]]<-get.upperRank(tree$tip)
    }
    rankOTU<-extractOTUbyRankName_C(topRank,tree$tip,OTUrankData[[2]])
    #when the dropping OTU rank has one or more than two OTUs, drop the OTU with the highest score
    if(length(rankOTU)==1||length(rankOTU)>2)
    {
      return(list(score[1,1][[1]],drop.tip(tree,score[1,1][[1]])))
    }
    #when the dropping OTU rank has two OTUs, compare thier intruder score and drop one with the higher score
    else
    {
      int1<-score[rankOTU[1]==score[,1],4]
      int2<-score[rankOTU[2]==score[,1],4]
      if(int1>int2)
      {
        if(drop)
        {
          return(list(rankOTU[1],drop.tip(tree,rankOTU[1])))
        }
        else
        {
          return(list(rankOTU[1],tree))
        }
      }
      else if(int2>int1)
      {
        if(drop)
        {
          return(list(rankOTU[2],drop.tip(tree,rankOTU[2])))
        }
        else
        {
          return(list(rankOTU[2],tree))
        }
      }
      else
      {
        if(drop)
        {
          return(list(rankOTU,drop.tip(tree,rankOTU)))
        }
        else
        {
          return(list(rankOTU,tree))
        }
      }
    }
  }
  #when there are multiple OTUs with the highest score
  else
  {
    #find OTUs whose clade number has the smallest number of OTUs
    topScoreCladeNumber<-score[topOTUScore==score[,2],6]
    minCladeNumber<-Inf
    minCladeIndex<-0
    for(i in 1:length(topScoreCladeNumber))
    {
      #count the same clade number OTUs
      temp<-sum(topScoreCladeNumber[i]==topScoreCladeNumber)
      if(minCladeNumber>temp)
      {
        minCladeIndex<-i
        minCladeNumber<-temp
      }
      else if(minCladeNumber==temp)
      {
        minCladeIndex<-c(minCladeIndex,i)
      }
    }
    #cladeNumber with the smallest number of OTUs
    candidateCladeNumber<-topScoreCladeNumber[minCladeIndex]
    #choose ones with the smallest clade number
    candidateCladeNumber<-candidateCladeNumber[candidateCladeNumber==min(candidateCladeNumber)]
    candidateCladeNumber<-candidateCladeNumber[1]
    deletingOTU<-score[score[,6]==candidateCladeNumber,1]
    if(drop)
    {
      return(list(deletingOTU,drop.tip(tree,deletingOTU)))
    }
    else
    {
      return(list(deletingOTU,tree))
    }
  }
}

autoDeletion<-function(tree,OTUrankData=NULL,show_progress=TRUE,num_threads=1)
{
  if(length(tree$tip)<=3)
  {
    return("The tree includes only three or less OTUs.")
  }
  if(show_progress)
  {
    calcTime<-proc.time()
  }
  totalScore<-list()
  droppedOTUs<-character()
  dropIndex<-integer()
  allRankNames<-getAllRankNames(tree,OTUrankData)
  if(!is.null(OTUrankData))
  {
    rankList<-OTUrankData[[2]]
  }
  else
  {
    rankList<-get.upperRank(tree$tip)
  }
  allCentroids<-getAllCentroids(tree,OTUrankData,show_progress,num_threads)
  counter<-1
  progress = 0
  firstPositiveScoreOTU = -1
  currentPositiveScoreOTU = 0
  while(TRUE)
  {
    if(length(tree$tip)<=3)
    {
      break
    }
    score<-calc.Score(tree,OTUrankData,allRankNames,allCentroids,dropIndex,show_progress=show_progress,num_threads=num_threads)
    #check the score reached 0
    if(score[1,2][[1]]==0)
    {
      break
    }
    if(firstPositiveScoreOTU==-1)
    {
      firstPositiveScoreOTU<-sum(score[,2]>0)
    }
    currentPositiveScoreOTU<-sum(score[,2]>0)
    if(show_progress)
    {
      print(paste0("auto-deletion loop",counter))
      if(firstPositiveScoreOTU == -1)
      {
        print(paste0("progress: 0%"))
      }
      else if(firstPositiveScoreOTU == 0)
      {
        print(paste0("progress: 100%"))
      }
      else
      {
        print(paste0("progress: " ,100*(firstPositiveScoreOTU-currentPositiveScoreOTU)/firstPositiveScoreOTU,"%"))
      }
      counter<-counter+1
    }
    totalScore<-append(totalScore,list(score))
    temp<-deleteAnomaly(tree,score,OTUrankData)
    droppedOTUs<-append(droppedOTUs,as.vector(temp[[1]]))
    #index of dropped OTU
    index<-match(temp[[1]],tree$tip)
    dropIndex<-c(dropIndex,index)
    lastDropRank<-unique(get.upperRank(tree$tip[index[1]],OTUrankData))
    lastDropRankIndex<-match(lastDropRank,allRankNames)
    #renew centroids
    if(show_progress)
    {
      print("renewing a centroid")
      centroidTime<-proc.time()
    }
    allCentroids[[lastDropRankIndex]]<-getRankCentroid_C(lastDropRank,dropIndex,tree$tip,tree[[1]][,1],tree[[1]][,2],rankList,show_progress,num_threads)
    if(show_progress)
    {
      print(proc.time()-centroidTime)
    }
  }
  tree<-drop.tip(tree,dropIndex)
  if(show_progress)
  {
    print("total time")
    print(proc.time()-calcTime)
  }
  return(list(resultantTree=tree,droppedOTU=droppedOTUs,scoreTransition=totalScore))
}

getAllRankNames<-function(tree,OTUrankData=NULL)
{
  dataRank<-get.upperRank(tree$tip,OTUrankData)
  return(unique(dataRank))
}

getAllCentroids<-function(tree,OTUrankData=NULL,show_progress=FALSE,num_threads=1)
{
  allRankNames<-getAllRankNames(tree,OTUrankData)
  allCentroid<-vector("list",length(allRankNames))
  if(show_progress)
  {
    print("calculating all centroids")
    calcTime<-proc.time()
  }
  if(is.null(OTUrankData))
  {
    OTUrankData<-list()
    OTUrankData[[2]]<-get.upperRank(tree$tip)
  }
  allCentroid<-getAllCentroids_C(tree$tip,allRankNames,tree[[1]][,1],tree[[1]][,2],OTUrankData[[2]],show_progress,num_threads)
  if(show_progress)
  {
    print(proc.time()-calcTime)
  }
  return(allCentroid)
}

#obtain upper node of a given node/tip
findUpperNode<-function(tree,node)
{
  treeMatrix<-tree[[1]]
  return(treeMatrix[treeMatrix[,2]==node,1])
}

#return TRUE if all the given nodes belong to the same rank
is.monophyleticByRank<-function(tree,nodeIndex,OTUrankData)
{
  nodes<-tree$tip[nodeIndex]
  rank<-OTUrankData[[2]][nodeIndex]
  return(all(rank[1]==rank))
}

calc.Score<-function(tree,OTUrankData=NULL,allRankNames=NULL,allCentroids=NULL,dropIndex=NULL,sort=TRUE,show_progress=TRUE,num_threads=1)
{
  OTUList<-tree$tip
  intScore<-numeric(length(OTUList))
  outScore<-numeric(length(OTUList))
  OTUScore<-numeric(length(OTUList))
  #numeric group which has the identical score and is monophyletic
  cladeNumber<-numeric(length(OTUList))
  #total score for the monophyletic group
  cladeScore<-numeric(length(OTUList))
  #cladeScore/number of OTUs in the clade
  perCladeScore<-numeric(length(OTUList))
  if(is.null(allCentroids))
  {
    allCentroids<-getAllCentroids(tree,OTUrankData,show_progress,num_threads)
  }
  if(is.null(dropIndex))
  {
    dropIndex<-integer(0)
  }
  if(is.null(allRankNames))
  {
    allRankNames<-getAllRankNames(tree,OTUrankData)
  }
  if(is.null(OTUrankData))
  {
    OTUrankData<-list()
    OTUrankData[[1]]<-tree$tip
    OTUrankData[[2]]<-get.upperRank(tree$tip)
  }

  isScored<-logical(length(OTUList))
  scoreCounter<-1
  if(show_progress)
  {
    #for new line
    print("calculating score")
    calcTime<-proc.time()
    pb<-txtProgressBar(style=3)
  }
  range<-1:length(OTUList)
  range<-setdiff(range,dropIndex)
  for(i in range)
  {
    if(show_progress)
    {
      setTxtProgressBar(pb,i/length(range))
    }
    if(isScored[i])
    {
      next
    }
    isScored[i]<-TRUE
    cladeIndex<-i
    nextCladeIndex<-i
    checkingNode<-i
    while(is.monophyleticByRank(tree,nextCladeIndex,OTUrankData))
    {
      cladeIndex<-nextCladeIndex
      checkingNode<-findUpperNode(tree,checkingNode)
      #when the checking node is the root node
      if(length(checkingNode)==0)
      {
        break
      }
      nextCladeIndex<-findSubTips_C(tree$tip,tree[[1]][,1],tree[[1]][,2],checkingNode)
    }
    intscore<-calcIntScore_C(tree$tip,tree[[1]][,1],tree[[1]][,2],OTUList[i],allCentroids,allRankNames,OTUrankData[[2]])
    outscore<-calcOutScore_C(tree$tip,tree[[1]][,1],tree[[1]][,2],OTUList[i],allCentroids,allRankNames,OTUrankData[[2]],dropIndex=dropIndex)
    intScore[cladeIndex]<-intscore
    outScore[cladeIndex]<-outscore
    OTUScore[cladeIndex]<-intscore+outscore
    perCladeScore[cladeIndex]<-(intscore+outscore)/length(cladeIndex)
    isScored[cladeIndex]<-TRUE
    cladeNumber[cladeIndex]<-scoreCounter
    scoreCounter<-scoreCounter+1
  }
  if(show_progress)
  {
    close(pb)
    print(proc.time()-calcTime)
  }
  score<-array(c(OTUList,perCladeScore,OTUScore,intScore,outScore,cladeNumber),dim=c(length(OTUList),6))
  score<-score[order(as.numeric(score[,6])),]
  colnames(score)<-c("OTU","perCladeOTUScore","sum","intruder","outlier","#clade")
  if(sort)
  {
    return(score[order(as.numeric(score[,2]),decreasing=TRUE),])
  }
  else
  {
    return(score)
  }
}
