#get MRCAs for all upper ranks in a tree.
#if a rank is monophyletic, return the tip node number instead of NULL.
getAllMRCAs<-function(tree,OTUrankData=NULL)
{
  allRankNames<-getAllRankNames(tree,OTUrankData)
  allTipUpperRank<-get.upperRank(tree$tip,OTUrankData)
  MRCAlist<-vector("list",length(allRankNames))
  for(i in 1:length(MRCAlist))
  {
    #use regular expressions to conduct exact matching
    rankTips<-grep(paste0("^",allRankNames[i],"$"),allTipUpperRank)
    temp<-getMRCA(tree,rankTips)
    if(is.null(temp))
    {
      MRCAlist[[i]]<-rankTips
    }
    else
    {
      MRCAlist[[i]]<-temp
    }
  }
  return(MRCAlist)
}

#get MRCA of rankName in tree
getRankMRCA<-function(rankName,tree,dropIndex,rankList)
{
  #use regular expressions to conduct exact matching
  rankTips<-grep(paste0("^",rankName,"$"),rankList)
  rankTips<-setdiff(rankTips,dropIndex)
  MRCA<-getMRCA(tree,rankTips)
  if(is.null(MRCA))
  {
    MRCA<-rankTips
  }
  return(MRCA)
}

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
#load a centroid or MRCA score and return dropping OTU(s)
deleteSubFunc<-function(tree,score,OTUrankData=NULL)
{
  topOTUScore<-score[1,2][[1]]
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
      return(score[1,1][[1]])
    }
    #when the dropping OTU rank has two OTUs, compare thier intruder score and drop one with the higher score
    else
    {
      int1<-score[rankOTU[1]==score[,1],4]
      int2<-score[rankOTU[2]==score[,1],4]
      if(int1>int2)
      {
        return(rankOTU[1])
      }
      else if(int2>int1)
      {
        return(rankOTU[2])
      }
      else
      {
        return(rankOTU)
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
    candidateCladeNumber<-as.integer(topScoreCladeNumber[minCladeIndex])
    #choose ones with the smallest clade number
    candidateCladeNumber<-candidateCladeNumber[candidateCladeNumber==min(candidateCladeNumber)]
    candidateCladeNumber<-candidateCladeNumber[1]
    deletingOTU<-score[score[,6]==candidateCladeNumber,1]
    return(deletingOTU)
  }
}

deleteAnomaly<-function(tree,scores,OTUrankData=NULL,drop=FALSE,prior="MRCA")
{
  droppingIndex<-list()
  scores[[1]]<-scores[[1]][order(as.numeric(scores[[1]][,2]),decreasing=T),]
  scores[[2]]<-scores[[2]][order(as.numeric(scores[[2]][,2]),decreasing=T),]

  topOTUScore<-max(scores[[1]][1,2][[1]],scores[[2]][1,2][[1]])
  #when the top OTU score is 0, return
  if(topOTUScore==0)
  {
    return("No anomaly OTU in the tree.")
  }
  #when centroidScore > MRCAScore, use centroid score
  topScore1<-as.numeric(scores[[1]][1,2][[1]])
  topScore2<-as.numeric(scores[[2]][1,2][[1]])
  if(topScore1>topScore2)
  {
    dropOTU<-deleteSubFunc(tree,scores[[1]],OTUrankData)
  }
  #when centroidScore < MRCAScore, use MRCA score
  else if(topScore1<topScore2)
  {
    dropOTU<-deleteSubFunc(tree,scores[[2]],OTUrankData)
  }
  #when centroidScore == MRCAScore, compare the number of dropping OTUs and use fewer score
  else
  {
    dropOTUcand1<-deleteSubFunc(tree,scores[[1]],OTUrankData)
    dropOTUcand2<-deleteSubFunc(tree,scores[[2]],OTUrankData)
    if(length(dropOTUcand1)<length(dropOTUcand2))
    {
      dropOTU<-dropOTUcand1
    }
    else if(length(dropOTUcand1)>length(dropOTUcand2))
    {
      dropOTU<-dropOTUcand2
    }
    #when the number of dropping OTUs are equal, use one defined by the prior argument
    else
    {
      if(prior=="centroid")
      {
        dropOTU<-dropOTUcand1
      }
      else if (prior=="MRCA")
      {
        dropOTU<-dropOTUcand2
      }
      else
      {
        print("'prior' argument must be 'MRCA' or 'centroid'.")
        return()
      }
    }
  }
  if(drop)
  {
    return(list(dropOTU,drop.tip(tree,dropOTU)))
  }
  else
  {
    return(list(dropOTU,tree))
  }
}

autoDeletion<-function(tree,OTUrankData=NULL,show_progress=TRUE,num_threads=1,prior="MRCA")
{
  if(prior!="MRCA"&&prior!="centroid")
  {
    print("'prior' argument must be 'MRCA' or 'centroid'.")
    return()
  }
  if(length(tree$tip)<=3)
  {
    return("The tree includes only three or less OTUs.")
  }
  if(show_progress)
  {
    calcTime<-proc.time()
  }
  totalCentroidScore<-list()
  totalMRCAScore<-list()
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
  allMRCAs<-getAllMRCAs(tree,OTUrankData)
  counter<-1
  progress = 0
  firstPositiveScoreOTU = -1
  currentPositiveScoreOTU = 0
  scores<-vector("list",2)
  while(TRUE)
  {
    if(length(tree$tip)<=3)
    {
      break
    }
    scores[[1]]<-calc.Score(tree,OTUrankData,allRankNames,allCentroids,dropIndex,show_progress=show_progress,num_threads=num_threads)[[1]]
    scores[[2]]<-calc.Score(tree,OTUrankData,allRankNames,allMRCAs,dropIndex,show_progress=show_progress,num_threads=num_threads)[[1]]
    #check the score reached 0
    #somewhat if(scores[[1]][1,2][[1]]==0||scores[[2]][1,2][[1]]==0) gives an error
    if(scores[[1]][1,2][[1]]==0)
    {
      break
    }
    if(scores[[2]][1,2][[1]]==0)
    {
      break
    }
    if(firstPositiveScoreOTU==-1)
    {
      firstPositiveScoreOTU<-max(sum(scores[[1]][,2]>0),sum(scores[[2]][,2]>0))
    }
    currentPositiveScoreOTU<-max(sum(scores[[1]][,2]>0),sum(scores[[2]][,2]>0))
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
    totalCentroidScore<-append(totalCentroidScore,list(scores[[1]]))
    totalMRCAScore<-append(totalMRCAScore,list(scores[[2]]))
    temp<-deleteAnomaly(tree,scores,OTUrankData)
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
    allMRCAs[[lastDropRankIndex]]<-getRankMRCA(lastDropRank,tree,dropIndex,rankList)
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
  return(list(resultantTree=tree,droppedOTU=droppedOTUs,centroidScoreTransition=totalCentroidScore,MRCAScoreTransition=totalMRCAScore))
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

calc.Score<-function(tree,OTUrankData=NULL,allRankNames=NULL,allCores=NULL,dropIndex=NULL,sort=TRUE,show_progress=TRUE,num_threads=1)
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
  if(is.null(allCores))
  {
    scores<-vector("list",2)
    names(scores)<-c("CentroidScore","MRCAScore")
    allCentroids<-getAllCentroids(tree,OTUrankData,show_progress,num_threads)
    allMRCAs<-getAllMRCAs(tree,OTUrankData)
    allCores<-list(allCentroids,allMRCAs)
  }
  else
  {
    allCores<-list(allCores)
    scores<-vector("list",1)
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
  range<-1:length(OTUList)
  range<-setdiff(range,dropIndex)
  for(j in 1:length(scores))
  {
    isScored<-logical(length(OTUList))
    scoreCounter<-1
    if(show_progress)
    {
      #for new line
      print("calculating score")
      calcTime<-proc.time()
      pb<-txtProgressBar(style=3)
    }
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
      intscore<-calcIntScore_C(tree$tip,tree[[1]][,1],tree[[1]][,2],OTUList[i],allCores[[j]],allRankNames,OTUrankData[[2]])
      outscore<-calcOutScore_C(tree$tip,tree[[1]][,1],tree[[1]][,2],OTUList[i],allCores[[j]],allRankNames,OTUrankData[[2]],dropIndex=dropIndex)
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
    scores[[j]]<-array(c(OTUList,perCladeScore,OTUScore,intScore,outScore,cladeNumber),dim=c(length(OTUList),6))
    scores[[j]]<-scores[[j]][order(as.numeric(scores[[j]][,6])),]
    colnames(scores[[j]])<-c("OTU","perCladeOTUScore","sum","intruder","outlier","#clade")
    if(sort)
    {
      scores[[j]]<-scores[[j]][order(as.numeric(scores[[j]][,2]),decreasing=T),]
    }
  }
  return(scores)
}
