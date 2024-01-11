#include <Rcpp.h>
#include <omp.h>
#include <progress.hpp>
#include <vector>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::export]]
std::vector<int> getTipNo_C (std::vector<std::string> treeTip, std::vector<std::string> OTU)
{
  std::vector<int> result(OTU.size());
  for(unsigned int i = 0, size=OTU.size(); i < size; ++i)
  {
    result[i] = std::distance(treeTip.begin(),std::find(treeTip.begin(),treeTip.end(),OTU[i]));
    result[i]++;
  }
  return(result);
}

inline int getRootNode_C (std::vector<std::string> treeTip)
{
  return(treeTip.size()+1);
}

std::vector<int> tip2Root_C (std::vector<std::string> treeTip, std::vector<int> treeMatCol0,std::vector<int> treeMatCol1, int tip)
{
  int rootNode = getRootNode_C(treeTip);
  std::vector<int> path(1);
  path[0] = tip;
  //reserve the size assuming that the tree is bifurcating
  path.reserve((int)log2(treeTip.size()));
  int upperNode = tip;
  int lowerNode = tip;
  while(upperNode!=rootNode)
  {
    int index = std::distance(treeMatCol1.begin(),std::find(treeMatCol1.begin(),treeMatCol1.end(),lowerNode));
    upperNode = treeMatCol0[index];
    path.push_back(upperNode);
    lowerNode = upperNode;
  }
  return(path);
}

std::vector<int> intersect_C (std::vector<int> vec1,std::vector<int> vec2)
{
  std::sort(vec1.begin(),vec1.end());
  std::sort(vec2.begin(),vec2.end());
  std::vector<int> shared(vec1.size()+vec2.size());
  std::set_intersection(vec1.begin(),vec1.end(),vec2.begin(),vec2.end(),shared.begin());
  //remove extra 0 element
  std::vector<int>::iterator it = std::remove(shared.begin(),shared.end(),0);
  shared.erase(it,shared.end());
  return(shared);
}

//generic function to connect two vectors (somewhat insert does not work in Rcpp)
std::vector<int> connect_C (std::vector<int> vec1,std::vector<int> vec2)
{
  std::vector<int> result(vec1.size()+vec2.size());
  //declare i here to use i in the next for loop
  unsigned int i = 0;
  for(unsigned size = vec1.size(); i < size; ++i)
  {
    result[i] = vec1[i];
  }
  for(unsigned int j = 0, size2 = vec2.size(); j < size2; ++j)
  {
    result[i+j] = vec2[j];
  }
  return(result);
}
//find intersect of two nodepaths to the root node without sorting them
std::vector<int> rootPathIntersect(std::vector<int> vec1,std::vector<int> vec2)
{
  std::vector<int> result;
  result.reserve(std::min(vec1.size(),vec2.size()));
  for(unsigned int i = 0, size = vec1.size(); i < size; ++i)
  {
    if(vec2.size() - 1 -i < 0)
    {
      break;
    }
    //check vec1 and vec2 from the last element (root node)
    if(vec1[vec1.size() - 1 - i]==vec2[vec2.size() -1 - i])
    {
      result.push_back(vec1[vec1.size() - 1 - i]);
    }
    else
    {
      break;
    }
  }
  std::reverse(result.begin(),result.end());
  return(result);
}

std::vector<int> findNodePath_C (std::vector<std::string> treeTip, std::vector<int> treeMatCol0,std::vector<int> treeMatCol1, int startTip, int goalTip)
{
  if(startTip==goalTip)
  {
    return(std::vector<int>{startTip});
  }
  std::vector<int> start2Root = tip2Root_C(treeTip,treeMatCol0, treeMatCol1,startTip);
  std::vector<int> goal2Root = tip2Root_C(treeTip,treeMatCol0, treeMatCol1,goalTip);
  std::vector<int> shared = rootPathIntersect(start2Root,goal2Root);
  std::vector<int> result;
  //when the two path only share the root node
  if(shared.size()==1)
  {
    std::vector<int> head(start2Root);
    head.pop_back();
    std::vector<int> tail(goal2Root);
    std::reverse(tail.begin(),tail.end());
    result = connect_C(head,tail);
    return(result);
  }
  else
  {
    int MRCA = shared[0];
    int startPathMRCAIndex = std::distance(start2Root.begin(),std::find(start2Root.begin(),start2Root.end(),MRCA));
    int goalPathMRCAIndex = std::distance(goal2Root.begin(),std::find(goal2Root.begin(),goal2Root.end(),MRCA));
    std::vector<int> head(start2Root.begin(),start2Root.begin()+startPathMRCAIndex);
    std::vector<int> tail(goal2Root.begin(), goal2Root.begin()+goalPathMRCAIndex+1);
    std::reverse(tail.begin(),tail.end());
    result = connect_C(head,tail);
    return(result);
  }
}

int findMRCA_C (std::vector<std::string> treeTip, std::vector<int> treeMatCol0,std::vector<int> treeMatCol1, std::vector<std::string> tipNames)
{
  std::vector<int> tips = getTipNo_C(treeTip, tipNames);
  int rootNode = getRootNode_C(treeTip);
  std::vector<int> shared;
  for(unsigned int i = 0, size = tips.size(); i < size; ++i)
  {
    if(i==0)
    {
      shared = findNodePath_C(treeTip,treeMatCol0, treeMatCol1,tips[i],rootNode);
    }
    else
    {
      shared = intersect_C(shared,findNodePath_C(treeTip,treeMatCol0, treeMatCol1,tips[i],rootNode));
    }
  }
  //find node fartherest from root
  int longestPath = 0;
  int result = rootNode;
  for(unsigned int i = 0, size = shared.size(); i < size; ++i)
  {
    int rootwayLength = findNodePath_C(treeTip,treeMatCol0, treeMatCol1,shared[i],rootNode).size() - 1;
    if(rootwayLength > longestPath)
    {
      longestPath = rootwayLength;
      result = shared[i];
    }
  }
  return(result);
}

std::vector<std::string> getUpperRank_C(std::vector<std::string> name, std::vector<std::string> treeTip, std::vector<std::string> rankList)
{
  std::vector<std::string> result(name.size());
  for(unsigned int i = 0, size = name.size(); i < size; ++i)
  {
    int index = std::distance(treeTip.begin(),std::find(treeTip.begin(),treeTip.end(),name[i]));
    result[i] = rankList[index];
  }
  return(result);
}

// [[Rcpp::export]]
std::vector<std::string> extractOTUbyRankName_C(std::string rankName, std::vector<std::string> treeTip,std::vector<std::string> rankList)
{
  std::vector<std::string> x;
  x.reserve(rankList.size()/10);
  for(unsigned int i = 0, size=rankList.size(); i < size; ++i)
  {
    if(rankName==rankList[i])
    {
      x.emplace_back(treeTip[i]);
    }
  }
  std::sort(x.begin(),x.end());
  x.erase(std::unique(x.begin(),x.end()),x.end());
  return(x);
}


//get set difference without sorting
std::vector<int> setdiff_C (std::vector<int> vec1,std::vector<int> vec2)
{
  std::vector<int> diff;
  diff.reserve(vec1.size());
  unsigned int size1 = vec1.size();
  unsigned int size2 = vec2.size();
  bool found;
  for(unsigned int i = 0; i < size1; ++i)
  {
    found = false;
    //when an element in vec1 is not found in vec2
    for(unsigned int j = 0; j < size2; ++j)
    {
      if(vec1[i]==vec2[j])
      {
        found = true;
        break;
      }
    }
    if(!found)
    {
      diff.push_back(vec1[i]);
    }
  }
  return(diff);
}

//return subset of vector x of index correcting index by one
inline std::vector<std::string> multiIndicesString_C(std::vector<std::string> x,std::vector<int> index)
{
  unsigned int size = index.size();
  std::vector<std::string> result(size);
  for(unsigned int i = 0; i < size; ++i)
  {
    result[i] = x[index[i]-1];
  }
  return(result);
}

//obtain all nodes that are lower than MRCA of OTUs (including MRCA node)
std::vector<int> getRankSubNodes_C (std::vector<std::string> treeTip,
                                    std::vector<int> treeMatCol0,std::vector<int> treeMatCol1,
                                    std::vector<int> rankTips, std::vector<std::string> rankList)
{
  std::vector<std::string> rankOTU = multiIndicesString_C(treeTip,rankTips);
  int MRCA = findMRCA_C(treeTip,treeMatCol0, treeMatCol1,rankOTU);
  //the index of this OTU name vector equals to the tip number of the OTU
  std::vector<int> lowerTip;
  for(unsigned int i = 0, size=rankOTU.size(); i < size; ++i)
  {
    std::vector<std::string> tempVec(1);
    tempVec[0] = rankOTU[i];
    std::vector<int> tempIndex = getTipNo_C(treeTip,tempVec);
    std::vector<int> temp = findNodePath_C(treeTip,treeMatCol0, treeMatCol1,MRCA,tempIndex[0]);
    lowerTip.insert(lowerTip.end(),temp.begin(),temp.end());
  }
  std::sort(lowerTip.begin(),lowerTip.end());
  lowerTip.erase(std::unique(lowerTip.begin(),lowerTip.end()),lowerTip.end());
  return(lowerTip);
}

std::vector<double> castIntVec2DoubleVec(std::vector<int> vec)
{
  std::vector<double> result(vec.size());
  for(unsigned int i = 0, size=vec.size(); i < size;++i)
  {
    result[i] = (double)vec[i];
  }
  return(result);
}

//get centroid node of rank. When there are multiple nodes, it means the centroid is the branch between the nodes.
// [[Rcpp::export]]
std::vector<double> getRankCentroid_C (std::string rankName, std::vector<int> dropIndex,
                                       std::vector<std::string> treeTip, std::vector<int> treeMatCol0,std::vector<int> treeMatCol1,
                                       std::vector<std::string> rankList,bool show_progress=0, int num_threads=1)
{
  //get tip numbers of the rank
  std::vector<int> rankTips = getTipNo_C(treeTip,extractOTUbyRankName_C(rankName,treeTip,rankList));
  rankTips = setdiff_C(rankTips,dropIndex);
  //when the rank is dropped and no longer remaining, return 0
  if(rankTips.size() == 0)
  {
    std::vector<double> temp(1);
    return(temp);
  }
  //when the rank is monotypic
  if(rankTips.size() == 1)
  {
    return(castIntVec2DoubleVec(rankTips));
  }
  std::vector<int> rankSubNodes = getRankSubNodes_C(treeTip,treeMatCol0, treeMatCol1,rankTips,rankList);
  std::vector<int> subNodesWithoutTips = setdiff_C(rankSubNodes,rankTips);
  subNodesWithoutTips = setdiff_C(subNodesWithoutTips,dropIndex);
  //the total distances from the subnode to all tips of the rank
  std::vector<int> totalDistFromNode(subNodesWithoutTips.size());
  unsigned int size = subNodesWithoutTips.size();
  unsigned int size2 = rankTips.size();
  Progress p(size,show_progress);
  bool stopFlag = false;
#if defined(_OPENMP)
#pragma omp parallel for num_threads(num_threads)
#endif
  for(unsigned int i = 0; i < size; ++i)
  {
    for(unsigned int j = 0; j < size2; ++j)
    {
      //the distance equals the length of nodepath - 1
      if(!Progress::check_abort())
      {
        totalDistFromNode[i] += -1 + (int)findNodePath_C(treeTip,treeMatCol0, treeMatCol1,subNodesWithoutTips[i],rankTips[j]).size();

      }
      else
      {
        stopFlag = true;
      }
    }
    p.increment();
  }
  if(stopFlag)
  {
    stop("Parallel calculation aborted.");
  }
  int minDist = *std::min_element(totalDistFromNode.begin(),totalDistFromNode.end());
  std::vector<int>::iterator it = std::find(totalDistFromNode.begin(),totalDistFromNode.end(),minDist);
  std::vector<int> indices;
  while(it!=totalDistFromNode.end())
  {
    indices.push_back(std::distance(totalDistFromNode.begin(),it));
    it = std::find(++it,totalDistFromNode.end(),minDist);
  }
  std::vector<double> minNode(indices.size());
  unsigned int size3 = indices.size();
#if defined(_OPENMP)
#pragma omp parallel for num_threads(num_threads)
#endif
  for(unsigned int i = 0; i < size3; ++i)
  {
    if(!Progress::check_abort())
    {
      minNode[i] = (double)subNodesWithoutTips[indices[i]];
    }
    else
    {
      stopFlag = true;
    }
  }
  if(stopFlag)
  {
    stop("Parallel calculation aborted.");
  }
  return(minNode);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> getAllCentroids_C (std::vector<std::string> treeTip, std::vector<std::string> allRankNames,std::vector<int> treeMatCol0,std::vector<int> treeMatCol1, std::vector<std::string> rankList, bool show_progress=0, int num_threads=1)
{
  std::vector<std::vector<double>> allCentroid(allRankNames.size());
  unsigned int size = allRankNames.size();
  std::vector<int> dropIndex(0);
  for(unsigned int i = 0; i < size; ++i)
  {
    allCentroid[i] = getRankCentroid_C(allRankNames[i],dropIndex,treeTip,treeMatCol0, treeMatCol1,rankList,true,num_threads);
    if(show_progress)
    {
      double progressPercent = 100*((double)i+1)/(double)size;
      Rprintf("All centroids calculation progress: %f%%\n", progressPercent);
    }
  }
  return(allCentroid);
}

int findRankIndex_C (std::vector<std::string> treeTip, std::string OTU, std::vector<std::string> allRankNames, std::vector<std::string> rankList)
{
  std::vector<std::string> temp(1);
  temp[0] = OTU;
  std::string upperRank = getUpperRank_C(temp,treeTip,rankList)[0];
  return(std::distance(allRankNames.begin(),std::find(allRankNames.begin(),allRankNames.end(),upperRank)));
}

//return shortest nodepath between node and one of the nodes
std::vector<int> minnodepath_C(std::vector<std::string> treeTip, std::vector<int> treeMatCol0,std::vector<int> treeMatCol1,int node, std::vector<int> nodes)
{
  std::vector<int> minpath;
  unsigned int size = nodes.size();
  for(unsigned int i = 0; i < size; ++i)
  {
    std::vector<int> temp = findNodePath_C(treeTip,treeMatCol0,treeMatCol1,node,nodes[i]);
    if(temp.size()<minpath.size()||minpath.size()==0)
    {
      minpath = temp;
    }
  }
  return(minpath);
}

//return centroid which is nearest to the root
int getBaseCentroid_C (std::vector<std::string> treeTip, std::vector<int> treeMatCol0,std::vector<int> treeMatCol1,std::vector<int> centroid)
{
  int rootNode = getRootNode_C(treeTip);
  std::vector<int> temp = minnodepath_C(treeTip,treeMatCol0,treeMatCol1,rootNode,centroid);
  return((int)temp[temp.size() - 1]);
}

//calculate intruder score
// [[Rcpp::export]]
int calcIntScore_C (std::vector<std::string> treeTip, std::vector<int> treeMatCol0,std::vector<int> treeMatCol1,
                    std::string OTU,std::vector<std::vector<int>> allCentroids,std::vector<std::string> allRankNames,std::vector<std::string> rankList)
{
  int score = 0;
  //nodes between the OTU and root
  std::vector<std::string> temp(1);
  temp[0] = OTU;
  std::vector<int> path = findNodePath_C(treeTip,treeMatCol0,treeMatCol1,getRootNode_C(treeTip),getTipNo_C(treeTip,temp)[0]);
  //index of the rank of the OTU
  int skipIndex = findRankIndex_C(treeTip,OTU,allRankNames,rankList);
  unsigned int size = allCentroids.size();
  for(unsigned int i = 0; i < size; ++i)
  {
    if((int)i == skipIndex)
    {
      continue;
    }
    //add one score when an exotic-rank centroid exists on the path from the OTU to the root
    if(intersect_C(path,allCentroids[i]).size()!=0)
    {
      ++score;
    }
  }
  return(score);
}

//return all indices x which fulfill vec[x] == target
std::vector<int> subsetIndex_C(std::vector<int> vec,int target)
{
  std::vector<int> x;
  x.reserve(vec.size());
  unsigned int size = vec.size();
  for(unsigned int i = 0; i < size; ++i)
  {
    if(vec[i]==target)
    {
      x.push_back((int)i);
    }
  }
  return(x);
}
//return subset of vector x of index
inline std::vector<int> multiIndices_C(std::vector<int> x,std::vector<int> index)
{
  unsigned int size = index.size();
  std::vector<int> result(size);
  for(unsigned int i = 0; i < size; ++i)
  {
    result[i] = x[index[i]];
  }
  return(result);
}


std::vector<int> subsetIndexbeta_C(std::vector<int> vec,int target)
{
  std::vector<int> x;
  x.reserve(vec.size());
  unsigned int size = vec.size();
  for(unsigned int i = 0; i < size; ++i)
  {
    if(vec[i]<=target)
    {
      x.push_back((int)i);
    }
  }
  return(x);
}

std::vector<int> subsetIndegamma_C(std::vector<int> vec,int target)
{
  std::vector<int> x;
  unsigned int size = vec.size();
  x.reserve(size);
  for(unsigned int i = 0; i < size; ++i)
  {
    if(vec[i]>target)
    {
      x.push_back((int)i);
    }
  }
  return(x);
}

//obtain tips below a given node
// [[Rcpp::export]]
std::vector<int> findSubTips_C (std::vector<std::string> treeTip,std::vector<int> treeMatCol0, std::vector<int> treeMatCol1,int node)
{
  int tipNo = treeTip.size();
  std::vector<int> belowNode(1);
  belowNode[0] = node;
  std::vector<int> tip;
  if(node <= tipNo)
  {
    return(belowNode);
  }
  while(true)
  {
    std::vector<int> newBelowNode;
    unsigned int size = belowNode.size();
    for(unsigned int i = 0;i < size; ++i)
    {
      newBelowNode = connect_C(newBelowNode,multiIndices_C(treeMatCol1,subsetIndex_C(treeMatCol0,belowNode[i])));
    }
    belowNode = newBelowNode;
    if(belowNode.size()==0)
    {
      break;
    }
    tip = connect_C(tip,multiIndices_C(belowNode,subsetIndexbeta_C(belowNode,tipNo)));
    belowNode = multiIndices_C(belowNode,subsetIndegamma_C(belowNode,tipNo));
    if(belowNode.size()==0)
    {
      break;
    }
  }
  return(tip);
}

inline int countNE (std::vector<std::string> vec, std::string target)
{
  unsigned int size = vec.size();
  int result = 0;
  for(unsigned int i = 0; i < size; ++i)
  {
    if(vec[i] != target)
    {
      ++result;
    }
  }
  return(result);
}

//calculate outlier score
// [[Rcpp::export]]
int calcOutScore_C (std::vector<std::string> treeTip,
                    std::vector<int> treeMatCol0,std::vector<int> treeMatCol1,std::string OTU,
                    std::vector<std::vector<int>> allCentroids,
                    std::vector<std::string> allRankNames,std::vector<std::string> rankList, std::vector<int> dropIndex)
{
  //find the node just below the centroid
  std::vector<std::string> temp(1);
  temp[0] = OTU;
  int otuTip = getTipNo_C(treeTip,temp)[0];
  int rankIndex = findRankIndex_C(treeTip,OTU,allRankNames,rankList);
  //centroid of the rank of OTU
  std::vector<int> rankCentroid = allCentroids[rankIndex];
  //when the centroid is virtually absent due to drop, return 0 (although this situation should not happen)
  if(rankCentroid.size() == 1 && rankCentroid[0] == 0)
  {
    return(0);
  }
  int baseCentroid = getBaseCentroid_C(treeTip,treeMatCol0,treeMatCol1,rankCentroid);
  std::vector<int> centroidToOTU = findNodePath_C(treeTip,treeMatCol0,treeMatCol1,otuTip,baseCentroid);
  centroidToOTU.erase(std::remove(centroidToOTU.begin(),centroidToOTU.end(),baseCentroid));
  //when centroid == otuTip, namely monotypic, return(0)
  if(centroidToOTU.size()==0)
  {
    return(0);
  }
  int centroidBelowNode = centroidToOTU[centroidToOTU.size() - 1];
  std::vector<int> checkingnodes = findNodePath_C(treeTip,treeMatCol0,treeMatCol1,otuTip,centroidBelowNode);
  std::vector<int> maxSubTip;
  unsigned int size = checkingnodes.size();
  for(unsigned int i = 0; i < size; ++i)
  {
    std::vector<int> temp2 = findSubTips_C(treeTip,treeMatCol0,treeMatCol1,checkingnodes[i]);
    temp2 = setdiff_C(temp2,dropIndex);
    if(maxSubTip.size()<temp2.size())
    {
      //maxSubTip.resize(temp2.size());
      //std::copy(temp2.begin(),temp2.end(),maxSubTip.begin());
      maxSubTip = temp2;
    }
  }
  //return(std::accumulate(maxSubTip.begin(),maxSubTip.end(),0));
  std::vector<std::string> belowNames = multiIndicesString_C(treeTip,maxSubTip);
  //check the below tips belong to the same rank as OTU
  std::vector<std::string> belowRank = getUpperRank_C(belowNames,treeTip,rankList);
  std::string OTUupperRank = getUpperRank_C(temp,treeTip,rankList)[0];
  return(countNE(belowRank,OTUupperRank));
}
