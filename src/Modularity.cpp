#pragma once

#include <chrono>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <chrono>

#include <RcppEigen.h>
#include <Rcpp.h>
#include <progress.hpp>

typedef std::vector<int> IVector;
typedef std::vector<double> DVector;

namespace ModularityOptimizer {

class JavaRandom {
private:
  uint64_t seed;
  int next(int bits);
public:
  JavaRandom(uint64_t seed);
  int nextInt(int n);
  void setSeed(uint64_t seed);
};

namespace Arrays2 {
  IVector generateRandomPermutation(int nElements);
  IVector generateRandomPermutation(int nElements, JavaRandom& random);
}


class Clustering {
private:
  int nNodes;
public:
  // Note: These two variables were "protected" in java, which means it is accessible to the whole package/public.
  // Although we could have used friend classes, this allows for better mirroring of the original code.
  int nClusters;
  IVector cluster;
  
  Clustering(int nNodes);
  Clustering(IVector cluster);
  int getNNodes() const {return nNodes;};
  int getNClusters() const {return nClusters;};
  IVector getClusters() const {return cluster;};
  int getCluster(int node) const {return cluster[node];};
  IVector getNNodesPerCluster() const;
  std::vector<IVector> getNodesPerCluster() const;
  void setCluster(int node, int cluster);
  void initSingletonClusters();
  void orderClustersByNNodes();
  void mergeClusters(const Clustering& clustering);

};

class Network {
  friend class VOSClusteringTechnique;
protected:
  int nNodes;
  int nEdges;
  DVector nodeWeight;
  IVector firstNeighborIndex;
  IVector neighbor;
  DVector edgeWeight;
  double totalEdgeWeightSelfLinks;
public:
  Network();
  Network(int nNodes, DVector* nodeWeight, std::vector<IVector>& edge, DVector* edgeWeight);
  Network(int nNodes, std::vector<IVector>& edge) : 
    Network(nNodes, nullptr, edge, nullptr) { };
  Network(int nNodes, DVector* nodeWeight, std::vector<IVector> edge) : 
    Network(nNodes, nodeWeight, edge, nullptr) {};
  Network(int nNodes, std::vector<IVector>& edge, DVector* edgeWeight) :
    Network(nNodes, nullptr, edge, edgeWeight) {};
  
  Network(int nNodes, DVector* nodeWeight, IVector& firstNeighborIndex, IVector& neighbor, DVector* edgeWeight);
  Network(int nNodes, IVector& firstNeighborIndex, IVector& neighbor) : 
    Network(nNodes, nullptr, firstNeighborIndex, neighbor, nullptr) {};
  
  Network(int nNodes, DVector* nodeWeight, IVector& firstNeighborIndex, IVector& neighbor) :
    Network(nNodes, nodeWeight, firstNeighborIndex, neighbor, nullptr){};
  
  Network(int nNodes, IVector& firstNeighborIndex, IVector& neighbor, DVector* edgeWeight) :
    Network(nNodes, nullptr, firstNeighborIndex, neighbor, edgeWeight) {};
  

 int getNNodes() {return nNodes;};
 double getTotalNodeWeight();
 DVector getNodeWeights();
 double getNodeWeight(int node) { return nodeWeight.at(node);};
 int getNEdges() {return nEdges / 2;};
 int getNEdges(int node) {return firstNeighborIndex.at(node + 1) - firstNeighborIndex.at(node);};
 IVector getNEdgesPerNode();
 std::vector<IVector> getEdges();
 IVector getEdges(int node);
 std::vector<IVector> getEdgesPerNode();
 double getTotalEdgeWeight();
 double getTotalEdgeWeight(int node);
 DVector getTotalEdgeWeightPerNode();
 DVector getEdgeWeights() {return edgeWeight;};
 DVector getEdgeWeights(int node);
 std::vector<DVector> getEdgeWeightsPerNode();
 double getTotalEdgeWeightSelfLinks()
 {
   return totalEdgeWeightSelfLinks;
 };
 // Added these to avoid making these values public
 int getFirstNeighborIndexValue(int i) const {
   return firstNeighborIndex.at(i);
 };
 int getNeighborValue(int index) const {
   return neighbor.at(index);
 }

 Network createReducedNetwork(const Clustering& clustering) const;
private:
  double generateRandomNumber(int node1, int node2, const IVector& nodePermutation);
};


class VOSClusteringTechnique {
private:
  std::shared_ptr<Network> network;
  std::shared_ptr<Clustering> clustering;
  double resolution;

public:
  VOSClusteringTechnique(std::shared_ptr<Network> network, double resolution);
  VOSClusteringTechnique(std::shared_ptr<Network> network, std::shared_ptr<Clustering> clustering, double resolution);
  std::shared_ptr<Network> getNetwork() { return network;}
  std::shared_ptr<Clustering> getClustering()  { return clustering; }
  double getResolution() {return resolution; }
  void setNetwork(std::shared_ptr<Network> network) {this->network = network;}
  void setClustering(std::shared_ptr<Clustering> clustering) {this->clustering = clustering;}
  void setResolution(double resolution) {this->resolution = resolution;}
  double calcQualityFunction();

  bool runLocalMovingAlgorithm(JavaRandom& random);
  bool runLouvainAlgorithm(JavaRandom& random);
  //bool runIteratedLouvainAlgorithm(int maxNIterations, JavaRandom& random);

  int removeCluster(int cluster);
};

std::shared_ptr<Network> matrixToNetwork(IVector& node1, IVector& node2, DVector& edgeWeight1, int modularityFunction, int nNodes);
};


using namespace ModularityOptimizer;
using namespace std::chrono;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
//' lovain modularity
//' @export
//' 
// [[Rcpp::export]]
IntegerVector RunModularityClusteringCpp(Eigen::SparseMatrix<double> SNN,
    int modularityFunction,
    double resolution,
    int algorithm,
    int nRandomStarts,
    int nIterations,
    int randomSeed,
    bool printOutput) {

  
try {
  bool update;
  double modularity, maxModularity, resolution2;
  int i, j;

  // Load netwrok
  std::shared_ptr<Network> network;
  // Load lower triangle
  int network_size = (SNN.nonZeros() / 2) + 3;
  IVector node1;
  IVector node2;
  DVector edgeweights;
  node1.reserve(network_size);
  node2.reserve(network_size);
  edgeweights.reserve(network_size);
  for (int k=0; k < SNN.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, k); it; ++it){
      if(it.col() >= it.row()){
        continue;
      }
      node1.emplace_back(it.col());
      node2.emplace_back(it.row());
      edgeweights.emplace_back(it.value());
    }
  }
  if(node1.size() == 0) {
    stop("Matrix contained no network data.  Check format.");
  }
  int nNodes = std::max(SNN.cols(), SNN.rows());
  network = matrixToNetwork(node1, node2, edgeweights, modularityFunction, nNodes);
  Rcpp::checkUserInterrupt();

  resolution2 = ((modularityFunction == 1) ? (resolution / (2 * network->getTotalEdgeWeight() + network->getTotalEdgeWeightSelfLinks())) : resolution);

  std::shared_ptr<Clustering> clustering;
  maxModularity = -std::numeric_limits<double>::infinity();
  JavaRandom random(randomSeed);

  Progress p(nRandomStarts, printOutput);
  for (i = 0; i < nRandomStarts; i++)
  {

    VOSClusteringTechnique vosClusteringTechnique(network, resolution2);

    j = 0;
    update = true;
    do
    {

      if (algorithm == 1)
        update = vosClusteringTechnique.runLouvainAlgorithm(random);
      j++;

      modularity = vosClusteringTechnique.calcQualityFunction();

      Rcpp::checkUserInterrupt();
    }
    while ((j < nIterations) && update);

    if (modularity > maxModularity)
    {
      clustering = vosClusteringTechnique.getClustering();
      maxModularity = modularity;
    }

    p.increment();
  }
  if(clustering == nullptr) {
    stop("Clustering step failed.");
  }

  // Return results
  clustering->orderClustersByNNodes();
  IntegerVector iv(clustering->cluster.cbegin(), clustering->cluster.cend());
  return iv;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return IntegerVector(1);
}

JavaRandom::JavaRandom(uint64_t seed) {
  setSeed(seed);
}

void JavaRandom::setSeed(uint64_t seed) {

  this->seed = (seed ^ uint64_t(0x5DEECE66D)) & ((uint64_t(1) << 48) - 1);
}

int JavaRandom::next(int bits) {
  // Only 31 bits ever used.
  seed = (seed * uint64_t(0x5DEECE66D) + uint64_t(0xB)) & ((uint64_t(1) << 48) - 1);
  return static_cast<int>(seed >> (48 - bits));
}


int JavaRandom::nextInt(int n) {
  if (n <= 0)
    throw std::out_of_range("n must be positive");
  if ((n & -n) == n) // i.e., n is a power of 2
    return static_cast<int>((static_cast<uint64_t>(n) * static_cast<uint64_t>(next(31))) >> 31);
  int bits, val;
  do
  {
    bits = next(31);
    val = bits % n;
  }
  while (bits - val + (n - 1) < 0);
  return val;
}


IVector Arrays2::generateRandomPermutation(int nElements, JavaRandom& random)
{
    IVector permutation(nElements, 0);
    for (int i = 0; i < nElements; i++)
      permutation[i] = i;
    for (int i = 0; i < nElements; i++)
    {
      int j = random.nextInt(nElements);
      int k = permutation[i];
      permutation[i] = permutation[j];
      permutation[j] = k;
    }
    return permutation;
}


Clustering::Clustering(int nNodes): 
  nNodes(nNodes),
  nClusters(1),
  cluster(nNodes)
{};

Clustering::Clustering(IVector cluster) :
  nNodes(cluster.size()),
  cluster(cluster.cbegin(), cluster.cend()) 
  {
    nClusters = *std::max_element(cluster.cbegin(), cluster.cend()) + 1;  
  }


IVector Clustering::getNNodesPerCluster() const {
  IVector nNodesPerCluster(nClusters, 0);
  for(const int& clust: cluster) {
    nNodesPerCluster.at(clust)++;
  }
  return nNodesPerCluster;
}

std::vector<IVector> Clustering::getNodesPerCluster() const {
  std::vector<IVector> nodePerCluster(nClusters);
  IVector nNodesPerCluster = getNNodesPerCluster();
  for(int i =0; i < nClusters; i++)
  {
    const int cnt = nNodesPerCluster.at(i);
    nodePerCluster.at(i).reserve(cnt);
  }
  for(int i=0; i< nNodes; i++) {
    nodePerCluster.at(cluster.at(i)).push_back(i);
  }
  return nodePerCluster;
}

void Clustering::setCluster(int node, int cluster) {
  this->cluster.at(node) = cluster;
  nClusters = std::max(nClusters, cluster+1);
}

void Clustering::initSingletonClusters() {
  for(int i=0; i < nNodes; i++) {
    cluster.at(i) = i;
  }
  nClusters = nNodes;
}

void Clustering::orderClustersByNNodes() {
  typedef std::pair<int, int> ipair; // holds numNodes, cluster
  std::vector<ipair> clusterNNodes;
  clusterNNodes.reserve(nClusters);
  IVector nNodesPerCluster = getNNodesPerCluster();
  for(int i=0; i<nClusters; i++) {
    clusterNNodes.push_back(std::make_pair(nNodesPerCluster.at(i), i));
  }
  // Note order is descending
  stable_sort(clusterNNodes.begin(), clusterNNodes.end(), 
       [](const std::pair<int, int>&a, const std::pair<int, int>& b) {
         return b.first < a.first;
       });
       //std::greater<ipair>());
  
  // now make a map from old to new names
  IVector newCluster(nClusters, 0);
  int i=0;
  do {
    newCluster[clusterNNodes[i].second] = i;
    i++;
  } while (i < nClusters && clusterNNodes[i].first > 0);
  nClusters = i;
  for(int i=0; i<nNodes; i++) {
    cluster[i] = newCluster[cluster[i]];
  }
}

void Clustering::mergeClusters(const Clustering& clustering) {
  for (int i = 0; i < nNodes; i++)
    cluster.at(i) = clustering.cluster.at(cluster.at(i));
  nClusters = clustering.nClusters;
}


Network::Network() {};

Network::Network(int nNodes, DVector* nodeWeight, IVector& firstNeighborIndex, IVector& neighbor, DVector* edgeWeight) :
  nNodes(nNodes),
  nEdges(neighbor.size()),
  nodeWeight(nNodes),
  firstNeighborIndex(firstNeighborIndex.cbegin(), firstNeighborIndex.cend()),
  neighbor(neighbor.cbegin(), neighbor.cend()),
  edgeWeight(nEdges, 1.0),
  totalEdgeWeightSelfLinks(0)
  {
  
  if (edgeWeight != nullptr)
    std::copy(edgeWeight->cbegin(), edgeWeight->cend(), this->edgeWeight.begin());

  if (nodeWeight != nullptr) {
    std::copy(nodeWeight->cbegin(), nodeWeight->cend(), this->nodeWeight.begin());
  } else {
    this->nodeWeight = getTotalEdgeWeightPerNode();
  }
}


Network::Network(int nNodes, DVector* nodeWeight, std::vector<IVector>& edge, DVector* edgeWeight) :
  nNodes(nNodes),
  nEdges(0),
  nodeWeight(),
  firstNeighborIndex(nNodes + 1, 0),
  neighbor(),
  edgeWeight(),
  totalEdgeWeightSelfLinks(0)
{
  if(edge.size() != 2 || edge[0].size() != edge[1].size()) {
    throw std::length_error("Edge was supposed to be an array with 2 columns of equal size.");
  }
  IVector neighbor(edge.at(0).size(), 0);
  DVector edgeWeight2(edge.at(0).size(), 0.0);
  
  int i = 1;
  for (size_t j = 0; j < edge[0].size(); j++)
    if (edge[0][j] != edge[1][j])
    {
      if (edge[0][j] >= i)
        for (; i <= edge[0][j]; i++)
          firstNeighborIndex.at(i) = nEdges;
      neighbor[nEdges] = edge[1][j];
      edgeWeight2[nEdges] = (edgeWeight != nullptr) ? (*edgeWeight)[j] : 1.0;
      nEdges++;
    }
    else
      totalEdgeWeightSelfLinks += (edgeWeight != nullptr) ? (*edgeWeight)[j] : 1;
    for (; i <= nNodes; i++)
      firstNeighborIndex.at(i) = nEdges;
    
    this->neighbor.resize(nEdges);
    std::copy(neighbor.begin(), neighbor.begin() + nEdges, this->neighbor.begin());
    this->edgeWeight.resize(nEdges);
    std::copy(edgeWeight2.begin(), edgeWeight2.begin() + nEdges, this->edgeWeight.begin());
    if(nodeWeight == nullptr) {
      this->nodeWeight = getTotalEdgeWeightPerNode();
    } else {
      this->nodeWeight = *nodeWeight;
    }
}

double Network::getTotalNodeWeight() {
  return std::accumulate(nodeWeight.cbegin(), nodeWeight.cend(), 0.0);
}
DVector Network::getNodeWeights() {
  return nodeWeight;
}
IVector Network::getNEdgesPerNode() {
  IVector nEdgesPerNode(nNodes, 0);
  for(int i=0; i< nNodes; i++) {
    nEdgesPerNode.at(i) = firstNeighborIndex.at(i + 1) - firstNeighborIndex.at(i);
  }
  return nEdgesPerNode;
}
std::vector<IVector> Network::getEdges() {
  std::vector<IVector> edge(2);
  edge[0].resize(nEdges);
  for(int i=0; i < nNodes; i++) {
    std::fill(edge[0].begin() + firstNeighborIndex.at(i), edge[0].begin() + firstNeighborIndex.at(i + 1), i);
  }
  edge.at(1) = neighbor;
  return edge;
}
IVector Network::getEdges(int node) {
  return IVector(neighbor.begin() + firstNeighborIndex.at(node),
                 neighbor.begin() + firstNeighborIndex.at(node + 1));
}

std::vector<IVector> Network::getEdgesPerNode() {
  std::vector<IVector> edgePerNode(nNodes);
  for (int i = 0; i < nNodes; i++) {
    edgePerNode[i] = IVector(neighbor.begin() + firstNeighborIndex.at(i),
                             neighbor.begin() + firstNeighborIndex.at(i + 1));
  }
  return edgePerNode;
}

double Network::getTotalEdgeWeight() {
  return std::accumulate(edgeWeight.cbegin(), edgeWeight.cend(), 0.0) / 2.0;
}

double Network::getTotalEdgeWeight(int node) {
  return std::accumulate(edgeWeight.cbegin() + firstNeighborIndex.at(node),
                         edgeWeight.cbegin() + firstNeighborIndex.at(node + 1),
                         0.0);
}

DVector Network::getTotalEdgeWeightPerNode() {
  DVector totalEdgeWeightPerNode(nNodes, 0.0);
  for (int i = 0; i < nNodes; i++) {
    totalEdgeWeightPerNode[i] = getTotalEdgeWeight(i);
  }
  return totalEdgeWeightPerNode;
}

DVector Network::getEdgeWeights(int node) {
  return DVector(edgeWeight.cbegin() + firstNeighborIndex.at(node),
                 edgeWeight.cbegin() + firstNeighborIndex.at(node+1));
}

std::vector<DVector> Network::getEdgeWeightsPerNode() {
  std::vector<DVector> edgeWeightPerNode(nNodes);
  for (int i = 0; i < nNodes; i++)
    edgeWeightPerNode[i] = getEdgeWeights(i);
  return edgeWeightPerNode;
}


// Skipping unused Network creators
// Network createNetworkWithoutNodeWeights()
// Network createNetworkWithoutEdgeWeights()
// Network createNetworkWithoutNodeAndEdgeWeights()
// Network createNormalizedNetwork1()
// Network createNormalizedNetwork2()
// Network createPrunedNetwork(int nEdges)
// Network createPrunedNetwork(int nEdges, Random random)
// Network createSubnetwork(int[] node)
// Network createSubnetwork(boolean[] nodeInSubnetwork)
// Network createSubnetwork(Clustering clustering, int cluster)
// Network createSubnetworkLargestComponent()
// Network createReducedNetwork(Clustering clustering)
Network Network::createReducedNetwork(const Clustering& clustering) const {
  Network reducedNetwork;
  reducedNetwork.nNodes = clustering.nClusters;
  
  reducedNetwork.nEdges = 0;
  reducedNetwork.nodeWeight = DVector(clustering.nClusters);
  reducedNetwork.firstNeighborIndex = IVector(clustering.nClusters + 1);
  reducedNetwork.totalEdgeWeightSelfLinks = totalEdgeWeightSelfLinks;
  IVector reducedNetworkNeighbor1(nEdges);
  DVector reducedNetworkEdgeWeight1(nEdges);
  IVector reducedNetworkNeighbor2(clustering.nClusters - 1);
  DVector reducedNetworkEdgeWeight2(clustering.nClusters);
  std::vector<IVector> nodePerCluster = clustering.getNodesPerCluster();
  for (int i = 0; i < clustering.nClusters; i++)
  {
    int j = 0;
    for (size_t k = 0; k < nodePerCluster[i].size(); k++)
    {
      int l = nodePerCluster[i][k];
      
      reducedNetwork.nodeWeight[i] += nodeWeight[l];
      
      for (int m = firstNeighborIndex[l]; m < firstNeighborIndex[l + 1]; m++)
      {
        int n = clustering.cluster[neighbor[m]];
        if (n != i)
        {
          if (reducedNetworkEdgeWeight2[n] == 0)
          {
            reducedNetworkNeighbor2[j] = n;
            j++;
          }
          reducedNetworkEdgeWeight2[n] += edgeWeight[m];
        }
        else
          reducedNetwork.totalEdgeWeightSelfLinks += edgeWeight[m];
      }
    }
    
    for (int k = 0; k < j; k++)
    {
      reducedNetworkNeighbor1[reducedNetwork.nEdges + k] = reducedNetworkNeighbor2[k];
      reducedNetworkEdgeWeight1[reducedNetwork.nEdges + k] = reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[k]];
      reducedNetworkEdgeWeight2[reducedNetworkNeighbor2[k]] = 0;
    }
    reducedNetwork.nEdges += j;
    reducedNetwork.firstNeighborIndex[i + 1] = reducedNetwork.nEdges;
  }
  reducedNetwork.neighbor = IVector(reducedNetworkNeighbor1.cbegin(), reducedNetworkNeighbor1.cbegin() + reducedNetwork.nEdges);
  reducedNetwork.edgeWeight = DVector(reducedNetworkEdgeWeight1.cbegin(), reducedNetworkEdgeWeight1.cbegin() + reducedNetwork.nEdges);
  return reducedNetwork;
}

VOSClusteringTechnique::VOSClusteringTechnique(std::shared_ptr<Network> network, double resolution) :
  network(network),
  resolution(resolution)
  {
    clustering = std::make_shared<Clustering>(network->getNNodes());
    clustering->initSingletonClusters();
  };

VOSClusteringTechnique::VOSClusteringTechnique(std::shared_ptr<Network> network, std::shared_ptr<Clustering> clustering, double resolution) :
  network(network),
  clustering(clustering),
  resolution(resolution){};

double VOSClusteringTechnique::calcQualityFunction() {
  double qualityFunction = 0.0;
  for (int i = 0; i < network->getNNodes(); i++)
  {
    int j = clustering->cluster[i];
    for (int k = network->getFirstNeighborIndexValue(i); k < network->getFirstNeighborIndexValue(i + 1); k++)
      if (clustering->cluster[network->getNeighborValue(k)] == j)
        qualityFunction += network->edgeWeight[k];
  }
  qualityFunction += network->totalEdgeWeightSelfLinks;
  
  DVector clusterWeight(clustering->nClusters);
  for (int i = 0; i < network->nNodes; i++)
    clusterWeight[clustering->cluster[i]] += network->nodeWeight[i];
  for (int i = 0; i < clustering->nClusters; i++)
    qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution;
  
  qualityFunction /= 2 * network->getTotalEdgeWeight() + network->totalEdgeWeightSelfLinks;
  
  return qualityFunction;
}

bool VOSClusteringTechnique::runLocalMovingAlgorithm(JavaRandom& random){
  bool update = false;
  double maxQualityFunction, qualityFunction;
  DVector clusterWeight(network->getNNodes(), 0); 
  IVector nNodesPerCluster(network->getNNodes(), 0);
  
  int bestCluster, j, k, l, nNeighboringClusters, nStableNodes;
  if (network->getNNodes() == 1)
    return false;
  
  for (int i = 0; i < network->getNNodes(); i++)
  {
    clusterWeight[clustering->cluster[i]] += network->nodeWeight[i];
    nNodesPerCluster[clustering->cluster[i]]++;
  }
  
  int nUnusedClusters = 0;
  IVector unusedCluster(network->getNNodes(), 0);
  for (int i = 0; i < network->getNNodes(); i++) {
    if (nNodesPerCluster[i] == 0)
    {
      unusedCluster[nUnusedClusters] = i;
      nUnusedClusters++;
    }
  }
    
  IVector nodePermutation = Arrays2::generateRandomPermutation(network->nNodes, random);
  DVector edgeWeightPerCluster(network->getNNodes(), 0.0);
  IVector neighboringCluster(network->getNNodes() - 1, 0);
  nStableNodes = 0;
  int i = 0;
  do {
    j = nodePermutation[i];
    nNeighboringClusters = 0;
    for (k = network->firstNeighborIndex.at(j); k < network->firstNeighborIndex.at(j + 1); k++)
      {
        l = clustering->cluster[network->neighbor[k]];
        if (edgeWeightPerCluster[l] == 0)
        {
          neighboringCluster[nNeighboringClusters] = l;
          nNeighboringClusters++;
        }
        edgeWeightPerCluster[l] += network->edgeWeight[k];
      }
      
      clusterWeight[clustering->cluster[j]] -= network->nodeWeight[j];
      nNodesPerCluster[clustering->cluster[j]]--;
      if (nNodesPerCluster[clustering->cluster[j]] == 0)
      {
        unusedCluster[nUnusedClusters] = clustering->cluster[j];
        nUnusedClusters++;
      }
      
      bestCluster = -1;
      maxQualityFunction = 0;
      for (k = 0; k < nNeighboringClusters; k++)
      {
        l = neighboringCluster[k];
        qualityFunction = edgeWeightPerCluster[l] - network->nodeWeight[j] * clusterWeight[l] * resolution;
        if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) && (l < bestCluster)))
        {
          bestCluster = l;
          maxQualityFunction = qualityFunction;
        }
        edgeWeightPerCluster[l] = 0;
      }
      if (maxQualityFunction == 0)
      {
        bestCluster = unusedCluster[nUnusedClusters - 1];
        nUnusedClusters--;
      }
      
      clusterWeight[bestCluster] += network->nodeWeight[j];
      nNodesPerCluster[bestCluster]++;
      if (bestCluster == clustering->cluster[j])
        nStableNodes++;
      else
      {
        clustering->cluster[j] = bestCluster;
        nStableNodes = 1;
        update = true;
      }
      
      i = (i < network->nNodes - 1) ? (i + 1) : 0;
    }
    while (nStableNodes < network->nNodes);
    
    IVector newCluster(network->getNNodes());
    clustering->nClusters = 0;
    for (i = 0; i < network->nNodes; i++)
      if (nNodesPerCluster[i] > 0)
      {
        newCluster[i] = clustering->nClusters;
        clustering->nClusters++;
      }
    for (i = 0; i < network->nNodes; i++)
        clustering->cluster[i] = newCluster[clustering->cluster[i]];
      
      return update;
}

bool VOSClusteringTechnique::runLouvainAlgorithm(JavaRandom& random) {

  if (network->nNodes == 1)
    return false;
  bool update = runLocalMovingAlgorithm(random);
  if (clustering->nClusters < network->nNodes)
  {
    VOSClusteringTechnique vosClusteringTechnique(std::make_shared<Network>(network->createReducedNetwork(*clustering)), resolution);
    
    bool update2 = vosClusteringTechnique.runLouvainAlgorithm(random);
    
    if (update2)
    {
      update = true;
      
      clustering->mergeClusters(*vosClusteringTechnique.clustering);
    }
  }
  return update;
}

std::shared_ptr<Network> ModularityOptimizer::matrixToNetwork(IVector& node1, IVector& node2, DVector& edgeWeight1, int modularityFunction, int nNodes) {
  
  int n1_max = *std::max_element(node1.cbegin(), node1.cend());
  int n2_max = *std::max_element(node2.cbegin(), node2.cend());
  IVector nNeighbors(nNodes);
  for (size_t i = 0; i < node1.size(); i++)
    if (node1[i] < node2[i])
    {
      nNeighbors[node1[i]]++;
      nNeighbors[node2[i]]++;
    }
    
    IVector firstNeighborIndex(nNodes + 1);
    int nEdges = 0;
    for (int i = 0; i < nNodes; i++)
    {
      firstNeighborIndex[i] = nEdges;
      nEdges += nNeighbors[i];
    }
    firstNeighborIndex[nNodes] = nEdges;
    
    IVector neighbor(nEdges);
    DVector edgeWeight2(nEdges);
    std::fill(nNeighbors.begin(), nNeighbors.end(), 0);
    for (size_t i = 0; i < node1.size(); i++)
      if (node1[i] < node2[i])
      {
        int j = firstNeighborIndex[node1[i]] + nNeighbors[node1[i]];
        neighbor[j] = node2[i];
        edgeWeight2[j] = edgeWeight1[i];
        nNeighbors[node1[i]]++;
        j = firstNeighborIndex[node2[i]] + nNeighbors[node2[i]];
        neighbor[j] = node1[i];
        edgeWeight2[j] = edgeWeight1[i];
        nNeighbors[node2[i]]++;
      }
      
      if (modularityFunction == 1)
        return std::make_shared<Network>(nNodes, firstNeighborIndex, neighbor, &edgeWeight2);
}
#ifdef STANDALONE

template<typename T>
void input(std::string msg, T& value) {
  std::cout << msg << std::endl << std::endl;
  std::cin >> value;
}

#endif