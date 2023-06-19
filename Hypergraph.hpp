#pragma once

#include "Settings.hpp"


// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933#72073933
struct HashFunction {
  std::size_t operator()(std::vector<int> const& vec) const {
    std::size_t seed = vec.size();
    for(auto x : vec) {
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = (x >> 16) ^ x;
      seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

class Hypergraph {
private:
  int N; // number of vertices
  int M; // number of hyperedges
  int K; // max. degree of a given edge
  std::unordered_set< std::vector<int>, HashFunction > hashEdge;
  std::vector< std::vector<int> > incidenceMatrix; // "sparse" edge-vertex incidence matrix
  std::vector< std::vector<int> > edgeBySize; 
  void sortAndCheck(std::vector< std::vector<int> >&);
  void compress();
  void readIncidenceMatrix(std::istream&);
  void printIncidenceMatrix(std::ostream&);
  bool shuffleEdges(int, int);
public:
  Hypergraph();
  void randomHypergraph(int, int, int);
  std::vector< std::vector<int> > applyFunction(const std::vector<int>& permutation);
  bool isEqual(const std::vector< std::vector<int> >&);
  std::vector< std::vector<int> > buildEdgeGraph();
  std::vector< std::vector<int> > buildVertexGraph(int);
  Hypergraph filterEdge(int);
  std::vector< std::vector<int> > getIncidenceMatrix();
  int getNodeCount();
  int getEdgeCount();
  int getEdgeMaxDeg();
  bool is_two_connected();
  std::vector< std::vector<int> > getGraph();
  bool validEdge(std::vector<int>); 
  void printEdgeSubgraph(std::vector< std::pair<int, int> >&); 
  void setN(int);
  void setM(int);
  void setIncidenceMatrix(std::vector< std::vector<int> >& );
  std::vector<int> getEdge(int);
  void printIncidenceMatrix(); // print to stdout
  void saveToFile(std::string); // print to file in specified format
  void readFromStdin(); // read from stdin
  void readFromFile(std::string); // read from file
  std::vector<std::vector<int>> getDegreeSequence();
  void shuffleHypergraph (int);
  Hypergraph induceSubgraph(const std::vector<int>&);
  Hypergraph induceSubgraphNoComp(const std::vector<int>&, const std::vector<int>&);
};
