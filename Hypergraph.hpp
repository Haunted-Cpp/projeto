#pragma once

class Hypergraph {
private:
  int N; // number of vertices
  int M; // number of hyperedges
  std::set< std::vector<int> > hashEdge;
  std::vector< std::vector<int> > incidenceMatrix; // "sparse" edge-vertex incidence matrix
  
  void sortAndCheck(std::vector< std::vector<int> >&);
  void randomHypergraph();
  void compress();
  
  
public:
  Hypergraph();
  
  
  std::vector< std::vector<int> > applyFunction(const std::vector<int>& permutation);
  bool isEqual(const std::vector< std::vector<int> >&);
  
  void readIncidenceMatrix();
  std::vector< std::vector<int> > buildEdgeGraph();
  std::vector< std::vector<int> > buildVertexGraph(int);
  Hypergraph filterEdge(int);
  
  
  std::vector< std::vector<int> > getIncidenceMatrix();
  int getNodeCount();
  int getEdgeCount();
  int getEdgeMaxDeg();
  bool validEdge(std::vector<int>); 
  void printIncidenceMatrix();
  void printEdgeSubgraph(std::vector< std::pair<int, int> >&); 
  
  
  
  void setN(int);
  void setM(int);
  void setIncidenceMatrix(std::vector< std::vector<int> >& );
  
  std::vector<int> getEdge(int);
  
  Hypergraph induceSubgraph(const std::vector<int>&);
};
