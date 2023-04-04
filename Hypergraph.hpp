#pragma once

class Hypergraph {
private:
  int N; // number of vertices
  int M; // number of hyperedges
  std::set< std::vector<int> > hashEdge;
  std::vector< std::vector<int> > incidenceMatrix; // "sparse" edge-vertex incidence matrix
public:
  Hypergraph();
  Hypergraph filterEdge(int);
  void readIncidenceMatrix();
  void randomHypergraph();
  void sortAndCheck(std::vector< std::vector<int> >&);
  std::vector< std::vector<int> > applyFunction(const std::vector<int>& permutation);
  bool isEqual(const std::vector< std::vector<int> >&) ;
  std::vector< std::vector<int> > getIncidenceMatrix() ;
  int getNodeCount() ;
  int getEdgeCount() ;
  void printIncidenceMatrix() ;
  void printEdgeSubgraph(std::vector< std::pair<int, int> >&); 
  std::vector< std::vector<int> > buildEdgeGraph();
  void setN(int);
  void setM(int);
  void setIncidenceMatrix(std::vector< std::vector<int> >& );
  Hypergraph induceSubgraph(const std::vector<int>&);
};
