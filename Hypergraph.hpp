#pragma once

class Hypergraph {
private:
  int N; // number of vertices
  int M; // number of hyperedges
  std::vector< std::vector<int> > incidenceMatrix; // "sparse" edge-vertex incidence matrix

  std::vector< std::vector<int> > edgeGraph;
  void buildEdgeGraph();
  
public:
  Hypergraph();
  void readIncidenceMatrix();
  void randomHypergraph();
  void sortAndCheck(std::vector< std::vector<int> >&) const;
  
  void ESU(int);
  
  std::vector< std::vector<int> > applyFunction(const std::vector<int>& permutation) const;
  
  bool isEqual(const std::vector< std::vector<int> >&) const;
  
  std::vector< std::vector<int> > getIncidenceMatrix() const;
  
  void enumerateSubgraphs(std::vector<int>, std::vector<int>, int, int);
  
  int getNodeCount() const;
  int getEdgeCount() const;
  void printIncidenceMatrix() const;
  //void printEdgeGraph() const;
};
