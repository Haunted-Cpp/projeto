#pragma once

class Hypergraph {
private:
  int N; // number of vertices
  int M; // number of hyperedges
  std::vector< std::vector<int> > incidenceMatrix; // "sparse" edge-vertex incidence matrix
public:
  Hypergraph();
  void readIncidenceMatrix();
  void randomHypergraph();
  void sortAndCheck(std::vector< std::vector<int> >&) const;
  std::vector< std::vector<int> > applyFunction(const std::vector<int>& permutation) const;
  bool isEqual(const std::vector< std::vector<int> >&) const;
  std::vector< std::vector<int> > getIncidenceMatrix() const;
  int getNodeCount() const;
  int getEdgeCount() const;
  void printIncidenceMatrix() const;
  std::vector< std::vector<int> > buildEdgeGraph();
};
