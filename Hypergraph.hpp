#pragma once

class Hypergraph {
private:
  int N; // number of vertices
  int M; // number of hyperedges
  std::set< std::vector<int> > hashEdge;
  std::vector< std::vector<int> > incidenceMatrix; // "sparse" edge-vertex incidence matrix
  void sortAndCheck(std::vector< std::vector<int> >&);
  void compress();
  void readIncidenceMatrix(std::istream&);
  void printIncidenceMatrix(std::ostream&);
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
  
  
  Hypergraph induceSubgraph(const std::vector<int>&);
  Hypergraph induceSubgraphNoComp(const std::vector<int>&);
};
