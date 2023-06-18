#pragma once

#include "Settings.hpp"

class ESU {
private:
  static int pos[MAX_INPUT_N];
  static std::vector< std::pair<int, int> > edgeList;
  static std::vector<int> subgraph;
  static std::vector<int> subgraph_compressed;
  static std::vector< std::vector< std::pair<int, int> > > subgraphs;
  static std::set<int> extension;
  static std::set<int> neiSubgraph;
  static std::stack< std::set<int> > saveExtension;
  static std::vector< std::vector<int> > g;
  static std::vector< std::vector<int> > removeExt;
  static std::vector< std::vector<int> > removeNei;
  static int K; // Subgraph size
  static int V; // Initial node
  static std::map<std::string, int> counter;
  //static int getNodeCount();
  //static void setGraph(const std::vector< std::vector<int> >&, int);
  
  static void enumerateSubgraphs(std::vector<int>);
  static void setupAndRun(const std::vector< std::vector<int> >& , int);
  static Hypergraph h;
  static std::map< std::vector<graph>, int> counterHyper;
  //static std::map< std::vector<graph>, int> counterHyperK4;
  //static std::map< std::vector<graph>, int> counterHyperBF4;
  
  
  static std::map< std::vector<graph>, int> counterHyperK4;
  static std::map< std::vector<graph>, int> counterHyperK3;
  static std::map< std::vector<graph>, int> counterHyperK3D;
  static std::map< std::vector<graph>, int> counterHyperBF4;
  static std::map< std::vector<graph>, int> counterHyperBF3;
  

  static int f[MAX_INPUT_N];
  //static int ESU_FIXED_extension[MAX_INPUT_N];
  //static int ptr;


  
  
  static std::set< std::vector<int> > visited;
public:
  static void clearDataStruct();
  static std::vector< std::vector< std::pair<int, int> > > startEdgeGraphSubgraphs(Hypergraph&, int);
  static std::vector< std::pair<int, std::string> > getEquivalenceClass(const std::vector< std::vector<int> >&, int);
  static std::vector< std::vector< std::pair<int, int> > > getAllSubgraphs(const std::vector< std::vector<int> >&, int);
  static std::map< std::vector<graph>, int> k3(Hypergraph&);
  static std::map< std::vector<graph>, int> k4(Hypergraph&);
  static std::map< std::vector<graph>, int> bruteForce3(Hypergraph& inputGraph);
  static std::map< std::vector<graph>, int> bruteForce4(Hypergraph& inputGraph);
  
  static std::map< std::vector<graph>, int> k3Modified(Hypergraph& inputGraph);
};
