#pragma once

#include "Settings.hpp"

class ESU {
private:
  static int pos[MAX_INPUT_N];
  static std::vector< std::pair<int, int> > edgeList;
  static std::set<int> subgraph;
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
  static int getNodeCount();
  static void setGraph(const std::vector< std::vector<int> >&);
  static void clearDataStruct();
  static void enumerateSubgraphs();
  static void setupAndRun(const std::vector< std::vector<int> >& , int);
  static Hypergraph h;
  static std::map< std::vector<graph>, int> counterHyper;
  static std::set< std::vector<int> > visited;
public:
  static std::vector< std::vector< std::pair<int, int> > > startEdgeGraphSubgraphs(Hypergraph&, int);
  static std::vector< std::pair<int, std::string> > getEquivalenceClass(const std::vector< std::vector<int> >&, int);
  static std::vector< std::vector< std::pair<int, int> > > getAllSubgraphs(const std::vector< std::vector<int> >&, int);
  static void k3(Hypergraph&);
  static void k4(Hypergraph&);
};
