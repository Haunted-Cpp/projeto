#pragma once

#include "Settings.hpp"

class ESU {
private:
  static int pos[MAX_INPUT_N];
  static int f[MAX_INPUT_N];
  static std::vector< std::pair<int, int> > edgeList;
  static std::vector<int> subgraph;
  static std::vector<int> subgraph_compressed;
  static std::vector< std::vector< std::pair<int, int> > > subgraphs;
  static std::vector< std::vector<int> > g;
  static int K; // Subgraph size
  static int V; // Initial node
  static std::map<std::string, long long> counter;
  static void enumerateSubgraphs(std::vector<int>);
  static void setupAndRun(const std::vector< std::vector<int> >& , int);
  static Hypergraph h;
  static std::map< int, long long> counterHyper;
  static std::unordered_set< std::vector<int>, HashFunction > visited;
  static std::map<std::string, long long> FaSE(const std::vector<std::pair<int, int> > edges, int k);
  
  
  static Hypergraph binaryToHyper(std::string str, int k);
  
  
  static void clearDataStruct();
  static std::vector< std::vector< std::pair<int, int> > > startEdgeGraphSubgraphs(Hypergraph&, int);
  static std::vector< std::pair<long long, std::string> > getEquivalenceClass(const std::vector< std::vector<int> >&, int);
  static std::vector< std::vector< std::pair<int, int> > > getAllSubgraphs(const std::vector< std::vector<int> >&, int);
  
  static std::map< int, long long> k3(Hypergraph&);
  static std::map< int, long long> k3Fase(Hypergraph&);
  static std::map< int, long long> k3Triangle(Hypergraph& inputGraph);
  static std::map< int, long long> k3ESU(Hypergraph& inputGraph);
  static std::map< int, long long> bruteForce3(Hypergraph& inputGraph);

  static std::map< int, long long> k4(Hypergraph&);
  static std::map< int, long long> k4Fase(Hypergraph&);
  static std::map< int, long long> k4ESU(Hypergraph& inputGraph);
  static std::map< int, long long> bruteForce4(Hypergraph& inputGraph);
  
  static void k3IntermediateForm(Hypergraph& inputGraph);
  static void k4IntermediateForm(Hypergraph& inputGraph);
  
  static void printResults(std::chrono::time_point<std::chrono::steady_clock>, std::chrono::time_point<std::chrono::steady_clock>, std::map< int, long long>, int, bool, std::ostream&);
  
public:

  
  static void networkCensus(Hypergraph&, int, bool, std::ostream&);
  
  static void findMotifs(Hypergraph&, int, bool, std::ostream&);
  
};
