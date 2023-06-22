#pragma once

#include "Settings.hpp"

class IsomorphismHyper {
private:
  static const int MAX_SZ = MAX_INPUT_N; 
  static const int MAX_N = MAX_MOTIF_SIZE;
  static const int MAX_M = SETWORDSNEEDED(MAX_N);
  static int lab[MAX_INPUT_N];
  static int ptn[MAX_INPUT_N];
  static int orbits[MAX_INPUT_N];
  static int counter;
  static graph g[MAX_INPUT_N * MAX_M];
  static graph cg[MAX_INPUT_N * MAX_M];
  static std::map< std::vector<graph>, int> found;
  static set* gv;
  static std::map< std::vector< std::pair<int, int> >, std::string> canonStrCache;
  static std::unordered_map< std::vector< std::vector<int> >, int, HashFunction> canonCache;
  static std::map< int, std::vector< std::vector<int> > > canonCacheReverse;
  static std::vector<graph> canonization(Hypergraph& h);
public:
   static bool isomorphismSlow( Hypergraph& h1,  Hypergraph& h2);
   static bool isomorphismNauty( Hypergraph& h1,  Hypergraph& h2);
   static int getLabel(Hypergraph& h);
   static int getLabel(const std::vector<std::vector<int> >&);
   static std::vector<graph> canonization( std::vector<std::vector<int> >&);
   static std::string canonStr( std::vector<std::pair<int, int> >&, int); 
   static std::vector< std::vector<int> > getHypergraph(int);
   static void precalc(int k);
};
