#pragma once

#include "Settings.hpp"

class Isomorphism {
private:
  static const int MAX_SZ = MAX_INPUT_N; 
  static const int MAX_N = MAX_MOTIF_SIZE;
  static const int MAX_M = SETWORDSNEEDED(MAX_N);
  
  
  static int lab[MAX_INPUT_N];
  static int ptn[MAX_INPUT_N];
  static int orbits[MAX_INPUT_N];
  static graph g[MAX_INPUT_N * MAX_M];
  static graph cg[MAX_INPUT_N * MAX_M];
  
  
  static set* gv;
  static std::map< std::vector< std::pair<int, int> >, std::string> canonStrCache;
  static std::map< std::vector< std::vector<int> >, std::vector<graph> > canonCache;
  static std::map< std::vector<graph>, std::vector< std::vector<int> > > canonCacheReverse;
public:
   static bool isomorphismSlow( Hypergraph& h1,  Hypergraph& h2);
   static bool isomorphismNauty( Hypergraph& h1,  Hypergraph& h2);
   static std::vector<graph> canonization( Hypergraph& h);
   static std::vector<graph> canonization( std::vector<std::vector<int> >&);
   
   static std::string canonStr( std::vector<std::pair<int, int> >&, int); 
   static std::string canonStrSparse( std::vector<std::pair<int, int> >&, int); 
   static std::vector< std::vector<int> > getHypergraph(const std::vector<graph>&);

};
