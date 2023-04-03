#pragma once

#include "Settings.hpp"

class Isomorphism {
private:
  static const int MAX_SZ = MAX_INPUT_N; 
  static const int MAX_N = MAX_MOTIF_SIZE;
  static const int MAX_M = SETWORDSNEEDED(MAX_N);
  static int lab[MAX_MOTIF_SIZE];
  static int ptn[MAX_MOTIF_SIZE];
  static int orbits[MAX_MOTIF_SIZE];
  static graph g[MAX_N * MAX_M];
  static graph cg[MAX_N * MAX_M];
  static set* gv;
  static std::map< std::vector< std::pair<int, int> >, std::string> canonStrCache;
public:
   static bool isomorphismSlow(const Hypergraph& h1, const Hypergraph& h2);
   static bool isomorphismNauty(const Hypergraph& h1, const Hypergraph& h2);
   static std::vector<graph> canonization(const Hypergraph& h);
   static std::string canonStr(const std::vector<std::pair<int, int> >&, int); 
   static std::string canonStrSparse(const std::vector<std::pair<int, int> >&, int); 
};
