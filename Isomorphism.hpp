#define MAX_MOTIF_SIZE 1000
#define MAX_INPUT_N 1000

#pragma once


class Isomorphism {
private:
  static const int MAX_SZ = MAX_INPUT_N; 
  static const int MAX_N = MAX_MOTIF_SIZE;
  static const int MAX_M = SETWORDSNEEDED(MAX_N);
  
  //static int lab[MAX_INPUT_N];
  //static int ptn[MAX_INPUT_N];
  //static int orbits[MAX_INPUT_N];
  //static graph g[MAX_INPUT_N * MAX_INPUT_N];
  //static graph cg[MAX_INPUT_N * MAX_INPUT_N];
  
  static int mapping[MAX_INPUT_N];
  
  //static int lab[MAXN];
  //DYNALLSTAT(int,lab,lab_sz);
    //DYNALLSTAT(int,ptn,ptn_sz);
    //DYNALLSTAT(int,orbits,orbits_sz);
    //DYNALLSTAT(graph,g,g_sz);
    //DYNALLSTAT(graph,cg,cg_sz);
  //static void compress(std::vector<std::pair<int, int> >&);
public:
  
   static bool isomorphismSlow(const Hypergraph& h1, const Hypergraph& h2);
   static bool isomorphismNauty(const Hypergraph& h1, const Hypergraph& h2);
   static std::vector<graph> canonization(const Hypergraph& h);
   static std::string canonStr(const std::vector<std::pair<int, int> >&, int); 
};
