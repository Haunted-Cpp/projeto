#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>

#include "nauty.h"
#include "Hypergraph.hpp"
#include "Isomorphism.hpp"

using namespace std;


/* 
 * Initialize static variables
 */
 
int Isomorphism::lab[MAX_MOTIF_SIZE];
int Isomorphism::ptn[MAX_MOTIF_SIZE];
int Isomorphism::orbits[MAX_MOTIF_SIZE];
graph Isomorphism::g[MAX_N * MAX_M];
graph Isomorphism::cg[MAX_N * MAX_M];
set* Isomorphism::gv;
map< vector< pair<int, int> >, string> Isomorphism::canonStrCache; // could maybe be replaced by hashing or pre-calc

 
  
bool Isomorphism::isomorphismSlow(const Hypergraph& h1, const Hypergraph& h2) {
  if (h1.getNodeCount() != h2.getNodeCount()) return false;
  if (h1.getEdgeCount() != h2.getEdgeCount()) return false;
  const int nodeCount = h1.getNodeCount();
  assert(nodeCount <= 10); // If value is greater than 10 this method is too slow!
  vector<int> perm(nodeCount);
  iota(perm.begin(), perm.end(), 0);
  do {
    vector< vector<int> > modifiedIncidenceMatrix = h2.applyFunction(perm);
    if (h1.isEqual(modifiedIncidenceMatrix)) return true;
  } while (next_permutation(perm.begin(), perm.end()));
  return false;
}

bool Isomorphism::isomorphismNauty(const Hypergraph& h1, const Hypergraph& h2) {
  return canonization(h1) == canonization(h2);
}


string Isomorphism::canonStr(const vector< pair<int, int> >& edgeList, int n) {
  if (canonStrCache.find(edgeList) != canonStrCache.end()) {
    // this step can be done also with a tree like structure ...
    return canonStrCache[edgeList];
  }
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  options.getcanon = TRUE;
  int m = SETWORDSNEEDED(n);
  //nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
  EMPTYGRAPH(g,m,n);
  for (auto& [a, b] : edgeList) ADDONEEDGE(g, a, b, m);
  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
  string adjMat = "";
  for (int i = 0; i < n; i++) {
    gv = GRAPHROW(cg,i,m); 
    for (int j = 0; j <= i; j++) adjMat += ISELEMENT(gv,j) ? '1' : '0';
    //adjMat += '\n';
  }
  return canonStrCache[edgeList] = adjMat;
}

/*
string Isomorphism::canonStrSparse(const vector< pair<int, int> >& edgeList, int n) {
  //DYNALLSTAT(int,lab,lab_sz);
  //DYNALLSTAT(int,ptn,ptn_sz);
  //DYNALLSTAT(int,orbits,orbits_sz);
  static DEFAULTOPTIONS_SPARSEGRAPH(options);
  statsblk stats;
  sparsegraph sg; 
  int n,m,i;
  options.getcanon = TRUE;
  SG_INIT(sg);
  m = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
  //DYNALLOC1(int,lab,lab_sz,n,"malloc");
  //DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
  //DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
  SG_ALLOC(sg,n,2*n,"malloc");
  sg.nv = n; 
  sg.nde = 2*n; 
  for (i = 0; i < n; ++i) {
    sg.v[i] = 2*i;
    sg.d[i] = 2;
    sg.e[2*i] = (i+n-1)%n;
    sg.e[2*i+1] = (i+n+1)%n; 
  }
  //printf("Generators for Aut(C[%d]):\n",n);
  sparsenauty(&sg,lab,ptn,orbits,&options,&stats,cg);
  //printf("Automorphism group size = ");
  //writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
  //printf("\n");
  string adjMat = "";
  for (int i = 0; i < n; i++) {
    gv = GRAPHROW(cg,i,m); 
    for (int j = 0; j <= i; j++) adjMat += ISELEMENT(gv,j) ? '1' : '0';
    adjMat += '\n';
  }
  return adjMat;

}*/

vector<graph> Isomorphism::canonization(const Hypergraph& h) {
  DYNALLSTAT(int,lab,lab_sz);
  DYNALLSTAT(int,ptn,ptn_sz);
  DYNALLSTAT(int,orbits,orbits_sz);
  DYNALLSTAT(graph,g,g_sz);
  DYNALLSTAT(graph,cg,cg_sz);
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  int n = h.getNodeCount() + h.getEdgeCount(); // number of nodes in our modified graph # node + # edges
  /* Select option for canonical labelling */
  options.getcanon = TRUE;
  /* Select option to use node colours */
  options.defaultptn = FALSE;
  int m = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
  DYNALLOC1(int,lab,lab_sz,n,"malloc");
  DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
  DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
  DYNALLOC2(graph,g,g_sz,n,m,"malloc");
  EMPTYGRAPH(g,m,n);
  vector< vector<int> > incidenceMatrix = h.getIncidenceMatrix();
  for (int i = 0; i < h.getEdgeCount(); i++) {
    for (auto& node : incidenceMatrix[i]) {
      ADDONEEDGE(g, node, h.getNodeCount() + i, m);
    }
  }
  /* Add Colors */
  for (int i = 0; i < n; i++) {
    lab[i] = i;
    ptn[i] = 1;
  }
  ptn[h.getNodeCount() - 1] = 0;
  ///* Create canonical graph */
  DYNALLOC2(graph,cg,cg_sz,n,m,"malloc");
  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
  vector<graph> labels;
  for (int i = 0; i < n; i++) labels.emplace_back(cg[i]);
  return labels;
}
