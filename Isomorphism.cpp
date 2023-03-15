#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>

#include "nauty.h"
#include "Hypergraph.hpp"
#include "Isomorphism.hpp"

using namespace std;

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
