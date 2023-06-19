#include <bits/stdc++.h>

#include "nauty.h"
#include "Hypergraph.hpp"
#include "IsomorphismHyper.hpp"

using std::string;
using std::vector;
using std::pair;
using std::map;
using std::vector;

/* 
 * Initialize static variables
 */
 
int IsomorphismHyper::lab[MAX_INPUT_N];
int IsomorphismHyper::ptn[MAX_INPUT_N];
int IsomorphismHyper::orbits[MAX_INPUT_N];
graph IsomorphismHyper::g[MAX_INPUT_N * MAX_M];
graph IsomorphismHyper::cg[MAX_INPUT_N * MAX_M];
set* IsomorphismHyper::gv;

map< vector< pair<int, int> >, string> IsomorphismHyper::canonStrCache; // could maybe be replaced by hashing or pre-calc

map< vector< vector<int> >, vector<graph> > IsomorphismHyper::canonCache;
map< vector<graph>, vector< vector<int> > > IsomorphismHyper::canonCacheReverse;


bool IsomorphismHyper::isomorphismSlow(Hypergraph& h1, Hypergraph& h2) {
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

bool IsomorphismHyper::isomorphismNauty(Hypergraph& h1, Hypergraph& h2) {
  return canonization(h1) == canonization(h2);
}

string IsomorphismHyper::canonStr( vector< pair<int, int> >& edgeList, int n) {
  if (canonStrCache.find(edgeList) != canonStrCache.end()) return canonStrCache[edgeList];
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  options.getcanon = TRUE;
  int m = SETWORDSNEEDED(n);
  EMPTYGRAPH(g,m,n);
  for (auto& [a, b] : edgeList) ADDONEEDGE(g, a, b, m);
  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
  string adjMat = "";
  for (int i = 0; i < n; i++) {
    gv = GRAPHROW(cg,i,m); 
    for (int j = 0; j <= i; j++) adjMat += ISELEMENT(gv,j) ? '1' : '0';
  }
  return canonStrCache[edgeList] = adjMat;
}


vector<graph> IsomorphismHyper::canonization(vector< vector<int> >& adj) {
  // CAREFUL ! CACHE CANNOT BE USED DIRECTLY
  //if (canonCache.find(adj) != canonCache.end()) return canonCache[adj];
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  
  int n = (int) adj.size(); // number of nodes in our modified graph # node + # edges
  int deg = 0;
  for (auto edge : adj) deg += (int) edge.size();
  deg /= 2;
  n += deg;
  
  /* Select option for canonical labelling */
  options.getcanon = TRUE;
  /* Select option to use node colours */
  options.defaultptn = FALSE;
  int m = SETWORDSNEEDED(n);
  EMPTYGRAPH(g,m,n);
  
  vector< vector<int> > incidenceMatrix;
  for (int i = 0; i < (int) adj.size(); i++) {
    for (auto& nei : adj[i]) {
      assert(i != nei); // no self-loop
      if (i < nei) {
        incidenceMatrix.push_back({i, nei});
      }
    }
  }
  assert( (int) incidenceMatrix.size() == deg );
  for (int i = 0; i < incidenceMatrix.size(); i++) {
    for (auto& node : incidenceMatrix[i]) {
      ADDONEEDGE(g, node, (int) adj.size() + i, m); // --- nodes SHOULD be numbered from 0 to n - 1 !!!
    }
  }
  /* Add Colors */
  for (int i = 0; i < n; i++) {
    lab[i] = i;
    ptn[i] = 1;
  }
  ptn[(int) adj.size() - 1] = 0;
  ///* Create canonical graph */
  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
  vector<graph> labels;
  for (int i = 0; i < n; i++) labels.emplace_back(cg[i]);
  return labels;
}

vector<graph> IsomorphismHyper::canonization(Hypergraph& h) {
  if (canonCache.find(h.getIncidenceMatrix()) != canonCache.end()) return canonCache[h.getIncidenceMatrix()];
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  int n = h.getNodeCount() + h.getEdgeCount(); // number of nodes in our modified graph # node + # edges
  /* Select option for canonical labelling */
  options.getcanon = TRUE;
  /* Select option to use node colours */
  options.defaultptn = FALSE;
  int m = SETWORDSNEEDED(n);
  EMPTYGRAPH(g,m,n);
  vector< vector<int> > incidenceMatrix = h.getIncidenceMatrix();
  for (int i = 0; i < h.getEdgeCount(); i++) {
    for (auto& node : incidenceMatrix[i]) {
      ADDONEEDGE(g, node, h.getNodeCount() + i, m); // --- nodes SHOULD be numbered from 0 to n - 1 !!!
    }
  }
  /* Add Colors */
  for (int i = 0; i < n; i++) {
    lab[i] = i;
    ptn[i] = 1;
  }
  ptn[h.getNodeCount() - 1] = 0;
  ///* Create canonical graph */
  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
  vector<graph> labels;
  for (int i = 0; i < n; i++) labels.emplace_back(cg[i]);
  
  
  canonCacheReverse[labels] = h.getIncidenceMatrix();
  return canonCache[h.getIncidenceMatrix()] = labels;
}

vector< vector<int> > IsomorphismHyper::getHypergraph(const vector<graph>& label) {
  assert(canonCacheReverse.find(label) != canonCacheReverse.end());
  return canonCacheReverse[label];
}
