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
 
 
int IsomorphismHyper::no_use = 0;
 
int IsomorphismHyper::lab[MAX_INPUT_N];
int IsomorphismHyper::ptn[MAX_INPUT_N];
int IsomorphismHyper::orbits[MAX_INPUT_N];
int IsomorphismHyper::counter = 0; // used to precalc every hash 
graph IsomorphismHyper::g[MAX_INPUT_N * MAX_M];
graph IsomorphismHyper::cg[MAX_INPUT_N * MAX_M];
set* IsomorphismHyper::gv;

map< vector< pair<int, int> >, string> IsomorphismHyper::canonStrCache; // could maybe be replaced by hashing or pre-calc

std::unordered_map< vector< vector<int> >, int, HashFunction> IsomorphismHyper::canonCache;

std::map< vector<graph>, int> IsomorphismHyper::found;

map< int, vector< vector<int> > > IsomorphismHyper::canonCacheReverse;


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
  assert(no_use == 0); // REMOVEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  // Make sure this is only called in the beginning ..............
  
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
  return labels;
}

int IsomorphismHyper::getLabel(Hypergraph& h) {
  auto adj = h.getIncidenceMatrix();
  assert(is_sorted(adj.begin(), adj.end()));
  if (canonCache.find(h.getIncidenceMatrix()) != canonCache.end()) return canonCache[h.getIncidenceMatrix()];
  assert(false);
}


void IsomorphismHyper::precalc(int N) {
  vector< vector<int> > edges;
  for (int mask = 0; mask < (1 << N); mask++) {
    if (__builtin_popcount(mask) < 2) continue;
    vector<int> edge;
    for (int i = 0; i < N; i++) {
      if ((mask >> i) & 1) {
        edge.emplace_back(i);
      }
    }
    edges.emplace_back(edge);
  } 
  const int M = (int) edges.size();
  for (int mask = 1; mask < (1 << M); mask++) {
    vector< vector<int> > adj;
    vector< vector<int> > g;
    for (int i = 0; i < M; i++) {
      if ((mask >> i) & 1) {
        adj.emplace_back(edges[i]);
        for (int e1 = 0; e1 < edges[i].size(); e1++) {
          for (int e2 = e1 + 1; e2 < edges[i].size(); e2++) {
            g.push_back({edges[i][e1], edges[i][e2]});
          }
        }
      }
    }
    sort(g.begin(), g.end());
    g.erase(unique(g.begin(), g.end()), g.end());
    Hypergraph h;
    h.setN(N);
    h.setIncidenceMatrix(g);
    if (!h.is_two_connected()) continue;
    Hypergraph h1;
    h1.setN(N);
    h1.setIncidenceMatrix(adj);
    vector<graph> label = canonization(h1);
    sort(adj.begin(), adj.end());
    if (found.find(label) != found.end()) {
      int id = found[label];
      canonCache[adj] = id;
      continue;
    } 
    canonCache[adj] = counter;
    canonCacheReverse[counter] = adj;
    found[label] = counter;
    ++counter;
  }
}

vector< vector<int> > IsomorphismHyper::getHypergraph(int mask) {
  assert(canonCacheReverse.find(mask) != canonCacheReverse.end());
  return canonCacheReverse[mask];
}
