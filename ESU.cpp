#include <bits/stdc++.h>

#include "Hypergraph.hpp"
#include "nauty.h"
#include "ESU.hpp"
#include "Isomorphism.hpp"

using namespace std;

/* 
 * Initialize static variables
 */

int ESU::pos[MAX_INPUT_N];
vector< std::pair<int, int> > ESU::edgeList;
std::set<int> ESU::subgraph;
std::set<int> ESU::extension;
std::set<int> ESU::neiSubgraph;
stack< std::set<int> > ESU::saveExtension;
vector< vector<int> > ESU::removeExt;
vector< vector<int> > ESU::removeNei;
int ESU::K; // Subgraph size
int ESU::V; // Initial node

enum option { COUNT_ALL, CLASS_ONLY, HYPERGRAPH};

// Used in classical graph
vector< vector<int> > ESU::g;
map<string, int> ESU::counter;
vector< vector< pair<int, int> > > ESU::subgraphs;
// 

// Used in Hypergraph
Hypergraph ESU::h;
map< vector<graph>, int> ESU::counterHyper;
std::set< vector<int> > ESU::visited;
option Search;
//

//8
//4327
//958
//801
//2652
//1346
//208
//8
//18


void ESU::enumerateSubgraphs() {
  if ((int) subgraph.size() == K) {
    if (Search == HYPERGRAPH) {
      vector<int> nodes;
      nodes.assign(subgraph.begin(), subgraph.end());
      assert(is_sorted(nodes.begin(), nodes.end()));
      Hypergraph motif = h.induceSubgraph(nodes);
      if (motif.getEdgeMaxDeg() != 2) return;
      assert(visited.find(nodes) == visited.end());
      assert(motif.getNodeCount() == 4);
      counterHyper[Isomorphism::canonization(motif)]++;
      return;
    } else if (Search == CLASS_ONLY) {
      // Join all topological equivalent subgraphs into a single class
      string mat = Isomorphism::canonStr(edgeList, K);
      ++counter[mat]; // we could use a trie like structure here
    } else {
      // Just add the subgraph created
      subgraphs.emplace_back(edgeList);
    }
    return;
  }
  saveExtension.push(extension);
  while (!extension.empty()) {
    int w = *extension.begin();
    extension.erase(extension.begin());
    pos[w] = (int) subgraph.size();
    int added = 0;
    for (int i = 0; i < (int) g[w].size(); i++) {
      int u = g[w][i];
      if (subgraph.find(u) != subgraph.end()) {
        if (Search == CLASS_ONLY) edgeList.emplace_back(pos[w], pos[u]);
        else edgeList.emplace_back(w, u);
        ++added;
      } else if (u > V && neiSubgraph.find(u) == neiSubgraph.end()) {
         extension.insert(u);
         removeExt[w][i] = 1;
      } 
      if (neiSubgraph.find(u) == neiSubgraph.end()) {
        neiSubgraph.insert(u);
        removeNei[w][i] = 1;
      }
    }
    subgraph.insert(w);
    enumerateSubgraphs(); // recursively add more nodes
    subgraph.erase(w); // restore data struct state
    while (added--) edgeList.pop_back();
    for (int i = 0; i < (int) g[w].size(); i++) {
      if (removeNei[w][i]) neiSubgraph.erase(g[w][i]);
      if (removeExt[w][i]) extension.erase(g[w][i]);
      removeExt[w][i] = removeNei[w][i] = 0;
    }
  }
  extension = saveExtension.top();
  saveExtension.pop();
};

void ESU::clearDataStruct() {
  edgeList.clear();
  subgraph.clear();
  extension.clear();
  neiSubgraph.clear();
}

void ESU::setupAndRun(const vector< vector<int> >& g, int k) {
  subgraphs.clear();
  //computeEquivalenceClass = merge;
  counter.clear();
  setGraph(g);
  clearDataStruct();
  K = k;
  for (int i = 0; i < getNodeCount(); i++) {
    clearDataStruct();
    V = i;
    pos[V] = 0;
    for (auto& u : g[i]) {
      if (u > V) extension.insert(u);
      neiSubgraph.insert(u);
    }
    subgraph.insert(V);
    enumerateSubgraphs();
  }
}

vector< pair<int, string> > ESU::getEquivalenceClass(const vector< vector<int> >& g, int k) {
  Search = CLASS_ONLY;
  setupAndRun(g, k);
  vector< pair<int, string> > census;
  for (auto [subgraph, cnt] : counter) census.emplace_back(cnt, subgraph);
  sort(census.rbegin(), census.rend());
  return census;
}

vector< vector< pair<int, int> > > ESU::getAllSubgraphs(const vector< vector<int> >& g, int k) {
  Search = COUNT_ALL;
  setupAndRun(g, k);
  return subgraphs;
}

vector< vector< pair<int, int> > > ESU::startEdgeGraphSubgraphs(Hypergraph& h, const int k) {
  return getAllSubgraphs(h.buildEdgeGraph(), k);
}

void ESU::setGraph(const vector< vector<int> >& inputGraph) {
  g = inputGraph;
  removeExt.resize(getNodeCount());
  removeNei.resize(getNodeCount());
  for (int i = 0; i < getNodeCount(); i++) {
    removeExt[i].resize((int)g[i].size());
    removeNei[i].resize((int)g[i].size());
    fill(removeExt[i].begin(), removeExt[i].end(), 0);
    fill(removeNei[i].begin(), removeNei[i].end(), 0);
  }
}

int ESU::getNodeCount() {
  return (int) g.size();
}

// implementar o caso K = 3 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

void ESU::k3(Hypergraph& inputGraph) {
  counterHyper.clear();
  visited.clear();
  for (auto edge : inputGraph.getIncidenceMatrix()) {
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    counterHyper[Isomorphism::canonization(motif)]++;
    visited.insert(edge);
  }
  vector< vector<int> > filter( inputGraph.getNodeCount() );
  for (auto& edge : inputGraph.filterEdge(2).getIncidenceMatrix()) {
    filter[ edge[0] ].emplace_back( edge[1] );
    filter[ edge[1] ].emplace_back( edge[0] );
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  setupAndRun(filter, 3);
  cout << counterHyper.size() << '\n';
  for (auto [x, cnt] : counterHyper) {
    cout << cnt << '\n';
  }
}

// implementar o caso K = 4 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf
void ESU::k4(Hypergraph& inputGraph) {
  counterHyper.clear();
  visited.clear();
  for (auto edge : inputGraph.getIncidenceMatrix()) {
    if (edge.size() != 4) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    counterHyper[Isomorphism::canonization(motif)]++;
    assert(visited.find(edge) == visited.end()); // assuming no duplicate edges
    //visited.insert(edge);
    assert(motif.getEdgeMaxDeg() == 4);
  }
  Hypergraph reducedGraph = inputGraph.filterEdge(3); // at most 3 edges
  vector< vector<int> > edgeGraph = reducedGraph.buildEdgeGraph(); // I don't like this ... It doesn't seem needed!
  for (int e = 0; e < (int) edgeGraph.size(); e++) {
    for (auto e1 : edgeGraph[e]) {
      vector<int> nodes;
      vector<int> ve = reducedGraph.getEdge(e);
      if ( (int) ve.size() != 3 ) continue;
      vector<int> ve1 = reducedGraph.getEdge(e1);
      nodes.insert(nodes.end(), ve1.begin(), ve1.end());
      nodes.insert(nodes.end(), ve.begin(), ve.end());
      sort(nodes.begin(), nodes.end());
      nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
      if ( (int) nodes.size() != 4) continue;
      if (visited.find(ve1) != visited.end()) continue;
      Hypergraph motif = inputGraph.induceSubgraph(nodes);
      if (motif.getEdgeMaxDeg() != 3) continue;
      assert(motif.getEdgeMaxDeg() == 3);
      counterHyper[Isomorphism::canonization(motif)]++;
      visited.insert(ve);
      assert(is_sorted(ve.begin(), ve.end()));
    }
  }
  vector< vector<int> > filter(inputGraph.getNodeCount());
  for (auto& edge : reducedGraph.filterEdge(2).getIncidenceMatrix()) {
    filter[ edge[0] ].emplace_back( edge[1] );
    filter[ edge[1] ].emplace_back( edge[0] );
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  setupAndRun(filter, 4);
  cout << counterHyper.size() << '\n';
  vector<int> xx;
  for (auto [x, cnt] : counterHyper) {
    //cout << cnt << '\n';
    xx.emplace_back(cnt);
  }
  sort(xx.rbegin(), xx.rend());
  for (auto val : xx) cout << val << '\n';
}











// Experimental code, ESU with GTRIE
// implementar o caso K = 3 do paper








/*

void ESU::enumerateSubgraphs() {
  if ((int) subgraph.size() == K) {
    // We are checking Hypergraphs here, only
    vector<int> nodes;
    nodes.assign(subgraph.begin(), subgraph.end());
    if (visited.find(nodes) != visited.end()) return;
    Hypergraph motif = h.induceSubgraph(nodes);
    assert(motif.getNodeCount() == 4);
    counterHyper[Isomorphism::canonization(motif)]++;
    return;
  }
  saveExtension.push(extension);
  while (!extension.empty()) {
    int w = *extension.begin();
    extension.erase(extension.begin());
    pos[w] = (int) subgraph.size();
    int added = 0;
    for (int i = 0; i < (int) g[w].size(); i++) {
      int u = g[w][i];
      if (subgraph.find(u) != subgraph.end()) {
        if (Search == CLASS_ONLY) edgeList.emplace_back(pos[w], pos[u]);
        else edgeList.emplace_back(w, u);
        ++added;
      } else if (u > V && neiSubgraph.find(u) == neiSubgraph.end()) {
         extension.insert(u);
         removeExt[w][i] = 1;
      } 
      if (neiSubgraph.find(u) == neiSubgraph.end()) {
        neiSubgraph.insert(u);
        removeNei[w][i] = 1;
      }
    }
    subgraph.insert(w);
    enumerateSubgraphs(); // recursively add more nodes
    subgraph.erase(w); // restore data struct state
    while (added--) edgeList.pop_back();
    for (int i = 0; i < (int) g[w].size(); i++) {
      if (removeNei[w][i]) neiSubgraph.erase(g[w][i]);
      if (removeExt[w][i]) extension.erase(g[w][i]);
      removeExt[w][i] = removeNei[w][i] = 0;
    }
  }
  extension = saveExtension.top();
  saveExtension.pop();
};

void ESU::setupAndRunGtrie(const vector< vector<int> >& g, int k) {
  subgraphs.clear();
  //computeEquivalenceClass = merge;
  counter.clear();
  setGraph(g);
  clearDataStruct();
  K = k;
  for (int i = 0; i < getNodeCount(); i++) {
    clearDataStruct();
    V = i;
    pos[V] = 0;
    for (auto& u : g[i]) {
      if (u > V) extension.insert(u);
      neiSubgraph.insert(u);
    }
    subgraph.insert(V);
    enumerateSubgraphsGTrie();
  }
}



void ESU::k3GTrie(Hypergraph& inputGraph) {
  counterHyper.clear();
  visited.clear();
  for (auto edge : inputGraph.getIncidenceMatrix()) {
    if (edge.size() != 3) continue;
    Hypergraph motif = h.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    //counterHyper[Isomorphism::canonization(motif)]++;
    visited.insert(edge); // notice how it's bounded by the number of edges ...
  }
  vector< vector<int> > filter( inputGraph.getNodeCount() );
  for (auto& edge : inputGraph.filterEdge(2).getIncidenceMatrix()) {
    filter[ edge[0] ].emplace_back( edge[1] );
    filter[ edge[1] ].emplace_back( edge[0] );
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  setupAndRun(filter, 3);
  cout << counterHyper.size() << '\n';
  for (auto [x, cnt] : counterHyper) {
    cout << cnt << '\n';
  }
}
*/
