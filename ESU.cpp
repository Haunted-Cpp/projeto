#include <bits/stdc++.h>

#include "Hypergraph.hpp"
#include "ESU.hpp"
#include "nauty.h"
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
vector< vector<int> > ESU::graph;
vector< vector<int> > ESU::removeExt;
vector< vector<int> > ESU::removeNei;
int ESU::K; // Subgraph size
int ESU::V; // Initial node
bool ESU::computeEquivalenceClass; // Compress every node into a single class?

map<string, int> ESU::counter;
vector< vector< pair<int, int> > > ESU::subgraphs;

void ESU::enumerateSubgraphs() {
  if ((int) subgraph.size() == K) {
    if (computeEquivalenceClass) {
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
    for (int i = 0; i < (int) graph[w].size(); i++) {
      int u = graph[w][i];
      if (subgraph.find(u) != subgraph.end()) {
        if (computeEquivalenceClass) edgeList.emplace_back(pos[w], pos[u]);
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
    for (int i = 0; i < (int) graph[w].size(); i++) {
      if (removeNei[w][i]) neiSubgraph.erase(graph[w][i]);
      if (removeExt[w][i]) extension.erase(graph[w][i]);
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

void ESU::setupAndRun(const vector< vector<int> >& g, int k, bool merge) {
  subgraphs.clear();
  computeEquivalenceClass = merge;
  counter.clear();
  setGraph(g);
  clearDataStruct();
  K = k;
  for (int i = 0; i < getNodeCount(); i++) {
    clearDataStruct();
    V = i;
    pos[V] = 0;
    for (auto& u : graph[i]) {
      if (u > V) extension.insert(u);
      neiSubgraph.insert(u);
    }
    subgraph.insert(V);
    enumerateSubgraphs();
  }
}

vector< pair<int, string> > ESU::getEquivalenceClass(const vector< vector<int> >& g, int k) {
  setupAndRun(g, k, true);
  vector< pair<int, string> > census;
  for (auto [subgraph, cnt] : counter) census.emplace_back(cnt, subgraph);
  sort(census.rbegin(), census.rend());
  return census;
}

vector< vector< pair<int, int> > > ESU::getAllSubgraphs(const vector< vector<int> >& g, int k) {
  setupAndRun(g, k, false);
  return subgraphs;
}

vector< vector< pair<int, int> > > ESU::startEdgeGraphSubgraphs(Hypergraph& h, const int k) {
  return getAllSubgraphs(h.buildEdgeGraph(), k);
}

void ESU::setGraph(const vector< vector<int> >& inputGraph) {
  graph = inputGraph;
  removeExt.resize(getNodeCount());
  removeNei.resize(getNodeCount());
  for (int i = 0; i < getNodeCount(); i++) {
    removeExt[i].resize((int)graph[i].size());
    removeNei[i].resize((int)graph[i].size());
    fill(removeExt[i].begin(), removeExt[i].end(), 0);
    fill(removeNei[i].begin(), removeNei[i].end(), 0);
  }
}

int ESU::getNodeCount() {
  return (int) graph.size();
}

// implementar o caso K = 3 do paper


// queria testar a minha variação mas com limite de nós -> problema (pode ficar com menos ...) -> adicionar enquanto for menor ... -> pode ter mts arestas!!!


// ESU Hypergraph version
