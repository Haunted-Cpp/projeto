#include <bits/stdc++.h>

#include "Hypergraph.hpp"
#include "nauty.h"
#include "ESU.hpp"
#include "Isomorphism.hpp"
#include "Settings.hpp"

/* 
 * Initialize static variables
 */


int ESU::f[MAX_INPUT_N];
int ESU::pos[MAX_INPUT_N];

vector< std::pair<int, int> > ESU::edgeList;

std::vector<int> ESU::subgraph;
std::vector<int> ESU::subgraph_compressed;

int ESU::K; // Subgraph size
int ESU::V; // Initial node

enum option { COUNT_ALL, CLASS_ONLY, HYPERGRAPH};

// Used in classical graph
vector< vector<int> > ESU::g;
map<string, long long> ESU::counter;
vector< vector< pair<int, int> > > ESU::subgraphs;
// 

// Used in Hypergraph
Hypergraph ESU::h;
map< vector<graph>, long long> ESU::counterHyper;
std::unordered_set< vector<int>, HashFunction> ESU::visited;
option Search;


vector<int> next_extension;

void ESU::enumerateSubgraphs(vector<int> extension) {
  if ( (int) subgraph.size() == K) {
    if (Search == HYPERGRAPH) {
      vector<int> subgraph_ordered = subgraph;
      sort(subgraph_ordered.begin(), subgraph_ordered.end());
      if (visited.find(subgraph_ordered) != visited.end()) {
        return;
      }
      Hypergraph motif = h.induceSubgraphNoComp(subgraph_ordered, subgraph_compressed);
      //Hypergraph motif = h.induceSubgraph(subgraph_ordered);
      counterHyper[Isomorphism::canonization(motif)]++;
    } else if (Search == CLASS_ONLY) {
      string mat = Isomorphism::canonStr(edgeList, K);
      ++counter[mat]; // we could use a trie like structure here
    } else {
      // Just add the subgraph created
      subgraphs.emplace_back(edgeList);
    }
    return;
  }
  while (!extension.empty()) {
    int current_node = extension.back();
    extension.pop_back();
    next_extension = extension;
    pos[current_node] = (int) subgraph.size();
    int added_nodes = 0;
    for (auto& to : g[current_node]) {
      if (to > V) {
        if (find(subgraph.begin(), subgraph.end(), to) != subgraph.end()) {
          ++added_nodes;
          if (Search != HYPERGRAPH) {
            if (Search == CLASS_ONLY) edgeList.emplace_back(pos[current_node], pos[to]);
            else edgeList.emplace_back(current_node, to);
          }
        } else if (f[to] == 0) next_extension.emplace_back(to);
        f[to]++;
      }
    }
    subgraph.emplace_back(current_node);
    subgraph_compressed.emplace_back(pos[current_node]);
    enumerateSubgraphs(next_extension);
    subgraph.pop_back();
    subgraph_compressed.pop_back();
    while (added_nodes--) edgeList.pop_back();
    for (auto& to : g[current_node]) if (to > V) --f[to];
  }
};

void ESU::clearDataStruct() {
  edgeList.clear();
  subgraph.clear();
}

void ESU::setupAndRun(const vector< vector<int> >& inputGraph, int k) {
  subgraphs.clear();
  //computeEquivalenceClass = merge;
  counter.clear();
  clearDataStruct();
  K = k;
  g = inputGraph;
  vector<int> extension;
  for (int i = 0; i < inputGraph.size(); i++) {
    V = i;
    clearDataStruct();
    extension = {i};
    enumerateSubgraphs(extension);
  }
}

vector< pair<long long, string> > ESU::getEquivalenceClass(const vector< vector<int> >& g, int k) {
  Search = CLASS_ONLY;
  setupAndRun(g, k);
  vector< pair<long long, string> > census;
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

// implementar o caso K = 3 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

std::map< std::vector<graph>, long long> ESU::k3(Hypergraph& inputGraph) {
  counterHyper.clear();
  visited.clear();
  for (auto edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    visited.insert(edge); // it will insert the 3 nodes just visited
    counterHyper[Isomorphism::canonization(motif)]++;
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  setupAndRun(h.getGraph(), 3);
  return counterHyper;
}

// implementar o caso K = 4 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

std::map< std::vector<graph>, long long> ESU::k4(Hypergraph& inputGraph) {
  counterHyper.clear();
  visited.clear();
  for (auto& edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 4) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    if (motif.is_two_connected()) {
      continue; // this occurence will be found by ESU!
    }
    counterHyper[Isomorphism::canonization(motif)]++;
    assert(motif.getEdgeMaxDeg() == 4);
  }
  Hypergraph reducedGraph = inputGraph.filterEdge(3); // at most 3 edges
  vector< vector<int> > nei = reducedGraph.buildVertexGraph(3);
  for (auto& edge : reducedGraph.getIncidenceMatrix()) {
    if (edge.size() != 3) continue;
    assert(is_sorted(edge.begin(), edge.end()));
    for (auto& node : edge) { // at most 3
      for (auto& add : nei[node]) { // can be many ... "node" to add 
        if (find(edge.begin(), edge.end(), add) != edge.end()) continue; // it must be a new node
        for (int mask = 1; mask < (1 << 3) - 1; mask++) { // all subsets of current edge
          if ( (mask >> node) & 1 == 0 ) continue; // must include the "node"
          vector<int> new_edge;
          for (int i = 0; i < 3; i++) {
            if ((mask >> i) & 1) {
              new_edge.emplace_back(edge[i]);
            }
          }
          new_edge.emplace_back(add);
          sort(new_edge.begin(), new_edge.end());
          if (! reducedGraph.validEdge(new_edge) ) continue;
          vector<int> nodes = edge;
          nodes.emplace_back(add);
          sort(nodes.begin(), nodes.end());
          if ( (int) nodes.size() != 4) continue;
          Hypergraph motif = inputGraph.induceSubgraph(nodes);
          if (motif.getEdgeMaxDeg() != 3) continue;
          if (visited.find(nodes) != visited.end()) {
            continue;
          }
          assert(motif.getEdgeMaxDeg() == 3);
          if (motif.is_two_connected()) {
            continue; // this occurence will be found later by ESU!
          }
          assert(is_sorted(nodes.begin(), nodes.end()));
          counterHyper[Isomorphism::canonization(motif)]++;
          visited.insert(nodes); // quadratic memory ... i don't like it
        }
      }
    }
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  setupAndRun(h.getGraph(), 4);
  return counterHyper;
}

// This method is used just for debug only!






/**
 * Description: Disjoint Set Union
 * Implementation Notes: 
   * DSU w/ Union by size
   * Find - O ( Log N )
   * Union - O ( 1 )
 * Verification: 
   * https://open.kattis.com/problems/unionfind
*/

struct DisjointSet {
  vector<int> parent, tamanho;
  DisjointSet (int n) {
    parent.clear(); parent.resize(n);
    tamanho.clear(); tamanho.resize(n);
    for (int i = 0; i < n; i++) {
      parent[i] = i;
      tamanho[i] = 1;
    }
  }
  int root (int x) {
    while (x != parent[x]) x = parent[x];
    return x;
  }
  void join (int a, int b) {
    a = root (a);
    b = root (b);
    if (a == b) return;
    if (tamanho[a] < tamanho[b]) swap (a, b);
    parent[b] = a;
    tamanho[a] += tamanho[b];
  }
  bool is_connected (int a, int b) {
    return root(a) == root(b);
  }
};

std::map< std::vector<graph>, long long> ESU::bruteForce3(Hypergraph& inputGraph) {
  
  counterHyper.clear();
  const int n = inputGraph.getNodeCount();
  vector<int> nodes;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      for (int z = j + 1; z < n; z++) {
        nodes = {i, j, z};
        Hypergraph motif = inputGraph.induceSubgraph(nodes);
        DisjointSet dsu(3);
        
        for (auto e : motif.getIncidenceMatrix()) {
          for (int a = 0; a < (int) e.size(); a++) {
            for (int b = a + 1; b < (int) e.size(); b++) {
              dsu.join(e[a], e[b]);
            }
          }
        }
        if (dsu.tamanho[dsu.root(0)] != 3) continue;
        counterHyper[Isomorphism::canonization(motif)]++;
      }
    }
  }
  return counterHyper;
}

 //This method is used just for debug only!
std::map< std::vector<graph>, long long> ESU::bruteForce4(Hypergraph& inputGraph) {
  counterHyper.clear();
  const int n = inputGraph.getNodeCount();
  vector<int> nodes;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      for (int z = j + 1; z < n; z++) {
        for (int x = z + 1; x < n; x++) {
          nodes = {i, j, z, x};
          Hypergraph motif = inputGraph.induceSubgraph(nodes);
          DisjointSet dsu(4);
          for (auto e : motif.getIncidenceMatrix()) {
            for (int a = 0; a < (int) e.size(); a++) {
              for (int b = a + 1; b < (int) e.size(); b++) {
                dsu.join(e[a], e[b]);
              }
            }
          }
          if (dsu.tamanho[dsu.root(0)] != 4) continue;
          counterHyper[Isomorphism::canonization(motif)]++;
        }
      }
    }
  }
  return counterHyper;
}


// Highly EXPERIMENTAL code


// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf


std::map< std::vector<graph>, long long> ESU::k3Modified(Hypergraph& inputGraph) {
  counterHyper.clear();
  visited.clear();
  for (auto edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    if (motif.is_two_connected()) {
      Hypergraph simpleMotif = motif.filterEdge(2);
      --counterHyper[Isomorphism::canonization(simpleMotif)];
    }
    assert(is_sorted(edge.begin(), edge.end()));
    counterHyper[Isomorphism::canonization(motif)]++;
    assert(motif.getEdgeMaxDeg() == 3);
  }
  h = inputGraph.filterEdge(2);
  Search = HYPERGRAPH;
  vector< vector<int> > g = h.getGraph();
  int n = (int) g.size();
  long long res = 0;
  for (int i = 0; i < n; i++) {
    const int sz = g[i].size();
    res += 1LL * sz * (sz - 1) / 2 ;
  }
  
  
  int triangle = 0;
  int line = 0;
  
  int count = 0;
  
  vector<int> deg(n);
  for (int i = 0; i < n; i++) {
    deg[i] = (int) g[i].size();
  }

  vector< tuple<int, int, int> > tt;
  
  
  for (int i = 0; i < n; i++) {
    vector<int> bigger;
    for (auto& to : g[i]) {
      if (make_pair(deg[to], to) > make_pair(deg[i], i)) bigger.emplace_back(to);
    }
    for (int j = 0; j < (int) bigger.size(); j++) {
      for (int z = j + 1; z < (int) bigger.size(); z++) {
        if (binary_search(g[bigger[j]].begin(), g[bigger[j]].end(), bigger[z])) ++triangle;
      }
    }
  }
  vector< vector<int> > g_line = { {0, 1}, {1, 2} };
  vector< vector<int> > g_trig = { {0, 1}, {1, 2}, {2, 0} };
  
  Hypergraph h_line; 
  h_line.setIncidenceMatrix(g_line); 
  h_line.setN(3);
  
  Hypergraph h_trig; 
  h_trig.setIncidenceMatrix(g_trig);
  h_trig.setN(3);

  counterHyper[Isomorphism::canonization(h_trig)] += triangle;
  counterHyper[Isomorphism::canonization(h_line)] += res - 3 * triangle;
  
  // Remove keys with value 0 (only useful for display only)
  for(auto it = counterHyper.begin(); it != counterHyper.end(); ) {
    if(it -> second == 0) it = counterHyper.erase(it);
    else ++it;
  }
  
  return counterHyper;
}
