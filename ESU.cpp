#include <bits/stdc++.h>

#include "Hypergraph.hpp"
#include "nauty.h"
#include "ESU.hpp"
#include "Isomorphism.hpp"

using namespace std;

/* 
 * Initialize static variables
 */


int ESU::f[MAX_INPUT_N];
//int ESU::ESU_FIXED_extension[MAX_INPUT_N];
//int ESU::ptr = 0;

int ESU::pos[MAX_INPUT_N];
vector< std::pair<int, int> > ESU::edgeList;

std::vector<int> ESU::subgraph;
std::vector<int> ESU::subgraph_compressed;

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

//int search_lim;
//

// Debug only
map< vector<graph>, int> ESU::counterHyperK4;
map< vector<graph>, int> ESU::counterHyperK3;
map< vector<graph>, int> ESU::counterHyperK3D;
map< vector<graph>, int> ESU::counterHyperBF4;
map< vector<graph>, int> ESU::counterHyperBF3;

//8
//4327
//958
//801
//2652
//1346
//208
//8
//18

int res = 0;

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
  extension.clear();
  neiSubgraph.clear();
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
    //break;
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

//void ESU::setGraph(const vector< vector<int> >& inputGraph) {
  //g = vector< vector<int> > (inputGraph.size());
  //for (int i = 0; i < inputGraph.size(); i++) {
    //for (auto to : inputGraph[i]) if (to >= v) g[i].emplace_back(to);
  //}
//}

//int ESU::getNodeCount() {
  //return (int) g.size();
//}

// implementar o caso K = 3 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

std::map< std::vector<graph>, int> ESU::k3(Hypergraph& inputGraph) {
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
  vector< pair<int, vector<graph> > > xx;
  for (auto [x, cnt] : counterHyper) {
    xx.emplace_back(cnt, x);
  }
  sort(xx.rbegin(), xx.rend());
  for (auto [a, b] : xx) {
    cout << a << '\n';
  }
  cout << "VISITED: " << visited.size() << '\n';
  return counterHyper;
}



//void ESU::k3(Hypergraph& inputGraph) {
  //search_lim = 4;
  //counterHyper.clear();
  //visited.clear();
  //for (auto edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    //if (edge.size() != 3) continue;
    //Hypergraph motif = inputGraph.induceSubgraph(edge);
    //if (motif.is_two_connected()) {
      //continue; // this occurence will be found by ESU!
    //}
    //assert(is_sorted(edge.begin(), edge.end()));
    //counterHyper[Isomorphism::canonization(motif)]++;
    //assert(motif.getEdgeMaxDeg() == 3);
  //}
  //h = inputGraph;
  //Search = HYPERGRAPH;
  //setupAndRun(h.getGraph(), 3);
  //cout << "Counter: " << counterHyper.size() << '\n';
  //vector<int> xx;
  //counterHyperK3 = counterHyper;
  //for (auto [x, cnt] : counterHyper) {
    //cout << cnt << '\n';
    //xx.emplace_back(cnt);
  //}
  //sort(xx.rbegin(), xx.rend());
  //for (auto cnt : xx) cout << cnt << '\n';
  //exit(0);
  //assert(counterHyperK3 == counterHyperBF3);
//}

// implementar o caso K = 4 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

std::map< std::vector<graph>, int> ESU::k4(Hypergraph& inputGraph) {
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
  //cout << "Counter: " << counterHyper.size() << '\n';
  vector<int> xx;
  
  counterHyperK4 = counterHyper;
  for (auto [x, cnt] : counterHyper) {
    //cout << cnt << '\n';
    xx.emplace_back(cnt);
  }
  sort(xx.rbegin(), xx.rend());
  //for (auto cnt : xx) cout << cnt << '\n';
  //assert(counterHyperK4 == counterHyperBF4);
  
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

std::map< std::vector<graph>, int> ESU::bruteForce3(Hypergraph& inputGraph) {
  
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
  
  //cout << "Counter: " << counterHyper.size() << '\n';
  vector<int> xx;
  counterHyperBF3 = counterHyper;
  for (auto [x, cnt] : counterHyper) {
    xx.emplace_back(cnt);
  }
  //cout << "YES" << endl;
  //exit(0);
  return counterHyper;
  //sort(xx.rbegin(), xx.rend());
  //for (auto cnt : xx) cout << cnt << '\n';
}

 //This method is used just for debug only!
std::map< std::vector<graph>, int> ESU::bruteForce4(Hypergraph& inputGraph) {
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
  //cout << "Counter: " << counterHyper.size() << '\n';
  vector<int> xx;
  counterHyperBF4 = counterHyper;
  for (auto [x, cnt] : counterHyper) {
    xx.emplace_back(cnt);
  }
  sort(xx.rbegin(), xx.rend());
  //for (auto cnt : xx) cout << cnt << '\n';
  return counterHyper;
}





// Highly EXPERIMENTAL code


// implementar o caso K = 3 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf









//int count_connected_fast() {
  
  //return res;
//}



//void count() {

//}


//3883887
//360367
//240888
//62507
//10203
//3526
//3




std::map< std::vector<graph>, int> ESU::k3Modified(Hypergraph& inputGraph) {
  counterHyper.clear();
  visited.clear();
  for (auto edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    if (motif.is_two_connected()) {
      Hypergraph simpleMotif = motif.filterEdge(2);
      --counterHyper[Isomorphism::canonization(simpleMotif)];
      //continue; // this occurence will be found by ESU!
    }
    assert(is_sorted(edge.begin(), edge.end()));
    counterHyper[Isomorphism::canonization(motif)]++;
    assert(motif.getEdgeMaxDeg() == 3);
  }
  h = inputGraph.filterEdge(2);
  //h.printIncidenceMatrix()
  //h.printIncidenceMatrix();
  
  Search = HYPERGRAPH;
  //setupAndRun(h.getGraph(), 3);
  
  // We must count the number of things
  
  
  
  vector< vector<int> > g = h.getGraph();
  int n = (int) g.size();
  long long res = 0;
  for (int i = 0; i < n; i++) {
    const int sz = g[i].size();
    res += 1LL * sz * (sz - 1) / 2 ;
  }
  
  
  auto jj = [&](int a, int b) -> bool {
    return find(g[a].begin(), g[a].end(), b) != g[a].end();
  };
  
  int triangle = 0;
  int line = 0;
  
  int count = 0;
  
  vector<int> deg(n);
  for (int i = 0; i < n; i++) {
    deg[i] = (int) g[i].size();
  }
  //int rem = 0;
  vector< tuple<int, int, int> > tt;
  
  
  
  
  for (int i = 0; i < n; i++) {
    vector<int> bigger;
    for (auto& to : g[i]) {
      //bigger.emplace_back(to);
      if ( make_pair(deg[to], to) > make_pair(deg[i], i) ) bigger.emplace_back(to);
      
    }
    for (int j = 0; j < (int) bigger.size(); j++) {
      for (int z = j + 1; z < (int) bigger.size(); z++) {
        //if (bigger[i] == bigger[z]) continue;
        if (jj(bigger[j], bigger[z])) {
          
          ++triangle;
        }
      }
    }
  }
  
  //sort(tt.begin(), tt.end());
  //tt.erase(unique(tt.begin(), tt.end()), tt.end());
  
  //for (auto& [a, b, c] : tt) {
    //cout << a << ' ' << b << ' ' << c << '\n';
    //assert(jj(a, b) && jj(b, c) && jj(c, a));
  //}
  
  //cout << "HERE: " << triangle << ' ' << (int) tt.size() << '\n';
  //triangle = tt.size();
  //cout << "HERE: " << triangle << ' ' << (int) tt.size() << '\n';
  
  //exit(0);
  //exit(0);
  
  //for (int i = 0; i < n; i++) {
    //cout << "S: " << i << ": ";
    //for (auto& nei : g[i]) {
      //cout << nei << ' ';
    //}
    //cout << '\n';
  //}
  
  //cout << "HERE:" << '\n';
  //cout << "T: " << triangle << '\n';
  //cout << "L: " << line << '\n';
  
  //assert( line == res - 3 * triangle);
  
  //exit(0);
  
  vector< vector<int> > g_line = { {1}, {0, 2}, {1} };
  vector< vector<int> > g_trig = { {1, 2}, {0, 2}, {1, 0} };
  
  counterHyper[Isomorphism::canonization(g_trig)] += triangle;
  counterHyper[Isomorphism::canonization(g_line)] += res - 3 * triangle;
  
  
  
  
  //cout << "Counter: " << counterHyper.size() << '\n';
  //vector<int> xx;
  
  // Remove keys with value 0 (only useful for display only)
  for(auto it = counterHyper.begin(); it != counterHyper.end(); ) {
    if(it->second == 0) it = counterHyper.erase(it);
    else ++it;
  }
  //return counterHyper;
  //counterHyperK3D = counterHyper;
  //for (auto [x, cnt] : counterHyper) {
    //cout << cnt << '\n';
    //xx.emplace_back(cnt);
  //}
  //sort(xx.rbegin(), xx.rend());
  //for (auto cnt : xx) cout << cnt << '\n';
  
  return counterHyper;
  
  //exit(0);
  
  //if (counterHyperK3D != counterHyperBF3) cout << "no" << '\n';
  //else cout << "yes" << '\n';
  //cout << counterHyperK3D.size() << ' ' << counterHyperBF3.size() << '\n';
  //assert( !counterHyperK3D.empty() );
  //assert( !counterHyperBF3.empty() );
  //assert(counterHyperK3D == counterHyperBF3);
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
