#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <set>

#include "Hypergraph.hpp"

using namespace std;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

int gen(int lo, int hi) { return uniform_int_distribution<int>(lo,hi)(rng); }

/*
 * Simple Hypergraph Class
 * No empty hyperedge
 * No duplicate hyperedge
 */

Hypergraph::Hypergraph() {
  //readIncidenceMatrix();
  randomHypergraph();
}

void Hypergraph::randomHypergraph() {
  N = gen(1, 8);
  M = min( (1 << N) - 1, gen(1, 5));
  vector<int> subset;
  for (int mask = 1; mask < (1 << N); mask++) {
    subset.emplace_back(mask);
  }
  shuffle(subset.begin(), subset.end(), rng);
  
  
  vector<int> perm(N);
  iota(perm.begin(), perm.end(), 0);
  shuffle(perm.begin(), perm.end(), rng);
  
  for (int i = 0; i < M; i++) {
    vector<int> edge;
    for (int j = 0; j < N; j++) {
      if ((subset[i] >> j) & 1) {
        edge.emplace_back(perm[j]);
      }
    }
    incidenceMatrix.emplace_back(edge);
  }
  sortAndCheck(incidenceMatrix);
  //vector<int> pickEdge(
}

void Hypergraph::readIncidenceMatrix() {
  /*
   * Format:
   * n # number of hypernodes
   * m # number of hyperedges
   * k_1 a1, a2, ..., ak_1
   * k_2 b1, b2, ..., bk_2
   * ...
   */
   cin >> N >> M;
   incidenceMatrix.clear();
   for (int i = 0; i < M; i++) {
     int k;
     cin >> k;
     vector<int> edge(k);
     for (int j = 0; j < k; j++) {
       int node;
       cin >> node; // number from 0 to n - 1
       edge[j] = node - 1;
       // Check that each node is numbered from 0 to n - 1
       assert(node >= 1 && node <= N);
     }
     incidenceMatrix.emplace_back(edge);
   }
   sortAndCheck(incidenceMatrix);
}
void Hypergraph::sortAndCheck(vector< vector<int> >& edge) const {
  for (int i = 0; i < M; i++) sort(edge[i].begin(), edge[i].end());
  sort(edge.begin(), edge.end());
  edge.erase(unique(edge.begin(), edge.end()), edge.end());
  // Each edge should be unique => List size without duplicates MUST be M
  assert ( (int) edge.size() == M );
}

vector< vector<int> > Hypergraph::applyFunction(const vector<int>& permutation) const {
  vector< vector<int> > modifiedEdgeList;
  for (int i = 0; i < M; i++) {
    vector<int> edge;
    for (auto& node : incidenceMatrix[i]) edge.emplace_back(permutation[node]);
    modifiedEdgeList.emplace_back(edge);
  }
  sortAndCheck(modifiedEdgeList);
  return modifiedEdgeList;
}

void Hypergraph::buildEdgeGraph() { 
  edgeGraph.resize(getEdgeCount());
  for (int i = 0; i < getEdgeCount(); i++) {
    for (int j = 0; j < i; j++) {
      
      int p1 = 0;
      int p2 = 0;
      
      while (p1 < (int) incidenceMatrix[i].size() && p2 < (int) incidenceMatrix[j].size()) {
        
        if (incidenceMatrix[i][p1] == incidenceMatrix[i][p2]) {
          edgeGraph[i].emplace_back(j);
          edgeGraph[j].emplace_back(i);
          //cout << i << ' ' << j << '\n';
          break;
        } else if (incidenceMatrix[i][p1] < incidenceMatrix[i][p2]) {
          ++p1;
        } else {
          ++p2;
        }
      }
      
    }
  }
}


void Hypergraph::enumerateSubgraphs(vector<int> vSubgraph, vector<int> vExtension, int k, int v) {
  if ( (int) vSubgraph.size() == k ) {
     //print subgraph
    for (auto& node : vSubgraph) cout << node + 1 << ' ';
    cout << '\n';
    cout << '\n';
    return;
  }
  
  set<int> chk;
  for (auto to : vSubgraph) {
    chk.insert(to);
    for (auto& nei : edgeGraph[to]) {
      chk.insert(nei);
    }
  }
  
  while (!vExtension.empty()) {
    int w = vExtension.back();
    vExtension.pop_back();
    
    vector<int> vExtensionX = vExtension;
    for (auto& nei : edgeGraph[w]) if (nei > v) {
      if (chk.find(nei) != chk.end()) continue;
      
      
      vExtensionX.emplace_back(nei);
    }
      
    sort(vExtensionX.begin(), vExtensionX.end());
    vExtensionX.erase(unique(vExtensionX.begin(), vExtensionX.end()), vExtensionX.end());
    
    vSubgraph.emplace_back(w);
    enumerateSubgraphs(vSubgraph, vExtensionX, k, v);
    vSubgraph.pop_back();
    
    
  }
};


void Hypergraph::ESU(int k) {
  buildEdgeGraph();
  for (int i = 0; i < getEdgeCount(); i++) {
    vector<int> vExtension;
    for (auto& nei : edgeGraph[i]) if (nei > i) vExtension.emplace_back(nei);
    vector<int> vSubgraph = {i};
    enumerateSubgraphs(vSubgraph, vExtension, k, i);
  }
  
}

bool Hypergraph::isEqual(const vector< vector<int> >& edgeList1) const {
  return incidenceMatrix == edgeList1;
}

vector< vector<int> > Hypergraph::getIncidenceMatrix() const {
  return incidenceMatrix;
}

int Hypergraph::getNodeCount() const {
  return N;
}

int Hypergraph::getEdgeCount() const {
  return M;
}

void Hypergraph::printIncidenceMatrix() const {
  for (int i = 0; i < M; i++) {
    cout << "Hyperedge " << i << ": " << '\n';
    for (auto& node : incidenceMatrix[i]) cout << node << ' ';
    cout << '\n';
  }
}

//void Hypergraph::printEdgeGraph() const {
  //for (int i = 0; i < M; i++) {
    //cout << "Hyperedge " << i << ": " << '\n';
    //for (auto& node : incidenceMatrix[i]) cout << node << ' ';
    //cout << '\n';
  //}
//}

