#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>

#include "Hypergraph.hpp"

using namespace std;

Hypergraph::Hypergraph() {
  readIncidenceMatrix();
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

