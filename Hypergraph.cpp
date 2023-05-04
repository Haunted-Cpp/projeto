#include <bits/stdc++.h>

#include "Settings.hpp"
#include "nauty.h"
#include "Hypergraph.hpp"
#include "Isomorphism.hpp"


using namespace std;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

int gen(int lo, int hi) { return uniform_int_distribution<int>(lo,hi)(rng); }

/*
 * Simple Hypergraph Class
 * No empty hyperedge
 * No duplicate hyperedge
 * Every edge vertex should be unique!!
 * No sefl-loops
 */

Hypergraph::Hypergraph() {
  //readIncidenceMatrix();
  //randomHypergraph();
  // Don't call here - Otherwise when a small, empty Hypergraph is created this methods will be called
}

void Hypergraph::randomHypergraph(int n, int m, int maxDegree) {
  N = n;
  M = m;
  
  vector<int> subset;
  for (int mask = 1; mask < (1 << N); mask++) {
    if (__builtin_popcount(mask) > maxDegree) continue;
    if (__builtin_popcount(mask) < 2) continue;
    subset.emplace_back(mask);
  }
  assert( (int) subset.size() >= M );
  shuffle(subset.begin(), subset.end(), rng);
  
  
  vector<int> perm(N);
  iota(perm.begin(), perm.end(), 0);
  shuffle(perm.begin(), perm.end(), rng);
  
  incidenceMatrix.clear();
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
}

void Hypergraph::readIncidenceMatrix(istream& in) {
  /*
   * Format:
   * n # number of hypernodes
   * m # number of hyperedges
   * k_1 a1, a2, ..., ak_1
   * k_2 b1, b2, ..., bk_2
   * ...
   */
   in >> N >> M;
   assert(N <= MAX_INPUT_N);
   incidenceMatrix.clear();
   for (int i = 0; i < M; i++) {
     int k;
     in >> k;
     vector<int> edge(k);
     for (int j = 0; j < k; j++) {
       int node;
       in >> node; // number from 0 to n - 1
       edge[j] = node - 1;
       // Check that each node is numbered from 0 to n - 1
       assert(node >= 1 && node <= N);
     }
     incidenceMatrix.emplace_back(edge);
   }
   sortAndCheck(incidenceMatrix);
}
void Hypergraph::sortAndCheck(vector< vector<int> >& edge) {
  hashEdge.clear();
  for (int i = 0; i < M; i++) {
    sort(edge[i].begin(), edge[i].end());
    if ((int) edge[i].size() <= MAX_HYPER_MOTIF_SIZE) hashEdge.insert(edge[i]); // insert into hash table
  }
  sort(edge.begin(), edge.end());
  edge.erase(unique(edge.begin(), edge.end()), edge.end());
  // Each edge should be unique => List size without duplicates MUST be M
  assert ( (int) edge.size() == M );
}

vector< vector<int> > Hypergraph::applyFunction(const vector<int>& permutation) {
  vector< vector<int> > modifiedEdgeList;
  for (int i = 0; i < M; i++) {
    vector<int> edge;
    for (auto& node : incidenceMatrix[i]) edge.emplace_back(permutation[node]);
    modifiedEdgeList.emplace_back(edge);
  }
  sortAndCheck(modifiedEdgeList);
  return modifiedEdgeList;
}

vector< vector<int> > Hypergraph::buildEdgeGraph() { 
  vector< vector<int> > edgeGraph(getEdgeCount());
  for (int i = 0; i < getEdgeCount(); i++) {
    for (int j = 0; j < i; j++) {
      int p1 = 0;
      int p2 = 0;
      while (p1 < (int) incidenceMatrix[i].size() && p2 < (int) incidenceMatrix[j].size()) {
        if (incidenceMatrix[i][p1] == incidenceMatrix[j][p2]) {
          edgeGraph[i].emplace_back(j);
          edgeGraph[j].emplace_back(i);
          break;
        } else if (incidenceMatrix[i][p1] < incidenceMatrix[j][p2]) {
          ++p1;
        } else {
          ++p2;
        }
      }
    }
  }
  return edgeGraph;
}

vector< vector<int> > Hypergraph::buildVertexGraph(int k) {
  vector< vector<int> > nei( getNodeCount() );
  for (auto& edge : incidenceMatrix) { // O (m * k ^ 2), but since k <= 10 ... ok
    if ( (int) edge.size() <= k ) {
      for (int i = 0; i < (int) edge.size(); i++) {
        for (int j = i + 1; j < (int) edge.size(); j++) {
          nei[ edge[i] ].emplace_back( edge[j] );
          nei[ edge[j] ].emplace_back( edge[i] );
        }
      }
    }
  }
  for (auto& node : nei) { // remove duplicates 
    sort(node.begin(), node.end());
    node.erase(unique(node.begin(), node.end()), node.end());
  }
  return nei;  
}



bool Hypergraph::isEqual(const vector< vector<int> >& edgeList1)  {
  return incidenceMatrix == edgeList1;
}

vector< vector<int> > Hypergraph::getIncidenceMatrix()  {
  return incidenceMatrix;
}

int Hypergraph::getNodeCount()  {
  return N;
}

int Hypergraph::getEdgeCount()  {
  return M;
}

vector<int> Hypergraph::getEdge(int n) {
  assert(n >= 0 && n < M);
  return incidenceMatrix[n];
}

void Hypergraph::printIncidenceMatrix(ostream& out)  {
  out << "---------" << '\n';
  out << getNodeCount() << '\n';
  out << getEdgeCount() << '\n';
  for (int i = 0; i < M; i++) {
    out << incidenceMatrix[i].size() << ' ';
    for (auto& node : incidenceMatrix[i]) out << node + 1 << ' ';
    out << '\n';
  }
  out << "---------" << '\n';
}

/*
 * If called with every possible produced subgraph it will print every non-induced subgraph
 * ------------ Make sure size is correct!!! ----
 * Note that number of nodes is not known!!
 */ 
void Hypergraph::printEdgeSubgraph(vector< pair<int, int> >& edgeList) {
  cout << "Selected subgraph" << '\n';
  vector<int> nodes;
  for (auto& [a, b] : edgeList) {
    cout << a + 1 << ' ' << b + 1 << '\n';
    nodes.emplace_back(a);
    nodes.emplace_back(b);
  }
  sort(nodes.begin(), nodes.end());
  nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
  cout << "Nodes in each edge" << '\n';
  for (auto& node : nodes) {
    for (auto& value : incidenceMatrix[node]) cout << value + 1 << ' ';
    cout << '\n';
  }
}


Hypergraph Hypergraph::filterEdge(int maximumSize) {
  Hypergraph h;
  vector< vector<int> > adj;
  for (auto& edge : incidenceMatrix) if ((int) edge.size() <= maximumSize) adj.emplace_back(edge);
  h.setN(getNodeCount());
  h.setIncidenceMatrix(adj);
  //h.compress(); // from 0 to n - 1 // This bugs K3 -> Notice how we need the vertices with the correct label to induce the subgraph properly!
  return h;
}


// Returns graph with ONLY order 2 links
vector< vector<int> > Hypergraph::getGraph() {
  vector< vector<int> > graph(getNodeCount());
  vector< vector<int> > adj = (this -> filterEdge(2)).getIncidenceMatrix();
  for (auto& edge : adj) {
    if ( (int) edge.size() != 2 ) continue; // we are ignoring self-loops
    graph[ edge[0] ].emplace_back( edge[1] );
    graph[ edge[1] ].emplace_back( edge[0] );
  }
  return graph;
}

bool Hypergraph::is_two_connected() {
  vector< vector<int> > graph = getGraph();
  vector<int> vis(getNodeCount());
  vis[0] = 1; // Notice how nodes should be numbered from 0 to n - 1 here ...
  queue<int> q;
  q.emplace(0);
  while (!q.empty()) {
    int node = q.front();
    q.pop();
    for (auto& to : graph[node]) {
      assert(to < getNodeCount());
      if (!vis[to]) {
        vis[to] = 1;
        q.emplace(to);
      }
    }
  }
  return count(vis.begin(), vis.end(), 1) == getNodeCount();
}

void Hypergraph::setN(int n) {
  N = n;
}
void Hypergraph::setM(int m) {
  M = m;
}

void Hypergraph::setIncidenceMatrix(std::vector< std::vector<int> >& adj) {
  setM( (int) adj.size() );
  incidenceMatrix = adj;
  sortAndCheck(incidenceMatrix);
}

void Hypergraph::compress() {
  vector<int> nodes;
  for (auto& edge : incidenceMatrix) {
    for (auto& vertex : edge) nodes.emplace_back(vertex);
  }
  sort(nodes.begin(), nodes.end());
  nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
  for (auto& edge : incidenceMatrix) {
    for (auto& vertex : edge) {
      vertex = lower_bound(nodes.begin(), nodes.end(), vertex) - nodes.begin();
    }
  }
  sortAndCheck(incidenceMatrix);
}

// 0 indexeddddddddddddddd !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// we assume subgraph is 0-indexed!!!
Hypergraph Hypergraph::induceSubgraph(const vector<int>& subgraph) {
  const int nodes = (int) subgraph.size();
  assert(nodes <= 4); // only for motifs of size 3 and 4!
  Hypergraph h;
  h.setN(nodes);
  vector< vector<int> > adj;
  for (int mask = 0; mask < (1 << nodes); mask++) {
    vector<int> edge;
    for (int i = 0; i < nodes; i++) {
      if ((mask >> i) & 1) edge.emplace_back(subgraph[i]); // convert to 0-indexed
    }
    sort(edge.begin(), edge.end());
    if (hashEdge.find(edge) != hashEdge.end()) {
      adj.emplace_back(edge);
    }
  }
  h.setIncidenceMatrix(adj);
  h.compress();
  return h;
}




Hypergraph Hypergraph::induceSubgraphNoComp(const vector<int>& subgraph) {
  const int nodes = (int) subgraph.size();
  assert(nodes <= 4); // only for motifs of size 3 and 4!
  Hypergraph h;
  h.setN(nodes);
  vector< vector<int> > adj;
  for (int mask = 0; mask < (1 << nodes); mask++) {
    vector<int> edge;
    for (int i = 0; i < nodes; i++) {
      if ((mask >> i) & 1) edge.emplace_back(subgraph[i]); // convert to 0-indexed
    }
    sort(edge.begin(), edge.end());
    if (hashEdge.find(edge) != hashEdge.end()) {
      adj.emplace_back(edge);
    }
  }
  h.setIncidenceMatrix(adj);
  //h.compress();
  return h;
}


void Hypergraph::printIncidenceMatrix() {
  printIncidenceMatrix(cout);
}

void Hypergraph::saveToFile(string filename) {
  ofstream fout;
  fout.open(filename);
  if (fout.fail()) {
    cout << filename << " could not be opened" << '\n';
    exit(EXIT_FAILURE);
  };
  printIncidenceMatrix(fout);
  fout.close();
}

void Hypergraph::readFromFile(string filename) {
  ifstream fin;
  fin.open(filename);
  if (fin.fail()) {
    cout << filename << " could not be opened" << '\n';
    exit(EXIT_FAILURE);
  };
  readIncidenceMatrix(fin);
  fin.close();
}

void Hypergraph::readFromStdin() {
  readIncidenceMatrix(cin);
}

bool Hypergraph::validEdge(std::vector<int> edge) { // it MUST be sorted
  assert(is_sorted(edge.begin(), edge.end()));
  assert(edge.size() <= MAX_HYPER_MOTIF_SIZE) ;
  return hashEdge.find(edge) != hashEdge.end();
}

int Hypergraph::getEdgeMaxDeg() {
  int mx = 0;
  for (auto edge : incidenceMatrix) mx = max(mx, (int) edge.size());
  return mx;
}
