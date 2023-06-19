#include <bits/stdc++.h>

#include "Settings.hpp"
#include "nauty.h"
#include "Hypergraph.hpp"
#include "IsomorphismHyper.hpp"


std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

int gen(int lo, int hi) { return std::uniform_int_distribution<int>(lo,hi)(rng); }

/*
 * Simple Hypergraph Class
 * |Hyperedge| >= 2
 * No duplicate hyperedge
 * Every vertex in an edge should be unique!!
 */

Hypergraph::Hypergraph() {
  K = 0; // initially the max. degree is 0
  edgeBySize.resize(MAX_EDGE_SIZE + 1);
  assert(MAX_HYPER_MOTIF_SIZE <= MAX_EDGE_SIZE);
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
   //hashEdge.reserve(M);
   int size = 0;
   incidenceMatrix.clear();
   for (int i = 0; i < M; i++) {
     int k;
     in >> k;
     //k = 2;
     K = max(K, k); // update the max degree
     assert(k <= MAX_EDGE_SIZE); // The max degree should not be bigger than 4.
     vector<int> edge(k);
     for (int j = 0; j < k; j++) {
       int node;
       in >> node; // number from 0 to n - 1
       edge[j] = node - 1;
       // Check that each node is numbered from 0 to n - 1
       assert(node >= 1);
       assert(node <= N);
     }
     //int trash;
     //in >> trash;
     incidenceMatrix.emplace_back(edge);
   }
   //cout << "YES" << endl;
   compress();
   sortAndCheck(incidenceMatrix);
}
void Hypergraph::sortAndCheck(vector< vector<int> >& edge) {
  hashEdge.clear();
  for (int i = 0; i < M; i++) {
    sort(edge[i].begin(), edge[i].end());
    if ((int) edge[i].size() <= MAX_EDGE_SIZE) hashEdge.insert(edge[i]); // insert into hash table
  }
  sort(edge.begin(), edge.end());
  edge.erase(unique(edge.begin(), edge.end()), edge.end());
  // Each edge should be unique => List size without duplicates MUST be M
  assert ( (int) edge.size() == M );
  //if (edge.size() != M) {
    //cout << "F" << '\n';
    //exit(0);
  //}
  //edgeBySize.clear();
  for (int i = 0; i < MAX_EDGE_SIZE; i++) {
    edgeBySize[i].clear();
  } 
  for (int i = 0; i < M; i++) {
    edgeBySize[ (int) edge[i].size() ].emplace_back(i); 
  }
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
  out << "-----------------------------------------------" << '\n';
  out << "Nodes: " << getNodeCount() << '\n';
  out << "Hyperedges:" << getEdgeCount() << '\n';
  for (int i = 0; i < M; i++) {
    out << incidenceMatrix[i].size() << ' ';
    for (auto& node : incidenceMatrix[i]) out << node + 1 << ' ';
    out << '\n';
  }
  out << "-----------------------------------------------" << '\n';
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




Hypergraph Hypergraph::induceSubgraphNoComp(const vector<int>& subgraph, const std::vector<int>& subgraph_compressed)  {
  const int nodes = (int) subgraph.size();
  assert(nodes <= 4); // only for motifs of size 3 and 4!
  Hypergraph h;
  h.setN(nodes);
  vector< vector<int> > adj;
  for (int mask = 0; mask < (1 << nodes); mask++) {
    vector<int> edge;
    vector<int> edgeComp;
    for (int i = 0; i < nodes; i++) {
      if ((mask >> i) & 1) {
        edge.emplace_back(subgraph[i]); // convert to 0-indexed
        edgeComp.emplace_back(subgraph_compressed[i]); // convert to 0-indexed
      }
    }
    sort(edge.begin(), edge.end());
    if (hashEdge.find(edge) != hashEdge.end()) {
      adj.emplace_back(edgeComp);
    }
  }
  h.setIncidenceMatrix(adj);
  h.setN(nodes);
  for (auto edge : adj) { // REMOVEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ME
    for (auto n : edge) assert(n >= 0 && n <= 3);
  }
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

/*
 * Get Degree Seq.
 * For each node we store (e_2, e_3, e_n) ... edges of size n
 */ 
 
vector< vector<int> > Hypergraph::getDegreeSequence() {
  vector< vector<int> > deg;
  for (int i = 0; i < N; i++) {
    deg.emplace_back(vector<int>());
    for (int j = 2; j <= K; j++) {
      int counter = 0;
      for (auto& e : incidenceMatrix) {
        if ( (int) e.size() != j ) continue;
        counter += find(e.begin(), e.end(), i) != e.end();
      }
      deg.back().emplace_back(counter);
    }
  }
  return deg;
}

/*
 * Creates a "random similar" hypergraph"
 */

bool Hypergraph::shuffleEdges(int e1, int e2) {
  assert( incidenceMatrix[e1].size() == incidenceMatrix[e2].size() );
  vector<int> nodes;
  for (auto node : incidenceMatrix[e1]) nodes.emplace_back(node); 
  for (auto node : incidenceMatrix[e2]) nodes.emplace_back(node); 
  shuffle(nodes.begin(), nodes.end(), rng);
  vector<int> n1;
  for (int i = 0; i < (int) nodes.size() / 2; i++) {
    n1.emplace_back(nodes[i]);
  }
  sort(n1.begin(), n1.end());
  n1.erase(unique(n1.begin(), n1.end()), n1.end()); // if there are duplicates, skip
  if ( (int) n1.size() != (int) incidenceMatrix[e1].size() ) return false;
  vector<int> n2;
  for (int i = (int) nodes.size() / 2; i < (int) nodes.size(); i++) {
    n2.emplace_back(nodes[i]);
  }
  sort(n2.begin(), n2.end());
  n2.erase(unique(n2.begin(), n2.end()), n2.end()); // if there are duplicates, skip
  if ( (int) n2.size() != (int) incidenceMatrix[e2].size() ) return false;
  // no overlaping edges ...
  if ( hashEdge.find(n1) != hashEdge.end() ) return false;
  if ( hashEdge.find(n2) != hashEdge.end() ) return false;
  assert( (int) n1.size() == (int) n2.size() );
  hashEdge.erase(incidenceMatrix[e1]);
  hashEdge.erase(incidenceMatrix[e2]);
  incidenceMatrix[e1] = n1;
  incidenceMatrix[e2] = n2;
  hashEdge.insert(incidenceMatrix[e1]);
  hashEdge.insert(incidenceMatrix[e2]);
  return true;
}


// careful - in very small graphs some shuffles might not be possible!!!
void Hypergraph::shuffleHypergraph (int iterations) {
  while (iterations > 0) {
    int edgeSize = gen(2, K); // select a random edge size
    if ( (int) edgeBySize[edgeSize].size() <= 1) continue; // we must have at least two edges two be able to shuffle
    int e1 = gen(0, (int) edgeBySize[edgeSize].size() - 1);
    int e2 = gen(0, (int) edgeBySize[edgeSize].size() - 1);
    if (edgeBySize[edgeSize][e1] == edgeBySize[edgeSize][e2]) continue; // we must shuffle two distinct edges ...
    if (!shuffleEdges(edgeBySize[edgeSize][e1], edgeBySize[edgeSize][e2])) continue;
    --iterations;
  }
  sortAndCheck(incidenceMatrix); // convert graph to "standard" form
}
