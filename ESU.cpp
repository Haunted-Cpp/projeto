#include <bits/stdc++.h>



#include "Hypergraph.hpp"
#include "nauty.h"
#include "ESU.hpp"
#include "IsomorphismHyper.hpp"
#include "Settings.hpp"

#include "FaSE/Fase.h"
#include "FaSE/DynamicGraph.h"
#include "FaSE/GraphMatrix.h"
#include "FaSE/GraphUtils.h"

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

 //Used in Hypergraph
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
      counterHyper[IsomorphismHyper::canonization(motif)]++;
    } else if (Search == CLASS_ONLY) {
      string mat = IsomorphismHyper::canonStr(edgeList, K);
      ++counter[mat]; // we could use a trie like structure here
    } else {
      // Just add the subgraph created
      //subgraphs.emplace_back(edgeList);
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
        if (f[to] == 0) next_extension.emplace_back(to);
        f[to]++;
      }
    }
    subgraph.emplace_back(current_node);
    subgraph_compressed.emplace_back(pos[current_node]);
    enumerateSubgraphs(next_extension);
    subgraph.pop_back();
    subgraph_compressed.pop_back();
    for (auto& to : g[current_node]) if (to > V) --f[to];
  }
};


void ESU::setupAndRun(const vector< vector<int> >& inputGraph, int k) {
  K = k;
  g = inputGraph;
  vector<int> extension;
  for (int i = 0; i < inputGraph.size(); i++) {
    V = i;
    assert(subgraph.empty());
    assert(subgraph_compressed.empty());
    assert(edgeList.empty());
    
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

/*
 * This two methods are used as "brute-force", correctness checker methods
 * Both are slow and SHOULD NOT be used in practice
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
  clearDataStruct();
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
        counterHyper[IsomorphismHyper::canonization(motif)]++;
      }
    }
  }
  return counterHyper;
}

 //This method is used just for debug only!
std::map< std::vector<graph>, long long> ESU::bruteForce4(Hypergraph& inputGraph) {
  clearDataStruct();
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
          counterHyper[IsomorphismHyper::canonization(motif)]++;
        }
      }
    }
  }
  return counterHyper;
}

/*
 * FaSE API to count subgraphs
 */
 
map<string, long long> ESU::FaSE(const vector<pair<int, int> > edges, int k) {
  // Set up FaSE usage
  Graph *G = new DynamicGraph(); // Assuming large scale ...
  bool zeroBased = false;
  GraphUtils::readFile(G, edges, false, false, zeroBased);
  G->sortNeighbours();
  G->makeArrayNeighbours();
  Random::init(time(NULL));
  Fase* fase = new Fase(G, false);
  // Network-census
  fase->runCensus(k);
  // Prepare output format
  map<string, long long> output;
  for (auto& element : fase->subgraphCount()) {
    output[element.second] = element.first;
  }
  // Delete allocated memory
  delete fase;
  delete G;
  
  return output;
}

Hypergraph ESU::binaryToHyper(string str, int k) {
  vector< vector<int> > adj;
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < i; j++) {
      if (str[i * k + j] == '1') {
        adj.push_back({min(i, j), max(i, j)});
      }
    }
  }
  Hypergraph h;
  h.setIncidenceMatrix(adj);
  h.setN(k);
  return h;
}

/*
 * Baseline method for K_3
 */

// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

// implementar o caso K = 3 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

std::map< std::vector<graph>, long long> ESU::k3(Hypergraph& inputGraph) {
  clearDataStruct();
  for (auto edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    visited.insert(edge); // it will insert the 3 nodes just visited
    counterHyper[IsomorphismHyper::canonization(motif)]++;
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  setupAndRun(h.getGraph(), 3);
  return counterHyper;
}


/*
 * The following are a series of proposed methods for K_3
 */

// Auxiliary method to compute intermediate form

void ESU::k3IntermediateForm(Hypergraph& inputGraph) {
  for (auto edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    if (motif.is_two_connected()) {
      Hypergraph simpleMotif = motif.filterEdge(2);
      --counterHyper[IsomorphismHyper::canonization(simpleMotif)];
      // this occurence will be found by ESU!
    }
    assert(is_sorted(edge.begin(), edge.end()));
    counterHyper[IsomorphismHyper::canonization(motif)]++;
    assert(motif.getEdgeMaxDeg() == 3);
  }
  h = inputGraph.filterEdge(2);
}


// TRIANGLE method

std::map< std::vector<graph>, long long> ESU::k3Triangle(Hypergraph& inputGraph) {
  clearDataStruct();
  k3IntermediateForm(inputGraph);
  vector< vector<int> > g = h.getGraph();
  int n = (int) g.size();
  long long res = 0;
  for (int i = 0; i < n; i++) {
    const int sz = g[i].size();
    res += 1LL * sz * (sz - 1) / 2 ;
  }
  long long triangle = 0;
  vector<int> deg(n);
  for (int i = 0; i < n; i++) {
    deg[i] = (int) g[i].size();
  }
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

  counterHyper[IsomorphismHyper::canonization(h_trig)] += triangle;
  counterHyper[IsomorphismHyper::canonization(h_line)] += res - 3 * triangle;
  
  return counterHyper;
}

// FASE method

std::map< std::vector<graph>, long long> ESU::k3Fase(Hypergraph& inputGraph) {
  clearDataStruct();
  k3IntermediateForm(inputGraph);
  vector< pair<int, int> > edges;
  for (auto edge : h.getIncidenceMatrix()) {
    if ((int) edge.size() == 2) edges.emplace_back(edge[0], edge[1]); 
  }
  auto census = FaSE(edges, 3);
  for (auto& [str, cnt] : census) {
    Hypergraph hyper = binaryToHyper(str, 3);
    counterHyper[IsomorphismHyper::canonization(hyper)] += cnt;
  }
  return counterHyper;
}

// SIMPLE ESU

std::map< std::vector<graph>, long long> ESU::k3ESU(Hypergraph& inputGraph) {
  clearDataStruct();
  k3IntermediateForm(inputGraph);
  h = inputGraph;
  Search = CLASS_ONLY;
  setupAndRun(h.getGraph(), 3);
  return counterHyper;
}


/*
 * Baseline method for K_4
 */

// implementar o caso K = 4 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

std::map< std::vector<graph>, long long> ESU::k4(Hypergraph& inputGraph) {
  clearDataStruct();
  for (auto& edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 4) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    visited.insert(edge); // it will insert the 3 nodes just visited
    counterHyper[IsomorphismHyper::canonization(motif)]++;
    assert(motif.getEdgeMaxDeg() == 4);
    
    //cout << "YES" << '\n';
    //exit(0);
  }
  Hypergraph reducedGraph = inputGraph.filterEdge(3); // at most 3 edges
  assert( reducedGraph.getEdgeCount() == reducedGraph.getIncidenceMatrix().size() );
  vector< vector<int> > g(reducedGraph.getNodeCount());
  int edge = 0;
  for (auto e : reducedGraph.getIncidenceMatrix()) {
    for (auto e_i : e) {
      g[e_i].emplace_back(edge);
    }
    ++edge;
  }
  edge = -1;
  for (auto e : reducedGraph.getIncidenceMatrix()) {
    ++edge;
    if ( (int) e.size() != 3 ) continue;
    for (auto node : e) {
      for (auto e_i : g[node]) {
        vector<int> e1 = reducedGraph.getEdge(e_i);
        for (auto add : reducedGraph.getEdge(edge)) {
          e1.emplace_back(add);
        }
        sort(e1.begin(), e1.end());
        e1.erase(unique(e1.begin(), e1.end()), e1.end());
        if (e1.size() == 4 && visited.find(e1) == visited.end()) {
          visited.insert(e1);
          Hypergraph motif = inputGraph.induceSubgraph(e1);
          assert(is_sorted(e1.begin(), e1.end()));
          counterHyper[IsomorphismHyper::canonization(motif)]++;
        }
      }
    }
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  setupAndRun(h.getGraph(), 4);
  return counterHyper;
}


 
/*
 * The following are a series of proposed methods for K_4
 */


// Auxiliary method to compute intermediate form

void ESU::k4IntermediateForm(Hypergraph& inputGraph) {
  for (auto& edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 4) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    assert(is_sorted(edge.begin(), edge.end()));
    if (motif.is_two_connected()) {
      Hypergraph simpleMotif = motif.filterEdge(2);
      --counterHyper[IsomorphismHyper::canonization(simpleMotif)];
      // this occurence will be found by ESU!
    }
    counterHyper[IsomorphismHyper::canonization(motif)]++;
    assert(motif.getEdgeMaxDeg() == 4);
  }
  Hypergraph reducedGraph = inputGraph.filterEdge(3); // at most 3 edges
  vector< vector<int> > g(reducedGraph.getNodeCount());
  int edge = 0;
  for (auto e : reducedGraph.getIncidenceMatrix()) {
    for (auto e_i : e) g[e_i].emplace_back(edge);
    ++edge;
  }
  edge = -1;
  vector<int> f(reducedGraph.getNodeCount());
  std::stack<int> rem;
  for (auto e : reducedGraph.getIncidenceMatrix()) {
    ++edge;
    if ( (int) e.size() != 3 ) continue;
    for (auto node : e) {
      for (auto e_i : g[node]) {
        vector<int> nodes = reducedGraph.getEdge(e_i);
        for (auto add : reducedGraph.getEdge(edge)) nodes.emplace_back(add);
        sort(nodes.begin(), nodes.end());
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
        if (nodes.size() != 4) continue;
        int added_node = -1;
        for (auto c : nodes) {
          if (find(e.begin(), e.end(), c) == e.end()) {
            added_node = c;
            break;
          }
        }
        assert(added_node != -1);
        if (f[added_node]) continue;
        f[added_node] = 1;
        rem.push(added_node);
        assert(is_sorted(nodes.begin(), nodes.end()));
        Hypergraph motif1 = inputGraph.induceSubgraphSkipComp(nodes);
        if (motif1.getEdgeMaxDeg() != 3) continue;
        int skip = 0;
        for (auto motif_edge : motif1.getIncidenceMatrix()) {
          if (motif_edge.size() != 3) continue;
          if (motif_edge <= e) continue;
          int missing = -1;
          for (auto n : nodes) {
            if (find(motif_edge.begin(), motif_edge.end(), n) == motif_edge.end()) {
              missing = n;
              break;
            }
          }
          assert(missing != -1);
          for (int mask1 = 1; mask1 < (1 << 3) - 1; mask1++) { // all subsets of current edge
            vector<int> new_edge1;
            for (int x = 0; x < 3; x++) {
              if ((mask1 >> x) & 1) {
                new_edge1.emplace_back(motif_edge[x]);
              }
            }
            new_edge1.emplace_back(missing);
            sort(new_edge1.begin(), new_edge1.end());
            if (! reducedGraph.validEdge(new_edge1) ) continue;
            skip = 1;
            break;
          }
        }
        if (skip) continue; 
        Hypergraph motif = inputGraph.induceSubgraph(nodes);
        assert(motif.getEdgeMaxDeg() == 3);
        if (motif.is_two_connected()) {
          Hypergraph simpleMotif = motif.filterEdge(2);
          --counterHyper[IsomorphismHyper::canonization(simpleMotif)];
          // this occurence will be found by ESU!
        }
        assert(is_sorted(nodes.begin(), nodes.end()));
        counterHyper[IsomorphismHyper::canonization(motif)]++;
      }
    }
    while (!rem.empty()) {
      f[rem.top()] = 0;
      rem.pop();
    }
  }
  
  h = inputGraph.filterEdge(2);
}

// FaSE
std::map< std::vector<graph>, long long> ESU::k4Fase(Hypergraph& inputGraph) {
  clearDataStruct();
  k4IntermediateForm(inputGraph);
  vector< pair<int, int> > edges;
  for (auto edge : h.getIncidenceMatrix()) {
    if ((int) edge.size() == 2) edges.emplace_back(edge[0], edge[1]); 
  }
  auto census = FaSE(edges, 4);
  for (auto& [str, cnt] : census) {
    Hypergraph hyper = binaryToHyper(str, 4);
    counterHyper[IsomorphismHyper::canonization(hyper)] += cnt;
  }
  return counterHyper;
  
  
}

// SIMPLE ESU

std::map< std::vector<graph>, long long> ESU::k4ESU(Hypergraph& inputGraph) {
  clearDataStruct();
  k4IntermediateForm(inputGraph);
  h = inputGraph;
  Search = CLASS_ONLY;
  setupAndRun(h.getGraph(), 4);
  return counterHyper;
}


/*
 * 
 * 
 */


void ESU::printResults(std::chrono::time_point<std::chrono::steady_clock> startTime, std::chrono::time_point<std::chrono::steady_clock> endTime, map< vector<graph>, long long> subgraph_count, int k, bool detailedOutput, ostream& out) {
  // Remove keys with value 0 (only useful for display only)
  for(auto it = counterHyper.begin(); it != counterHyper.end(); ) {
    if(it -> second == 0) it = counterHyper.erase(it);
    else ++it;
  }
  out << "-----------------------------------------------" << endl;
  out << "Network census completed in: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
  long long total_subgraph = 0;
  for (auto& [a, b] : subgraph_count) {
    total_subgraph += b;
  }
  out << total_subgraph << " hyper-subgraphs extracted" << endl;
  out << subgraph_count.size() << " types of hyper-subgraphs found" << endl;
  out << "-----------------------------------------------" << endl;
  int counter = 0;
  for (auto& [a, b] : subgraph_count) {
    out << "Type #" << ++counter << ": " << b << endl;
  }
  out << "-----------------------------------------------" << endl;
  if (detailedOutput) {
    counter = 0;
    for (auto& [a, b] : subgraph_count) {
      auto adj = IsomorphismHyper::getHypergraph(a);
      Hypergraph h;
      h.setIncidenceMatrix(adj);
      h.setN(k);
      out << "Hyper-subgraph #" << ++counter << endl;
      h.printIncidenceMatrix();
      out << "Number of occurences: " << b << endl;
      out << "-----------------------------------------------" << endl;
    }
  }
  out << flush;
}


void ESU::clearDataStruct() {
  edgeList.clear();
  subgraph.clear();
  counterHyper.clear();
  visited.clear();
  subgraphs.clear();
  counter.clear();
}


void ESU::networkCensus(Hypergraph& h, int motifSize, bool detailedOutput, ostream& out) {
  clearDataStruct();
  int k;
  auto startTime = steady_clock::now();
  if (motifSize == 3) { // Execute our fastest method K = 3, TRIANGLE
    k = 3;
    //k3Triangle(h);
    //k3Fase(h);
    k3(h);
  } else { // Execute our fastest method K=4, FASE
    k = 4;
    k4(h);
    //k4ESU(h);
    //k4Fase(h);
  }
  auto endTime = steady_clock::now();
  printResults(startTime, endTime, counterHyper, k, detailedOutput, out);
}
  
  
// ADD Alternative method  
  
void ESU::findMotifs(Hypergraph& h, int motifSize, bool detailedOutput, ostream& out) {
  
  auto census = (motifSize == 3 ? ESU::k3Triangle(h) : ESU::k4Fase(h));
  map< vector<graph>, vector<int> > sample;
  const int NUMBER_NETWORKS = 100;
  for (int i = 0; i < NUMBER_NETWORKS; i++) {
    h.shuffleHypergraph(100);
    auto count = (motifSize == 3 ? ESU::k3Triangle(h) : ESU::k4Fase(h));
    for (auto [a, b] : census) {
      sample[a].emplace_back(count[a]);
    }
    for (auto [a, b] : count) {
      out << b << endl;
    }
    out << "----" << endl;
  }
  double sum = 0;
  vector<double> sp;
  for (auto [a, b] : census) {
    double mean = 0;
    for (auto value : sample[a]) {
      mean += value;
    }
    mean /= NUMBER_NETWORKS;
    sp.emplace_back( (b - mean) / (b + mean + 4) );
    sum += sp.back() * sp.back();
  }
  sum = sqrt(sum);
  int counter = 0;
  for (auto [a, b] : census) {
    out << "Hyper-subgraph #" << ++counter << endl;
    if (detailedOutput) {
      auto adj = IsomorphismHyper::getHypergraph(a);
      Hypergraph h;
      h.setIncidenceMatrix(adj);
      h.setN(motifSize);
      h.printIncidenceMatrix();
    }
    out << "Number of occurences: " << b << endl;
    out << "Significance Profile: " << sp[counter - 1] / sum << endl;
    out << "-----------------------------------------------" << endl;
  }
}

