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

std::vector<int> ESU::subgraph;

int ESU::K; // Subgraph size
int ESU::V; // Initial node
std::chrono::time_point<std::chrono::steady_clock> ESU::startAlgo;

// Used in classical graph
vector< vector<int> > ESU::g;

 //Used in Hypergraph
Hypergraph ESU::h;
map< int, long long> ESU::counterHyper;
std::unordered_set< vector<int>, HashFunction> ESU::visited;

enum option { HYPERGRAPH, GRAPH };
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
      Hypergraph motif = h.induceSubgraphNoComp(subgraph_ordered);
      counterHyper[IsomorphismHyper::getLabel(motif)]++;
    } else {
      // Since this procedure is method is only used in graphs, we can void the induceSubgraph method and use a simpler one
      vector< vector<int> > adj;
      for (int i = 0; i < (int) subgraph.size(); i++) {
        for (int j = i + 1; j < (int) subgraph.size(); j++) {
          if (h.validEdge({ std::min(subgraph[i], subgraph[j]), std::max(subgraph[i], subgraph[j])})) {
            adj.push_back({i, j});
          }
        }
      }
      counterHyper[IsomorphismHyper::getLabel(adj)]++;
    }
    return;
  }
  while (!extension.empty()) {
    int current_node = extension.back();
    extension.pop_back();
    next_extension = extension;
    for (auto& to : g[current_node]) {
      if (to > V) {
        if (f[to] == 0) next_extension.emplace_back(to);
        f[to]++;
      }
    }
    subgraph.emplace_back(current_node);
    enumerateSubgraphs(next_extension);
    subgraph.pop_back();
    for (auto& to : g[current_node]) if (to > V) --f[to];
  }
};


void ESU::setupAndRun(const vector< vector<int> >& inputGraph, int k) {
  K = k;
  g = inputGraph;
  vector<int> extension;
  for (int i = 0; i < inputGraph.size(); i++) {
    V = i;
    extension = {i};
    ++f[i];
    enumerateSubgraphs(extension);
    --f[i];
  }
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

std::map< int, long long> ESU::bruteForce3(Hypergraph& inputGraph) {
  clearDataStruct();
  const int n = inputGraph.getNodeCount();
  vector<int> nodes;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      for (int z = j + 1; z < n; z++) {
        nodes = {i, j, z};
        Hypergraph motif = inputGraph.induceSubgraph(nodes);
        DisjointSet dsu(3);
        
        for (auto& e : motif.getIncidenceMatrix()) {
          for (int a = 0; a < (int) e.size(); a++) {
            for (int b = a + 1; b < (int) e.size(); b++) {
              dsu.join(e[a], e[b]);
            }
          }
        }
        if (dsu.tamanho[dsu.root(0)] != 3) continue;
        counterHyper[IsomorphismHyper::getLabel(motif)]++;
      }
    }
  }
  return counterHyper;
}

 //This method is used just for debug only!
std::map< int, long long> ESU::bruteForce4(Hypergraph& inputGraph) {
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
          for (auto& e : motif.getIncidenceMatrix()) {
            for (int a = 0; a < (int) e.size(); a++) {
              for (int b = a + 1; b < (int) e.size(); b++) {
                dsu.join(e[a], e[b]);
              }
            }
          }
          if (dsu.tamanho[dsu.root(0)] != 4) continue;
          counterHyper[IsomorphismHyper::getLabel(motif)]++;
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

std::map< int, long long> ESU::k3(Hypergraph& inputGraph) {
  clearDataStruct();
  for (auto& edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    visited.insert(edge); // it will insert the 3 nodes just visited
    counterHyper[IsomorphismHyper::getLabel(motif)]++;
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  vector< vector<int> > g = h.getGraph();
  startAlgo = steady_clock::now();
  setupAndRun(g, 3);
  return counterHyper;
}


/*
 * The following are a series of proposed methods for K_3
 */

// Auxiliary method to compute intermediate form

void ESU::k3IntermediateForm(Hypergraph& inputGraph) {
  for (auto& edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 3) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    if (motif.is_two_connected()) {
      Hypergraph simpleMotif = motif.filterEdge(2);
      --counterHyper[IsomorphismHyper::getLabel(simpleMotif)];
      // this occurence will be found by ESU!
    }
    counterHyper[IsomorphismHyper::getLabel(motif)]++;
  }
}


// TRIANGLE method

std::map< int, long long> ESU::k3Triangle(Hypergraph& inputGraph) {
  clearDataStruct();
  k3IntermediateForm(inputGraph);
  vector< vector<int> > g = inputGraph.getGraph();
  startAlgo = steady_clock::now();
  int n = (int) g.size();
  long long res = 0;
  long long triangle = 0;
  vector<int> deg(n);
  for (int i = 0; i < n; i++) {
    const int sz = g[i].size();
    res += 1LL * sz * (sz - 1) / 2 ;
    deg[i] = (int) g[i].size();
  }
  for (int i = 0; i < n; i++) {
    vector<int> bigger;
    for (auto& to : g[i]) {
      if (make_pair(deg[to], to) > make_pair(deg[i], i)) bigger.emplace_back(to);
    }
    for (int j = 0; j < (int) bigger.size(); j++) {
      for (int z = j + 1; z < (int) bigger.size(); z++) {
        if (inputGraph.validEdge({bigger[j], bigger[z]})) ++triangle;
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
  counterHyper[IsomorphismHyper::getLabel(h_trig)] += triangle;
  counterHyper[IsomorphismHyper::getLabel(h_line)] += res - 3 * triangle;
  return counterHyper;
}

// FASE method

std::map< int, long long> ESU::k3Fase(Hypergraph& inputGraph) {
  clearDataStruct();
  k3IntermediateForm(inputGraph);
  vector< pair<int, int> > edges;
  for (auto& edge : inputGraph.getIncidenceMatrix()) {
    if ((int) edge.size() == 2) edges.emplace_back(edge[0], edge[1]); 
  }
  startAlgo = steady_clock::now();
  auto census = FaSE(edges, 3);
  for (auto& [str, cnt] : census) {
    Hypergraph hyper = binaryToHyper(str, 3);
    counterHyper[IsomorphismHyper::getLabel(hyper)] += cnt;
  }
  return counterHyper;
}

// SIMPLE ESU

std::map< int, long long> ESU::k3ESU(Hypergraph& inputGraph) {
  clearDataStruct();
  k3IntermediateForm(inputGraph);
  Search = GRAPH;
  h = inputGraph.filterEdge(2);
  vector< vector<int> > g = inputGraph.getGraph();
  startAlgo = steady_clock::now();
  setupAndRun(g, 3);
  return counterHyper;
}


/*
 * Baseline method for K_4
 */

// implementar o caso K = 4 do paper
// https://github.com/HGX-Team/hypergraphx
// https://arxiv.org/pdf/2209.10241.pdf

std::map< int, long long> ESU::k4(Hypergraph& inputGraph) {
  clearDataStruct();
  for (auto& edge : inputGraph.getIncidenceMatrix()) { // assuming no duplicate edges ...
    if (edge.size() != 4) continue;
    Hypergraph motif = inputGraph.induceSubgraph(edge);
    visited.insert(edge); // it will insert the 3 nodes just visited
    counterHyper[IsomorphismHyper::getLabel(motif)]++;
  }
  Hypergraph reducedGraph = inputGraph.filterEdge(3); // at most 3 edges
  vector< vector<int> > g(reducedGraph.getNodeCount());
  int edge = 0;
  for (auto& e : reducedGraph.getIncidenceMatrix()) {
    for (auto& e_i : e) {
      g[e_i].emplace_back(edge);
    }
    ++edge;
  }
  edge = -1;
  for (auto& e : reducedGraph.getIncidenceMatrix()) {
    ++edge;
    if ( (int) e.size() != 3 ) continue;
    for (auto& node : e) {
      for (auto& e_i : g[node]) {
        vector<int> e1 = reducedGraph.getEdge(e_i);
        for (auto& add : reducedGraph.getEdge(edge)) {
          e1.emplace_back(add);
        }
        sort(e1.begin(), e1.end());
        e1.erase(unique(e1.begin(), e1.end()), e1.end());
        if (e1.size() == 4 && visited.find(e1) == visited.end()) {
          visited.insert(e1);
          Hypergraph motif = inputGraph.induceSubgraph(e1);
          counterHyper[IsomorphismHyper::getLabel(motif)]++;
        }
      }
    }
  }
  h = inputGraph;
  Search = HYPERGRAPH;
  g = h.getGraph();
  startAlgo = steady_clock::now();
  setupAndRun(g, 4);
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
    if (motif.is_two_connected()) {
      Hypergraph simpleMotif = motif.filterEdge(2);
      --counterHyper[IsomorphismHyper::getLabel(simpleMotif)];
      // this occurence will be found by ESU!
    }
    counterHyper[IsomorphismHyper::getLabel(motif)]++;
  }
  Hypergraph reducedGraph = inputGraph.filterEdge(3); // at most 3 edges
  vector< vector<int> > g(reducedGraph.getNodeCount());
  int edge = 0;
  for (auto& e : reducedGraph.getIncidenceMatrix()) {
    for (auto& e_i : e) g[e_i].emplace_back(edge);
    ++edge;
  }
  edge = -1;
  vector<int> f(reducedGraph.getNodeCount());
  std::stack<int> rem;
  for (auto& e : reducedGraph.getIncidenceMatrix()) {
    ++edge;
    if ( (int) e.size() != 3 ) continue;
    for (auto& node : e) {
      for (auto& e_i : g[node]) {
        vector<int> nodes = reducedGraph.getEdge(e_i);
        for (auto& add : reducedGraph.getEdge(edge)) nodes.emplace_back(add);
        sort(nodes.begin(), nodes.end());
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
        if (nodes.size() != 4) continue;
        int added_node = -1;
        for (auto& c : nodes) {
          if (find(e.begin(), e.end(), c) == e.end()) {
            added_node = c;
            break;
          }
        }
        if (f[added_node]) continue;
        f[added_node] = 1;
        rem.push(added_node);
        Hypergraph motif1 = inputGraph.induceSubgraphSkipComp(nodes);
        if (motif1.getEdgeMaxDeg() != 3) continue;
        int skip = 0;
        for (auto& motif_edge : motif1.getIncidenceMatrix()) {
          if (motif_edge.size() != 3) continue;
          if (motif_edge <= e) continue;
          int missing = -1;
          for (auto& n : nodes) {
            if (find(motif_edge.begin(), motif_edge.end(), n) == motif_edge.end()) {
              missing = n;
              break;
            }
          }
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
        if (motif.is_two_connected()) {
          Hypergraph simpleMotif = motif.filterEdge(2);
          --counterHyper[IsomorphismHyper::getLabel(simpleMotif)];
          // this occurence will be found by ESU!
        }
        counterHyper[IsomorphismHyper::getLabel(motif)]++;
      }
    }
    while (!rem.empty()) {
      f[rem.top()] = 0;
      rem.pop();
    }
  }
}

// FaSE
std::map< int, long long> ESU::k4Fase(Hypergraph& inputGraph) {
  clearDataStruct();
  k4IntermediateForm(inputGraph);
  vector< pair<int, int> > edges;
  for (auto& edge : inputGraph.getIncidenceMatrix()) {
    if ((int) edge.size() == 2) edges.emplace_back(edge[0], edge[1]); 
  }
  startAlgo = steady_clock::now();
  auto census = FaSE(edges, 4);
  for (auto& [str, cnt] : census) {
    Hypergraph hyper = binaryToHyper(str, 4);
    counterHyper[IsomorphismHyper::getLabel(hyper)] += cnt;
  }
  return counterHyper;
}

// SIMPLE ESU

std::map< int, long long> ESU::k4ESU(Hypergraph& inputGraph) {
  clearDataStruct();
  k4IntermediateForm(inputGraph);
  Search = GRAPH;
  h = inputGraph.filterEdge(2);
  vector< vector<int> > g = h.getGraph();
  startAlgo = steady_clock::now();
  setupAndRun(g, 4);
  return counterHyper;
}


/*
 * 
 * 
 */

void ESU::printResults(std::chrono::time_point<std::chrono::steady_clock> startTime, std::chrono::time_point<std::chrono::steady_clock> endTime, map< int, long long> subgraph_count, int k, bool detailedOutput, ostream& out) {
  // Remove keys with value 0 (only useful for display only)
  for(auto it = subgraph_count.begin(); it != subgraph_count.end(); ) {
    if(it -> second == 0) it = subgraph_count.erase(it);
    else ++it;
  }
  long long total_subgraph = 0;
  for (auto& [a, b] : subgraph_count) {
    total_subgraph += b;
  }
  out << "-----------------------------------------------" << endl;
  out << "Hyper-subgraphs: " << total_subgraph << " extracted" << endl;
  out << "Types: " << subgraph_count.size() << " found" << endl;
  out << "-----------------------------------------------" << endl;
  int counter = 0;
  for (auto& [a, b] : subgraph_count) {
    out << "Type #" << ++counter << ": " << std::setw(26) << std::right << b << " occurrences" << endl;
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
      h.printIncidenceMatrix(out);
      out << "Number of occurences: " << b << endl;
      out << "-----------------------------------------------" << endl;
    }
  }
}


void ESU::clearDataStruct() {
  subgraph.clear();
  counterHyper.clear();
  visited.clear();
}


void ESU::networkCensus(Hypergraph& h, int motifSize, bool detailedOutput, int algorithm, ostream& out) {
  cout << "Task: Network-census" << endl;
  cout << "Computing network-census on given network" << endl;
  cout << "-----------------------------------------------" << endl;
  clearDataStruct();
  int k;
  auto startTime = steady_clock::now();
  if (motifSize == 3) { // Execute our fastest method K = 3, TRIANGLE
    k = 3;
    switch (algorithm) {
      case 1: ESU::k3(h); break;
      case 2: ESU::k3ESU(h); break;
      case 3: ESU::k3Triangle(h); break;
      case 4: ESU::k3Fase(h); break;
      default: ESU::k3Triangle(h); break;
    }
  } else { // Execute our fastest method K=4, FASE
    k = 4;
    switch (algorithm) {
      case 1: ESU::k4(h); break;
      case 2: ESU::k4ESU(h); break;
      case 4: ESU::k4Fase(h); break;
      default: ESU::k4Fase(h); break;
    }
  }
  auto endTime = steady_clock::now();
  out << "Task completed in: " << duration_cast<duration<double>>(endTime - startAlgo).count() << " seconds" << endl;
  printResults(startTime, endTime, counterHyper, k, detailedOutput, out);
}
  
  
void ESU::findMotifs(Hypergraph& h, int motifSize, bool detailedOutput, bool significance_profile, int randomNetworks, int randomShuffles, int algorithm, ostream& out) {
  cout << "Task: Motif discovery" << endl;
  auto startTime = steady_clock::now();
  map<int, long long> census;
  if (motifSize == 3) {
    switch (algorithm) {
      case 1: census = ESU::k3(h); break;
      case 2: census = ESU::k3ESU(h); break;
      case 3: census = ESU::k3Triangle(h); break;
      case 4: census = ESU::k3Fase(h); break;
      default: census = ESU::k3Triangle(h); break;
    }
    // the numbers 0 - 5 are the internal labeling of all hyper-subgraphs with size = 3 
    // these numbers should cover all the possible values of IsomorphismHyper::precalc(3)
    for (int i = 0; i <= 5; i++) { 
      if (census.find(i) == census.end()) census[i] = 0; // if a subgraph didn't appear in the census => count is 0 
    }
  } else {
    switch (algorithm) {
      case 1: census = ESU::k4(h); break;
      case 2: census = ESU::k4ESU(h); break;
      case 4: census = ESU::k4Fase(h); break;
      default: census = ESU::k4Fase(h); break;
    }
    // the numbers 6 - 176 are the internal labeling of all hyper-subgraphs with size = 4 
    // these numbers should cover all the possible values of IsomorphismHyper::precalc(4)
    for (int i = 6; i <= 176; i++) {
      if (census.find(i) == census.end()) census[i] = 0; // if a subgraph didn't appear in the census => count is 0 
    }
  }
  map< int, vector<long long> > sample;
  cout << "Computing network-census on random networks" << endl;
  cout << "-----------------------------------------------" << endl;
  for (int i = 0; i < randomNetworks; i++) {
    cout << "Running on network: " << i + 1 << "/" << randomNetworks << endl;
    Hypergraph tmp = h;
    tmp.shuffleHypergraph(randomShuffles);
    map<int, long long> count;
    if (motifSize == 3) {
      switch (algorithm) {
        case 1: count = ESU::k3(tmp); break;
        case 2: count = ESU::k3ESU(tmp); break;
        case 3: count = ESU::k3Triangle(tmp); break;
        case 4: count = ESU::k3Fase(tmp); break;
        default: count = ESU::k3Triangle(tmp); break;
      }
    } else {
      switch (algorithm) {
        case 1: count = ESU::k4(tmp); break;
        case 2: count = ESU::k4ESU(tmp); break;
        case 4: count = ESU::k4Fase(tmp); break;
        default: count = ESU::k4Fase(tmp); break;
      }
    }
    for (auto& [a, b] : census) {
      sample[a].emplace_back(count[a]);
    }
  }
  double sum = 0;
  vector<double> score;
  cout << "-----------------------------------------------" << endl;
  if (significance_profile) { // sp
    cout << "Calculating the significance profile" << endl;
    cout << "-----------------------------------------------" << endl;
    for (auto& [a, b] : census) {
      double mean = 0;
      for (auto& value : sample[a]) mean += value;
      mean /= randomNetworks;
      score.emplace_back( (b - mean) / (b + mean + 4) );
      sum += score.back() * score.back();
    }
  } else {
    cout << "Calculating the z_score" << endl;
    cout << "-----------------------------------------------" << endl;
    vector<double> z_score;
    for (auto& [a, b] : census) {
      double mean = 0;
      for (auto& value : sample[a]) mean += value;
      mean /= randomNetworks;
      double std = 0;
      for (auto& value : sample[a]) std += (value - mean) * (value - mean);
      std /= randomNetworks - 1;
      std = sqrt(std);
      score.emplace_back((b - mean) / std);
      sum += score.back() * score.back();
    }
  }
  sum = sqrt(sum);
  for (auto& value : score) {
    if (sum > 0) value /= sum;
    else value = 0; // probably not enough networks were created / null model didn't work as intended ... set the value to 0
  }
  int counter = 0;
  out << "-----------------------------------------------" << endl;
  for (auto& [a, b] : census) {
    out << "Hyper-subgraph #" << counter + 1 << endl;
    if (detailedOutput) {
      auto adj = IsomorphismHyper::getHypergraph(a);
      Hypergraph h;
      h.setIncidenceMatrix(adj);
      h.setN(motifSize);
      h.printIncidenceMatrix(out);
    }
    out << "Number of occurences: " << b << endl;
    out << (significance_profile ? "Significance Profile: " : "Z score: ") << score[counter++] << endl;
    out << "-----------------------------------------------" << endl;
  }
  auto endTime = steady_clock::now();
  out << "Task completed in: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
  out << "-----------------------------------------------" << endl;
}
