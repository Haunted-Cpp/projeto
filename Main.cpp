#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set>
//#include <set>
#include <stack>
#include <map>
#include <memory>
#include <chrono>
#include <cmath>
#include <queue>

#include "nauty.h"

#include "Hypergraph.hpp"
#include "GTrie.hpp"
#include "IsomorphismHyper.hpp"
#include "ESU.hpp"



#include "FaSE/Fase.h"
#include "FaSE/DynamicGraph.h"
#include "FaSE/GraphMatrix.h"
#include "FaSE/GraphUtils.h"


void findMotifs3() {
  Hypergraph h;
  h.readFromStdin();
  auto census = ESU::bruteForce3(h);
  map< vector<graph>, vector<int> > sample;
  const int NUMBER_NETWORKS = 100;
  for (int i = 0; i < NUMBER_NETWORKS; i++) {
    h.shuffleHypergraph(20);
    auto count = ESU::bruteForce3(h);
    for (auto [a, b] : census) {
      sample[a].emplace_back(count[a]);
    }
  }
  for (auto [a, b] : census) {
    double mean = 0;
    for (auto value : sample[a]) {
      mean += value;
    }
    mean /= NUMBER_NETWORKS;
    double std = 0;
    for (auto value : sample[a]) {
      std += (value - mean) * (value - mean);
    }
    std /= NUMBER_NETWORKS - 1;
    std = sqrt(std);
    cout << b << ' ' << (b - mean) / std << '\n';
  }
}

void printResults
(
std::chrono::time_point<std::chrono::steady_clock> startTime, 
std::chrono::time_point<std::chrono::steady_clock> endTime, 
map< vector<graph>, long long> subgraph_count, 
int k,
int printDetail = 0
) {
  cout << "-----------------------------------------------" << '\n';
  cout << "Network census completed in: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
  cout << subgraph_count.size() << " types of hyper-subgraphs found" << '\n';
  cout << "-----------------------------------------------" << '\n';
  int counter = 0;
  for (auto [a, b] : subgraph_count) {
    cout << "Type #" << ++counter << ": " << b << '\n';
  }
  cout << "-----------------------------------------------" << '\n';
  if (printDetail) {
    counter = 0;
    for (auto [a, b] : subgraph_count) {
      auto adj = IsomorphismHyper::getHypergraph(a);
      Hypergraph h;
      h.setIncidenceMatrix(adj);
      h.setN(k);
      cout << "Hyper-subgraph #" << ++counter << '\n';
      h.printIncidenceMatrix();
      cout << "Number of occurences: " << b << '\n';
      cout << "-----------------------------------------------" << '\n';
    }
  }
  cout << flush;
}

void readHypergraph() {
  Hypergraph h;
  h.readFromFile("Dataset/dblp.in");
  int k = 3;
  {
    auto startTime = steady_clock::now();
    auto subgraph_count = ESU::k3Modified(h);
    auto endTime = steady_clock::now();
    printResults(startTime, endTime, subgraph_count, 3);
  }
  {
    auto startTime = steady_clock::now();
    auto subgraph_count = ESU::k3(h);
    auto endTime = steady_clock::now();
    printResults(startTime, endTime, subgraph_count, 3);
    
  }
  //{
    //auto startTime = steady_clock::now();
    //auto subgraph_count = ESU::bruteForce3(h);
    //auto endTime = steady_clock::now();
    //printResults(startTime, endTime, subgraph_count, 3);
    
  //}
}

void readNormal() {
  int n, m;
  cin >> n >> m;
  vector< vector<int> > g(n);
  
  //std::set< pair<int, int> > vis;
  for (int i = 0; i < m; i++) {
    int st, et, w;
    cin >> st >> et >> w;
    --st; --et;
    if (st == et) continue;
    assert(st != et);
    //if (vis.find({st, et}) == vis.end()) {
      //vis.insert({st, et});
      //vis.insert({et, st});
    //} else {
      //continue;
    //}
    g[st].emplace_back(et);
    g[et].emplace_back(st);
  }
  for (auto [a, b] : ESU::getEquivalenceClass(g, 3)) {
    cout << a << ' ' << b << '\n';
  }
}

//int main() {
  //std::ios::sync_with_stdio(0);
  //std::cin.tie(0);
  //readHypergraph();
  //return 0; 
//}

Graph *G;
int K = 0;
bool dir = false, largeScale = true;

void read()
{
  largeScale = true;
  if (largeScale) G = new DynamicGraph();
  else G = new GraphMatrix();

  bool zeroBased = false;
  
  vector<pair<int, int> > edges;
  
  int a, b;
  while (cin >> a >> b) {
    edges.emplace_back(a, b);
  }
  
  GraphUtils::readFile(G, edges, dir, false, zeroBased);
  G->sortNeighbours();
  G->makeArrayNeighbours();

  // Subgraph Size
  K = 3;
}

void output(Fase* fase)
{
  //printf("Finished Calculating\n");
  FILE *f = stdout;
  
  
  //if (largeScale) fprintf(f, "Graph Representation: Large Scale\n");
  //else fprintf(f, "Graph Representation: Adjacency Matrix\n");
  //fprintf(f, "\nExact Enumeration, no Sampling done\n");
  //fprintf(f, "\n\tDetailed Output:\n");
  for (auto element : fase->subgraphCount()) {
    fprintf(f, "%s: %d occurrences\n", element.second.c_str(), element.first);
  }
}

void finish(Fase* fase)
{
  delete fase;
  delete G;
}

int main(int argc, char **argv)
{
  read();
  Random::init(time(NULL));
  Fase* fase = new Fase(G, dir);
  auto startTime = steady_clock::now();
  fase->runCensus(K);
  auto endTime = steady_clock::now();
  cout << "Time: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
  output(fase);
  finish(fase);
  return 0;
}



