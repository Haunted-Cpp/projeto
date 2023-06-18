#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <set>
#include <stack>
#include <map>
#include <memory>
#include <chrono>
#include <cmath>

#include "nauty.h"

#include "Hypergraph.hpp"
#include "GTrie.hpp"
#include "Isomorphism.hpp"
#include "ESU.hpp"

using namespace std;
using namespace std::chrono;

void testIsomorphism(int tt = 100'000'000) {
  int diff = 0;
  int same = 0;
  for (int i = 0; i < tt; i++) {
    Hypergraph h1;
    Hypergraph h2;
    int a = Isomorphism::isomorphismNauty(h1, h2);
    int b = Isomorphism::isomorphismSlow(h1, h2);
    same += a == 1;
    diff += a == 0;
    if (same == 1) {
      exit(0);
    }
    if (a == b) cout << "Running on test: " << i << " " << a << ' ' << b << '\n';
    else {
     cout << "WA!!!" << '\n';
     cout << '\n';
     break; 
    }
  }
  cout << same << ' ' << diff << '\n';
}

//Counter: 5
//885
//363
//187
//83
// 10

void findMotifs3() {
  Hypergraph h;
  h.readFromStdin();
  auto census = ESU::bruteForce3(h);
  map< vector<graph>, vector<int> > sample;
  const int NUMBER_NETWORKS = 100;
  //for (auto [a, b] : census) {
    //cout << b << '\n';
  //}
  //exit(0);
  for (int i = 0; i < NUMBER_NETWORKS; i++) {
    h.shuffleHypergraph(20);
    auto count = ESU::bruteForce3(h);
    //exit(0);
    for (auto [a, b] : census) {
      //cout << b << '\n';
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
    //cout << mean << ' ' << std << '\n';
    cout << b << ' ' << (b - mean) / std << '\n';
  }
}

void readHypergraph() {
  Hypergraph h;
  h.readFromFile("Dataset/ps.in");
  
  //auto g = Isomorphism::canonization(h);
  //for (auto value : g) {
    //cout << value << '\n';
  //}
  //cout << h.getNodeCount() << ' ' << g.size() << '\n';
  //exit(0);
  
  {
    ESU::clearDataStruct();
    auto startTime = steady_clock::now();
    auto subgraph_count = ESU::k3Modified(h);
    auto endTime = steady_clock::now();
    cout << "-----------------------------------------------" << '\n';
    cout << "Network census completed in: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
    cout << subgraph_count.size() << " types of hyper-subgraphs found" << '\n';
    cout << "-----------------------------------------------" << '\n';
    int counter = 0;
    for (auto [a, b] : subgraph_count) {
      cout << "Type #" << ++counter << ": " << b << '\n';
    }
    cout << "-----------------------------------------------" << '\n';
    counter = 0;
    for (auto [a, b] : subgraph_count) {
      auto adj = Isomorphism::getHypergraph(a);
      Hypergraph h;
      h.setIncidenceMatrix(adj);
      cout << "Hyper-subgraph #" << ++counter << '\n';
      h.printIncidenceMatrix();
      cout << "Number of occurences: " << b << '\n';
      cout << "-----------------------------------------------" << '\n';
    }
    
    
  }
  
  exit(0);
  
  {
    ESU::clearDataStruct();
    auto startTime = steady_clock::now();
    ESU::k3(h);
    auto endTime = steady_clock::now();
    cout << "Time: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
  }
  //cout << "----" << '\n';
  //exit(0);
  
}

void readNormal() {
  int n, m;
  cin >> n >> m;
  vector< vector<int> > g(n);
  
  std::set< pair<int, int> > vis;
  for (int i = 0; i < m; i++) {
    int st, et, w;
    cin >> st >> et >> w;
    --st; --et;
    if (st == et) continue;
    assert(st != et);
    if (vis.find({st, et}) == vis.end()) {
      vis.insert({st, et});
      vis.insert({et, st});
    } else {
      continue;
    }
    g[st].emplace_back(et);
    g[et].emplace_back(st);
  }
  //exit(0);
  //cout << "YES" << '\n';
  
  for (auto [a, b] : ESU::getEquivalenceClass(g, 3)) {
    cout << a << ' ' << b << '\n';
  }
  
  //for (auto [a, b] : ESU::getEquivalenceClass(g, 3)) {
    //cout << a << ' ' << b << '\n';
  //}
}

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);
  //findMotifs3();
  //readNormal();
  readHypergraph();
  //readHypergraph();
  
  return 0; 
}

//55574541

//55574537
//duarte@duarte-pc:~/Desktop/Projeto/Código$ ./Main
//55574537
//16659790
//12
//4
//VISITED: 16659802
//Time: 247.627 seconds
//----
//duarte@duarte-pc:~/Desktop/Projeto/Código$ make
//g++ -O3 -std=c++17 -c -o Main.o Main.cpp
//g++ -O3 -std=c++17  -o Main nauty/nauty.h nauty/nauty.c nauty/nautil.c nauty/naugraph.c nauty/schreier.c nauty/naurng.c ESU.o GTrie.o Hypergraph.o Isomorphism.o Main.o
//duarte@duarte-pc:~/Desktop/Projeto/Código$ ./Main
//Counter: 4
//55574537
//16659790
//12
//4
//Time: 79.7647 seconds



