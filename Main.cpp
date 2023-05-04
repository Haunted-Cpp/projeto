#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <set>
#include <stack>
#include <map>
#include <memory>
#include <chrono>

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
//10

void readHypergraph() {
  Hypergraph h;
  //h.readFromFile("bug.in");
  //h.printIncidenceMatrix();
  //exit(0);
  //h.readFromFile("test_bug_hyper.in");
  h.readFromStdin();
  //for (int i = 0; i < 10000; i++) {
    //string filename = "test_input" + to_string(i) + ".in";
    //h.randomHypergraph(25, 10000, 6);
    //h.saveToFile("buggy");
    
    //cout << i << endl;
    
    //ESU::bruteForce3(h);
    //ESU::k3Modified(h);
    //ESU::k3(h);
  //}
  {
    auto startTime = steady_clock::now();
    ESU::k3(h);
    auto endTime = steady_clock::now();
    cout << "Time: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
  }
  {
    auto startTime = steady_clock::now();
    ESU::k3Modified(h);
    auto endTime = steady_clock::now();
    cout << "Time: " << duration_cast<duration<double>>(endTime - startTime).count() << " seconds" << endl;
  }
    //ESU::bruteForce4(h);
  //}
  //ESU::k3(h);
  //h.readFromFile("test_input.in");
  //h.saveToFile("test_output.out");
  //cout << h.getEdgeCount() << '\n';
  //h.randomHypergraph();
  //h.saveToFile("gen_hypergraph.txt");
  //h.printIncidenceMatrix();
  //h.readIncidenceMatrix();
  
  //Hypergraph subgraph;
  //subgraph.readIncidenceMatrix();
  //Hypergraph h1;
  //h1.readIncidenceMatrix();
  //GTrie trie(4);
  //trie.insert(subgraph);
  //Hypergraph h2;
  //h2.readIncidenceMatrix();
  
  //Hypergraph h3;
  //h3.readIncidenceMatrix();
  
  //trie.insert(h1);
  //trie.insert(h2);
  //trie.insert(h3);
  //trie.dfs();
  //trie.dfs();
  //trie.search(h);
  //ESU::k4(h);
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
    if (vis.find({st, et}) == vis.end()) {
      vis.insert({st, et});
      vis.insert({et, st});
    } else {
      assert(false);
    }
    g[st].emplace_back(et);
    g[et].emplace_back(st);
  }
  
  for (auto [a, b] : ESU::getEquivalenceClass(g, 3)) {
    cout << a << ' ' << b << '\n';
  }
}

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);
  readHypergraph();
  //readNormal();
  return 0; 
}

