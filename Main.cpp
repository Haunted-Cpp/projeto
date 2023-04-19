#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <set>
#include <stack>
#include <map>
#include <memory>

#include "nauty.h"

#include "Hypergraph.hpp"
#include "GTrie.hpp"
#include "Isomorphism.hpp"
#include "ESU.hpp"

using namespace std;

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

void readHypergraph() {
  Hypergraph h;
  h.readIncidenceMatrix();
  
  Hypergraph subgraph;
  subgraph.readIncidenceMatrix();
  //Hypergraph h1;
  //h1.readIncidenceMatrix();
  GTrie trie(2);
  trie.insert(subgraph);
  //Hypergraph h2;
  //h2.readIncidenceMatrix();
  
  //Hypergraph h3;
  //h3.readIncidenceMatrix();
  
  //trie.insert(h1);
  //trie.insert(h2);
  //trie.insert(h3);
  //trie.dfs();
  //trie.dfs();
  trie.search(h);
  //ESU::k4(h);
}

void readNormal() {
  int n, m;
  cin >> n >> m;
  vector< vector<int> > g(n);
  for (int i = 0; i < m; i++) {
    int st, et, w;
    cin >> st >> et >> w;
    --st; --et;
    g[st].emplace_back(et);
    g[et].emplace_back(st);
  }
  
  for (auto [a, b] : ESU::getEquivalenceClass(g, 6)) {
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

