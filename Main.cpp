#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <set>
#include <stack>
#include <map>

#include "nauty.h"

#include "Hypergraph.hpp"
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
      //h1.printGraph();
      //h2.printGraph();
      exit(0);
    }
    if (a == b) cout << "Running on test: " << i << " " << a << ' ' << b << '\n';
    else {
     cout << "WA!!!" << '\n';
     //h1.printGraph();
     cout << '\n';
     //h2.printGraph();
     break; 
    }
  }
  cout << same << ' ' << diff << '\n';
}

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);
  
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
  auto counter = ESU::getAllSubgraphs(g, 6);
  cout << (int) counter.size() << '\n';
  auto x = ESU::getAllSubgraphs(g, 5);
  cout << (int) x.size() << '\n';
  //for (auto sub) {
    //cout << a << ' ' << b << '\n';
  //}
  //Hypergraph h;
  //h.readIncidenceMatrix();
  //ESU::startEdgeGraphSubgraphs(h, 2);
  //ESU::startEdgeGraphSubgraphs(h, 3);
  //ESU solve;
  //solve.start(g, 1);
  //solve.start(g, 2);
  //solve.start(g, 3);
  //ESU::start(g, 4);
  //ESU::start(g, 5);
  //ESU::start(g, 6);
  //solve.start(g, 7);
  //Hypergraph h1;
  //h1.readIncidenceMatrix();
  //h1.ESU(2);
  //h1.buildEdgeGraph();
  return 0; 
}

