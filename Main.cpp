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

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);
  //Hypergraph h;
  //h.readIncidenceMatrix();
  //auto test = ESU::startEdgeGraphSubgraphs(h, 3);
  //cout << test.size() << '\n';
  //for (auto t : test) {
    //h.printEdgeSubgraph(t);
  //}
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
  
  for (auto [a, b] : ESU::getEquivalenceClass(g, 5)) {
    cout << a << ' ' << b << '\n';
  }
  return 0; 
}

