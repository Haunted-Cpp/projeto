#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>

#include "nauty.h"
#include "Hypergraph.hpp"
#include "Isomorphism.hpp"

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
  Hypergraph h1;
  h1.readIncidenceMatrix();
  h1.ESU(2);
  //h1.buildEdgeGraph();
  return 0; 
}

