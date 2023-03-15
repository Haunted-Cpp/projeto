#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>

#include "nauty.h"
#include "Hypergraph.hpp"
#include "Isomorphism.hpp"

using namespace std;

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);
  Hypergraph h1;
  Hypergraph h2;
  cout << Isomorphism::isomorphismNauty(h1, h2) << '\n';
  return 0;
}

