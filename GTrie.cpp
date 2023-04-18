#include <bits/stdc++.h>

#include "Hypergraph.hpp"
#include "GTrie.hpp"

using namespace std;


GTrie::GTrie(int k) : K(k) { // Limit on MAX EDGE size
  root = make_shared<Node>();
}

void GTrie::insert(Hypergraph& h) {
  assert(h.getEdgeMaxDeg() <= K);
  vector< vector<int> > edge = h.getIncidenceMatrix();
  shared_ptr<Node> node = root;
  for (int i = 0; i < h.getNodeCount(); i++) {
    string encode = "";
    
    { // Build encoding
      
      vector<int> previousNodes;
      auto findEncoding = [&](auto&& findEncoding, int len) {
        if (len == 0) {
          previousNodes.emplace_back(i);
          
          for (auto n : previousNodes) cout << n << ' ';
          cout << '\n';
          
          sort(previousNodes.begin(), previousNodes.end());
          if (h.validEdge(previousNodes)) encode += '1';
          else encode += '0';
          previousNodes.pop_back();
          return;
        }
        for (int j = (previousNodes.empty() ? 0 : previousNodes.back() + 1); j < i; j++) {
          previousNodes.emplace_back(j);
          findEncoding(findEncoding, len - 1);
          previousNodes.pop_back();
        }
      };
      findEncoding(findEncoding, 2 - 1);
      encode += (h.validEdge({i, i}) ? '1' : '0');
      encode += '#';
      for (int j = 3; j <= K; j++) {
        findEncoding(findEncoding, j - 1);
        encode += '#';
      }
    }
    //if (i == 0) encode += (h.validEdge({0, 0}) ? '1' : '0');
    //cout << encode << '\n';
    //if (i == 1) return;
    
    int j = 0;
    for (; j < (int) node -> edgeLink.size(); j++) {
      if (node -> edgeLink[j] == encode) {
        break;
      }
    }
    if (j == (int) node -> edgeLink.size() ) {
      node -> children.emplace_back(make_shared<Node>());
      node -> edgeLink.emplace_back(encode);
    }
    node = node -> children[j];
  }
}

void GTrie::dfs(shared_ptr<Node> node) {
  if (node == nullptr) return;
  for (int i = 0; i < (int) node -> edgeLink.size(); i++) {
    cout << "Moving along edge: " << node -> edgeLink[i] << '\n';
    dfs(node -> children[i]);
    cout << "Leaving edge: " << '\n';
  }
}

void GTrie::dfs() {
  dfs(root);
}
