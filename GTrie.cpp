#include <bits/stdc++.h>

#include "Hypergraph.hpp"
#include "GTrie.hpp"
#include "Settings.hpp"

using namespace std;


GTrie::GTrie(int k) : K(k) { // Limit on MAX EDGE size
  assert(k <= 10); // Limit the MAX EDGE SIZE
  root = make_shared<Node>();
  bijection.clear();
  bijection.emplace_back(vector< vector<int> >()); // level 0 is empty ...
  for (int i = 1; i <= MAX_HYPER_MOTIF_SIZE; i++) { // The Hyper Motif SHOULD NOT be too big ...
    bijection.emplace_back(vector< vector<int> >());
    vector<int> perm;
    auto Gen = [&](auto&& Gen, int len, int last) {
      if (len == 0) {
        if (perm.empty()) return;
        perm.emplace_back(i - 1);
        bijection[i].emplace_back(perm);
        perm.pop_back( );
        return;
      }
      for (int k = last + 1; k < i - 1; k++) {
        perm.emplace_back(k);
        Gen(Gen, len - 1, k);
        perm.pop_back();
      }
    };
    bijection[i].push_back({i - 1});
    for (int j = 0; j < i - 1; j++) {
      bijection[i].push_back({j, i - 1});
    }
    for (int j = 2; j <= k - 1; j++) Gen(Gen, j, -1); // For edge sizes greater than 2
  }
  for (int i = 1; i <= MAX_HYPER_MOTIF_SIZE; i++) {
    cout << "Level: " << i << '\n';
    for (auto tot : bijection[i]) {
      for (auto vec : tot) cout << vec << ' ';
      cout << '\n';
    }
    cout << '\n';
  }
}

/*
 * Given a set of nodes, produce the labelling at depth i.
 * The vector nodes is used to allow orderings different from 0 .. n - 1
 * The last vertex is always the one to be added
 */


string GTrie::getEncoding(Hypergraph& subHyperGraph, const vector<int>& nodes, int depth) {
  assert( (int) nodes.size() == depth );
  string encode;
  for (auto values : bijection[depth]) {
    vector<int> ordered_nodes;
    for (auto index : values) ordered_nodes.emplace_back(nodes[index]);
    sort(ordered_nodes.begin(), ordered_nodes.end());
    encode += (subHyperGraph.validEdge(ordered_nodes) ? '1' : '0');
  }
  return encode;
}

void GTrie::insert(Hypergraph& subHyperGraph) {
  assert( (int) subHyperGraph.getNodeCount() <= 10 );
  assert(subHyperGraph.getEdgeMaxDeg() <= K);
  shared_ptr<Node> node = root;
  for (int i = 0; i < subHyperGraph.getNodeCount(); i++) {
    vector<int> nodes(i + 1); 
    iota(nodes.begin(), nodes.end(), 0);
    string encode = getEncoding(subHyperGraph, nodes, i + 1);
    int j = 0;
    for (; j < (int) node -> children.size(); j++) {
      if (node -> children[j] -> edgeLink == encode) {
        break;
      }
    }
    if (j == (int) node -> children.size() ) {
      node -> children.emplace_back(make_shared<Node>());
      node -> children.back() -> edgeLink = encode;
    }
    node = node -> children[j];
  }
  node -> isLeaf = true;
}

vector<int> Vused;
vector<int> inVused;

void GTrie::match(Hypergraph& h, shared_ptr<Node> gtrie, int len) {
  vector<int> Vcand;
  if (Vused.empty()) { // every vertex is a valid candidate
    Vcand = vector<int>( (int) nei.size() );
    iota(Vcand.begin(), Vcand.end(), 0);
  } else {
    vector<int> Vconn = {}; // find all the nodes that the current one MUST be connected
    assert( (int) gtrie -> edgeLink.size() == (int) bijection[len + 1].size() );
    for (int i = 0; i < (int) gtrie -> edgeLink.size(); i++) {
      if (gtrie -> edgeLink[i] == '1') {
        for (int j = 0; j < (int) bijection[len + 1][i].size() - 1; j++) { // -1 because we don't consider the last node - it's not added yet
          assert(bijection[len + 1][i][j] < (int) Vused.size());
          Vconn.emplace_back(Vused[ bijection[len + 1][i][j] ]);
        }
      }
    }
    if (Vconn.empty()) {
      return; // the subhypergraph should be induced
    }
    // Find the Vconn with the SMALLEST deg.
    int m = Vconn[0];
    for (auto& node : Vconn) {
      if ( (int) nei[node].size() < (int) nei[m].size() ) {
        m = node;
      }
    }
    for (auto& node : nei[m]) if (!inVused[node]) Vcand.emplace_back(node);
  }
  for (auto& node : Vcand) {
    inVused[node] = true;
    Vused.emplace_back(node);
    string encode = getEncoding(h, Vused, len + 1);
    if (encode == gtrie -> edgeLink) {
      if (gtrie -> isLeaf) {
        cout << "FOUND GRAPH: " << '\n';
        for (auto& report : Vused) {
          cout << report << ' ';
        }
        cout << '\n';
      } else {
        for (auto& children : gtrie -> children) {
          match(h, children, len + 1);
        }
      }
      
    }
    inVused[Vused.back()] = false;
    Vused.pop_back();
  }
}

void GTrie::search(Hypergraph& h) { // Search in h all the subhypergraphs found so far
  nei = h.buildVertexGraph(K);
  for (auto child : root -> children) {
    inVused = vector<int>(h.getNodeCount());
    Vused.clear();
    match(h, child, 0);
  }
}

void GTrie::dfs(shared_ptr<Node> node) {
  if (node == nullptr) return;
  if (!node -> edgeLink.empty()) cout << "Current matrix: " << node -> edgeLink << '\n';
  for (auto nei : node -> children) {
    dfs(nei);
  }
  cout << "Leaving ..." << '\n';
}

void GTrie::dfs() {
  dfs(root);
}
