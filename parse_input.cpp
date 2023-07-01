#include <bits/stdc++.h>
using namespace std;

//#define DISPLAY_STATS

set<int> parse_str(string nodes) {
  string node;
  istringstream token(nodes);
  set<int> arr;
  while (token >> node) {
    int n = stoi(node);
    if (arr.find(n) != arr.end()) {
      continue;
    }
    assert(arr.find(n) == arr.end());
    arr.insert(n);
  }
  return arr;
}

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);
  string line;
  getline(cin, line);
  auto nodes = parse_str(line);
  vector<int> counter(10);
  vector< vector<int> > edges;
  set< vector<int> > is_unique;
  
  while (getline(cin, line)) {
    auto edge = parse_str(line);
    int size = (int) edge.size(); 
    if (size > 1 && size <= 4) {
      vector<int> hyperedge;
      hyperedge.assign(edge.begin(), edge.end());
      if(is_unique.find(hyperedge) != is_unique.end()) {
        continue;
      }
      assert(is_unique.find(hyperedge) == is_unique.end());
      is_unique.insert(hyperedge);
      ++counter[size];
      edges.emplace_back(hyperedge);
    }
  }
  //cout << (int) nodes.size() << '\n';
  //cout << (int) edges.size() << '\n';
  assert( (int) edges.size() == accumulate(counter.begin(), counter.end(), 0) );
  map<int, int> number;
  int current_number = 0;
  for (auto edge : edges) {
    //cout << (int) edge.size() << ' ';
    for (auto e : edge) {
      if (number.find(e) == number.end()) {
        number[e] = ++current_number;
      }
      cout << number[e] << ' ';
    }
    cout << '\n';
  }
  #ifdef DISPLAY_STATS
    cout << "------------------------------------------------------------" << '\n';
    cout << "Number of nodes: " << nodes.size() << '\n';
    cout << "Edge Sizes: " << '\n';
    for (int i = 2; i <= 4; i++) cout << counter[i] << '\n';
  #endif
  return 0;
}
