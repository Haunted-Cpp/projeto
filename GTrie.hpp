#pragma once

class GTrie {
private:
  
  struct Node {
    std::vector< std::shared_ptr<Node> > children;
    std::vector<std::string> edgeLink;
    bool isLeaf;
  };
  std::shared_ptr<Node> root;
  const int K;
  void dfs(std::shared_ptr<Node>);
public:
  GTrie(int k);
  void insert(Hypergraph&);
  void dfs();
};
