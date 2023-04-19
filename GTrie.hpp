#pragma once

class GTrie {
private:
  
  struct Node {
    std::vector< std::shared_ptr<Node> > children;
    std::string edgeLink;
    bool isLeaf;
  };
  
  std::shared_ptr<Node> root;
  const int K;
  void dfs(std::shared_ptr<Node>);
  std::string getEncoding(Hypergraph&, const std::vector<int>&, int);
  std::vector< std::vector<int> > nei;
  std::vector< std::vector< std::vector<int> > > bijection;
  void match(Hypergraph&, std::shared_ptr<Node>, int);
public:
  GTrie(int k);
  void insert(Hypergraph&);
  void dfs();
  void search(Hypergraph&);
};
