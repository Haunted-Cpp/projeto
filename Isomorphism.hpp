#pragma once

class Isomorphism {
private:
  
public:
   static bool isomorphismSlow(const Hypergraph& h1, const Hypergraph& h2);
   static bool isomorphismNauty(const Hypergraph& h1, const Hypergraph& h2);
   static std::vector<graph> canonization(const Hypergraph& h);
};
