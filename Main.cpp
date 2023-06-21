#include <bits/stdc++.h>

#include "nauty.h"
#include "Hypergraph.hpp"
#include "GTrie.hpp"
#include "IsomorphismHyper.hpp"
#include "ESU.hpp"
#include "FaSE/Fase.h"
#include "FaSE/DynamicGraph.h"
#include "FaSE/GraphMatrix.h"
#include "FaSE/GraphUtils.h"

int main(int argc, char* argv[]) {
  std::ios::sync_with_stdio(0);
  std::cin.tie(0);
  
  int motifSize = -1, n = -1;
  string inputFile, outputFile, task;
  bool detailedOutput = false;
  
  
  vector<string> args;
  for (int i = 1; i < argc; i++) {
    string arg = string(argv[i]);
    if (arg == "-d") detailedOutput = true;
    else args.emplace_back(arg);
  }
  
  if ( args.size() & 1 ) { // If number of arguments is odd, then a pair is missing
    cout << "Invalid number of arguments provided." << endl;
    cout << "Each flag must be followed by an argument" << endl;
    cout << "Use the format: -flag argument" << endl;
    cout << "More information provide in the readme file" << endl;
    return 0;
  }
  
  for (int i = 0; i < (int) args.size() - 1; i += 2) {
    if (args[i] == "-s") {
      string number = args[i + 1];
      if (any_of(number.begin(), number.end(), [](char c) { return !isdigit(c); })) {
        cout << "No valid motifSize provided with -s flag" << endl;
        cout << "Use \"-s <k>\", where 3 <= k <= 4" << endl;
        return 0;
      }
      motifSize = std::stoi(args[i + 1]);
    } else if (args[i] == "-i") {
      inputFile = args[i + 1];
    } else if (args[i] == "-o") {
      outputFile = args[i + 1];
    } else if (args[i] == "-m") {
      task = args[i + 1];
    } else if (args[i] == "-n") {
      string number = args[i + 1];
      if (any_of(number.begin(), number.end(), [](char c) { return !isdigit(c); })) {
        cout << "No valid number of nodes provided with -n flag" << endl;
        cout << "Use \"-n <k>\", \n k should be as large as the number of nodes used." << endl;
        return 0;
      }
      n = std::stoi(args[i + 1]);
    } else {
      cout << "Invalid argument flag provided." << endl;
      cout << "More information in the readme provided" << endl;
      return 0;
    }
  }
  
  // -s **size**, is mandatory
  if (motifSize == -1 || motifSize < 3 || motifSize > 4) {
    cout << "No valid motifSize provided with -s flag" << endl;
    cout << "Use \"-s <k>\", where 3 <= k <= 4" << endl;
    return 0; 
  } 
  
  // -i **inputFile**, is mandatory
  if (inputFile.empty()) {
    cout << "No input file provided with -i flag" << endl;
    return 0; 
  }
  
  // -m **method**, is mandatory
  if (task != "count" && task != "motif") {
    cout << "No valid method provided with -m flag." << endl;
    cout << "Valid options are: \"count\", \"motif\"" << endl;
    return 0; 
  }
  
  std::ofstream fout;
  std::ostream& out = !outputFile.empty() ? fout : std::cout;
  
  
  if (!outputFile.empty()) {
    fout.open(outputFile);
    if (fout.fail()) {
      cout << "output file: " << outputFile << " - could not be opened" << endl;
      exit(EXIT_FAILURE);
    };
  }
  
  // https://oeis.org/A323817
  IsomorphismHyper::precalc(3);
  IsomorphismHyper::precalc(4);
  IsomorphismHyper::no_use = 1;
  
  //cout << "-------" << endl;
  Hypergraph h(n); 
  h.readFromFile(inputFile);

  
  
  bool significanceProfile = true;
  int randomNetworks = 100;
  int randomShuffles = 1000;
  
  
  out << "Hypergraph read from file: " << inputFile << endl;
  out << "-----------------------------------------------" << endl;
  out << "Number of nodes: " << h.getNodeCount() << endl;
  out << "Number of hyperedges: " << h.getEdgeCount() << endl;
  out << "-----------------------------------------------" << endl;
  vector<int> size = h.getEdgeBySize();
  for (int i = 2; i <= MAX_EDGE_SIZE; i++) {
    out << "E" << i << ": " << size[i] << endl;
  }
  out << "-----------------------------------------------" << endl;
  
  if (task == "count") { // network-census of whole network
    ESU::networkCensus(h, motifSize, detailedOutput, out);
  } else { // find significance profile of each subgraph
    ESU::findMotifs(h, motifSize, detailedOutput, significanceProfile, randomNetworks, randomShuffles, out);
  }
  
  fout.close();
  
  return 0; 
}








