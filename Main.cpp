#include <bits/stdc++.h>

#include "nauty.h"
#include "Hypergraph.hpp"
#include "GTrie.hpp"
#include "IsomorphismHyper.hpp"
#include "CountingMethods.hpp"
#include "FaSE/Fase.h"
#include "FaSE/DynamicGraph.h"
#include "FaSE/GraphMatrix.h"
#include "FaSE/GraphUtils.h"

int main(int argc, char* argv[]) {
  std::ios::sync_with_stdio(0);
  std::cin.tie(0);
  
  int motifSize = -1, n = -1, randomNetworks = 1000, randomShuffles = 1000, algorithm = -1;
  string inputFile, outputFile = "report.txt";
  bool detailedOutput = false, significanceProfile = true, motif = false;
  
  vector<string> args;
  for (int i = 1; i < argc; i++) {
    string arg = string(argv[i]);
    if (arg == "-d") detailedOutput = true;
    else if (arg == "-z") significanceProfile = false; // use z_score instead
    else if (arg == "-m") motif = true;
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
    } else if (args[i] == "-a") {
      if (any_of(args[i + 1].begin(), args[i + 1].end(), [](char c) { return !isdigit(c); })) {
        cout << "No valid algorithm provided with -a flag" << endl;
        cout << "Use \"-a <k>\", \n 1 <= k <= 4" << endl;
        return 0;
      }
      algorithm = std::stoi(args[i + 1]);
      if (algorithm > 4 || algorithm < 1) {
        cout << "No valid algorithm provided with -a flag" << endl;
        cout << "Use \"-a <k>\", \n 1 <= k <= 4" << endl;
        return 0;
      }
    } else if (args[i] == "-n") {
      string number = args[i + 1];
      if (any_of(number.begin(), number.end(), [](char c) { return !isdigit(c); })) {
        cout << "No valid number of nodes provided with -n flag" << endl;
        cout << "Use \"-n <k>\", \n k should be as large as the number of nodes used." << endl;
        return 0;
      }
      n = std::stoi(args[i + 1]);
      if (n < 0) {
        cout << "No valid number of nodes provided with -n flag" << endl;
        cout << "Use \"-n <k>\", k integer > 0" << endl;
        return 0;
      }
    } else if (args[i] == "-r") { // number of random networks
      string number = args[i + 1];
      if (any_of(number.begin(), number.end(), [](char c) { return !isdigit(c); })) {
        cout << "No valid number of nodes provided with -r flag" << endl;
        cout << "Use \"-r <positive_integer>\"" << endl;
        return 0;
      }
      randomNetworks = std::stoi(args[i + 1]);
    } else if (args[i] == "-e") { // number of random shuffles
      string number = args[i + 1];
      if (any_of(number.begin(), number.end(), [](char c) { return !isdigit(c); })) {
        cout << "No valid number of nodes provided with -e flag" << endl;
        cout << "Use \"-e <positive_integer>\"" << endl;
        return 0;
      }
      randomShuffles = std::stoi(args[i + 1]);
    } else {
      cout << "Invalid argument flag provided." << endl;
      cout << "More information in the readme provided" << endl;
      return 0;
    }
  }
  // -r is optional, but should be bigger than 0
  if (randomNetworks < 0) {
    cout << "No valid number of random networks provided with -r flag" << endl;
    cout << "Use \"-r <k>\", k integer > 0" << endl;
    return 0;
  }
  
  // -e is optional, but should be bigger than 0
  if (randomShuffles < 0) {
    cout << "No valid number of edge swaps provided with -e flag" << endl;
    cout << "Use \"-e <k>\", k integer > 0" << endl;
    return 0;
  }
  
  // -s **size**, is mandatory
  if (motifSize == -1 || motifSize < 3 || motifSize > 4) {
    cout << "No valid motifSize provided with -s flag" << endl;
    cout << "Use \"-s <k>\", 3 <= k <= 4" << endl;
    return 0; 
  } 
  
  // -i **inputFile**, is mandatory
  if (inputFile.empty()) {
    cout << "No input file provided with -i flag" << endl;
    return 0; 
  }
  
  if (motifSize != 3 && algorithm == 3) {
    cout << "Sorry, the triangle method is only available for Size = 3" << endl;
    cout << "Please choose another algorithm using -a <k>, 1 <= k <= 4" << endl;
    return 0; 
  }
  
  std::ofstream fout;
  std::ostream& out = fout;
  
  
  fout.open(outputFile);
  if (fout.fail()) {
    cout << "output file: " << outputFile << " - could not be opened" << endl;
    exit(EXIT_FAILURE);
  };
  // https://oeis.org/A323817
  IsomorphismHyper::precalc(3);
  IsomorphismHyper::precalc(4);
  Hypergraph h(n); 
  cout << "-----------------------------------------------" << endl;
  cout << "Reading hypergraph from file: " << inputFile << endl;
  h.readFromFile(inputFile);
  cout << "Hypergraph read from file: " << inputFile << endl;
  cout << "-----------------------------------------------" << endl;
  cout << "Number of nodes: " << h.getNodeCount() << endl;
  cout << "Number of hyperedges: " << h.getEdgeCount() << endl;
  cout << "-----------------------------------------------" << endl;
  vector<int> size = h.getEdgeBySize();
  for (int i = 2; i <= MAX_EDGE_SIZE; i++) {
    cout << "E" << i << ": " << size[i] << endl;
  }
  
  if (algorithm == -1) algorithm = motifSize; // Set default algorithm - If size 3 => Triangle, size 4 => FaSE
  
  cout << "-----------------------------------------------" << endl;
  switch (algorithm) {
    case 1: cout << "Method: ESU Baseline" << endl; out << "Method: ESU Baseline" << endl; break;
    case 2: cout << "Method: ESU Modified" << endl; out << "Method: ESU Modified" << endl; break;
    case 3: cout << "Method: Triangle" << endl; out << "Method: Triangle" << endl; break;
    case 4: cout << "Method: FaSE" << endl; out << "Method: FaSE" << endl; break;
  }
  
  if (!motif) { // network-census of whole network
    CountingMethods::networkCensus(h, motifSize, detailedOutput, algorithm, out);
  } else { // find motifs - significance profile of each subgraph
    CountingMethods::findMotifs(h, motifSize, detailedOutput, significanceProfile, randomNetworks, randomShuffles, algorithm, out);
  }
  
  cout << "Results saved to file: " << outputFile << '\n';
  fout.close();
  
  return 0; 
}







