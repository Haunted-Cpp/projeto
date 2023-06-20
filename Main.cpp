#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <stack>
#include <map>
#include <memory>
#include <chrono>
#include <cmath>
#include <queue>

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
  
  
  int motifSize = -1;
  string inputFile, outputFile, task;
  bool detailedOutput = false;
  
  vector<string> args;
  for (int i = 1; i < argc; i++) {
    string arg = string(argv[i]);
    if (arg == "-d") detailedOutput = true;
    else args.emplace_back(arg);
  }
  
  if ( args.size() & 1 ) { // If number of arguments is odd, then a pair is missing
    cout << "Invalid number of arguments provided." << '\n';
    cout << "Each flag must be followed by an argument" << '\n';
    cout << "Use the format: -flag argument" << '\n';
    cout << "More information provide in the readme file" << '\n';
    return 0;
  }
  
  for (int i = 0; i < (int) args.size() - 1; i += 2) {
    if (args[i] == "-s") {
      string number = args[i + 1];
      if (any_of(number.begin(), number.end(), [](char c) { return !isdigit(c); })) {
        cout << "No valid motifSize provided with -s flag" << '\n';
        cout << "Use \"-s k\", where 3 <= k <= 4" << '\n';
        return 0;
      }
      motifSize = std::stoi(args[i + 1]);
    } else if (args[i] == "-i") {
      inputFile = args[i + 1];
    } else if (args[i] == "-o") {
      outputFile = args[i + 1];
    } else if (args[i] == "-m") {
      task = args[i + 1];
    } else if (args[i] == "-d") {
      detailedOutput = true;
    } else {
      cout << "Invalid argument flag provided." << '\n';
      cout << "More information in the readme provided" << '\n';
      return 0;
    }
  }
  
  // -s **size**, is mandatory
  if (motifSize == -1 || motifSize < 3 || motifSize > 4) {
    cout << "No valid motifSize provided with -s flag" << '\n';
    cout << "Use \"-s k\", where 3 <= k <= 4" << '\n';
    return 0; 
  } 
  
  // -i **inputFile**, is mandatory
  if (inputFile.empty()) {
    cout << "No input file provided with -i flag" << '\n';
    return 0; 
  }
  
  // -m **method**, is mandatory
  if (task != "count" && task != "motif") {
    cout << "No valid method provided with -m flag." << '\n';
    cout << "Valid options are: \"count\", \"motif\"" << '\n'; 
    return 0; 
  }
  
  //cout << motifSize << '\n';
  //cout << inputFile << '\n';
  //cout << outputFile << '\n';
  //cout << task << '\n';
  
  Hypergraph h; 
  h.readFromFile(inputFile);
  
  //cout << detailedOutput << '\n';
  if (task == "count") { // network-census of whole network
    ESU::networkCensus(h, motifSize, outputFile, detailedOutput);
  } else { // find significance profile of each subgraph
    ESU::findMotifs(h, motifSize, outputFile, detailedOutput);
  }
  
  return 0; 
}








