#pragma once

#define MAX_EDGE_SIZE 4
#define MAX_HYPER_MOTIF_SIZE 4 // Should be SMALLER than MAX_EDGE_SIZE [CHANGE TO 4]
#define MAX_MOTIF_SIZE 4 
#define MAX_INPUT_N 10'000'005

using std::string;
using std::vector;
using std::pair;
using std::map;
using std::vector;
using std::cin;
using std::cout;
using std::flush;
using std::endl;
using std::ios;
using std::ofstream;
using std::istream;
using std::ostream;
using std::ifstream;
using std::max;
using std::queue;
using std::swap;
using std::tuple;
using std::make_pair;

using namespace std::chrono;

// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933#72073933
struct HashFunction {
  std::size_t operator()(std::vector<int> const& vec) const {
    std::size_t seed = vec.size();
    for(auto xx : vec) {
      unsigned int x = static_cast<unsigned int>(xx);
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = (x >> 16) ^ x;
      seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
  std::size_t operator()(std::vector< vector<int> > const& vec) const {
    std::size_t seed = vec.size();
    for(auto a : vec) {
      for (auto xx : a) {
        unsigned int x = static_cast<unsigned int>(xx);
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = (x >> 16) ^ x;
        seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
    }
    return seed;
  }
};
