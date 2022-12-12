#ifndef VRP_H
#define VRP_H

using namespace std;
#include <queue>
#include <unordered_map>
#include <set>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>

#include "timing.h"
#include "mpi.h"
#include "omp.h"

#define INF numeric_limits<int>::max()

struct StartupOptions {
  int nodes;
  int vehicles;
  int seed;
  bool printPaths;
  std::string inputFile;
};

inline StartupOptions parseOptions(int argc, char *argv[]) {
  StartupOptions rs;
  rs.printPaths = false;
  rs.seed = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0)
      rs.nodes = atoi(argv[i + 1]);
    else if (strcmp(argv[i], "-v") == 0)
      rs.vehicles = atoi(argv[i + 1]);
    else if (strcmp(argv[i], "-p") == 0)
      rs.printPaths = true;
    else if (strcmp(argv[i], "-s") == 0)
      rs.seed = atoi(argv[i + 1]);
  }
  return rs;
}


// TSP DEFINITIONS
class Node
{
public:
    vector<pair<int, int>> path;
    int **matrix_reduced;
    int time;
    int vertex;
    int level;
};

class comp {
public:
    bool operator()(const Node* l, const Node* r) const
    {
        return l->time > r->time;
    }
};




// VRP DEFINITIONS
struct VRP {
  vector<int> list;
  int numVehicles;
  int master;
  int proc;
  int pid;
  int originalP;
  int originalV;
  int printPaths;

  // Can you compare vectors with ==?
  friend bool operator==(const VRP& a, const VRP& b) {
    return a.list == b.list && a.numVehicles == b.numVehicles;
  }
};

struct VRPsolution {
    vector<vector<int>> routes;
    int time;
};

// Our custom std::hash specialization for VRP
template <>
struct std::hash<VRP> {
  size_t operator()(const VRP& p) const noexcept {
    int sum = 0;
    for (auto &n : p.list) {
        sum += std::hash<int>{}(n);
    }
    return std::hash<int>{}(p.numVehicles) + sum;
  }
};

enum MSG_TAG {LEVEL_SYNC = 0, REQUEST, ANSWER};


// TSP Function
VRPsolution tspSolve(int **adjacencyMatrix, int size, vector<int> nodes, bool printOn);

#endif
