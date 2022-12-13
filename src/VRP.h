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

// Int max to be used as infinity
#define INF numeric_limits<int>::max()

// Struct for storing startup options based on commandline input
struct StartupOptions {
  int nodes;
  int vehicles;
  int seed;
  bool printPaths;
  bool reduceComm;
  bool timeLevels;
  bool printMatrix;
  bool ringReduce;
  std::string inputFile;
};

// Function to parse commandline options into corresponding struct
inline StartupOptions parseOptions(int argc, char *argv[]) {
  StartupOptions rs;
  rs.printPaths = false;
  rs.seed = -1;
  rs.reduceComm = true;
  rs.timeLevels = false;
  rs.printMatrix = false;
  rs.ringReduce = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0)
      rs.nodes = atoi(argv[i + 1]);
    else if (strcmp(argv[i], "-v") == 0)
      rs.vehicles = atoi(argv[i + 1]);
    else if (strcmp(argv[i], "-p") == 0)
      rs.printPaths = true;
    else if (strcmp(argv[i], "-s") == 0)
      rs.seed = atoi(argv[i + 1]);
    else if (strcmp(argv[i], "-l") == 0)
      rs.reduceComm = false;
    else if (strcmp(argv[i], "-t") == 0)
      rs.timeLevels = true;
    else if (strcmp(argv[i], "-M") == 0)
      rs.printMatrix = true;
    else if (strcmp(argv[i], "-RR") == 0)
      rs.ringReduce = true;
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
struct VRPspecs {
  int master;
  int proc;
  int pid;
  int originalP;
  int originalV;
  int printPaths;
  bool reduceComm;
  bool ringReduce;
};

struct VRP {
  vector<int> list;
  int numVehicles;
  
  // Comparison function for hash table lookup
  friend bool operator==(const VRP& a, const VRP& b) {
    return a.list == b.list && a.numVehicles == b.numVehicles;
  }
};

// The solution struct of a VRP containing routes and cost
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

// Enum to distinguish message types, REDUCE must be the last tag in the enum
enum MSG_TAG {LEVEL_SYNC = 0, 
              REQUEST = 1, 
              REQUEST_COST = 2, 
              REQUEST_ROUTES = 3, 
              ANSWER = 4, 
              ANSWER_COST = 5, 
              ANSWER_ROUTES = 6, 
              REDUCE = 7,
              REDUCE_DATA = 8};

// TSP Function to be called at base level in VRP
void printMatrix(int **adjacencyMatrix, int size);
VRPsolution tspSolve(int **adjacencyMatrix, int size, vector<int> nodes, bool printOn);

#endif
