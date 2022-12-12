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

#include "timing.h"
#include "mpi.h"
#include "omp.h"

#define INF numeric_limits<int>::max()
#define GRANULARITY 1

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


// TSP Function
VRPsolution tspSolve(int **adjacencyMatrix, int size, vector<int> nodes, bool printOn);

#endif
