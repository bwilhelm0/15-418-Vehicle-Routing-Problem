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
#define GRANULARITY 128


struct VRP {
  vector<int> list;
  int numVehicles;

  // Can you compare vectors with ==?
  friend bool operator==(const VRP& a, const VRP& b) {
    return a.list == b.list && a.numVehicles == b.numVehicles;
  }
};

struct VRPsolution {
    vector<vector<int>> routes;
    int cost;
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

void printMatrix(int **adjacencyMatrix, int size) {
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size - 1; col++) {
            cout << adjacencyMatrix[row][col] << ", ";
        }
        cout << adjacencyMatrix[row][size - 1] << endl;
    }
}

int** read_uber_data()
{

    int size = 64;//number of nodes 
    string fname;
    vector <string> row;
    string line, word, temp;
    vector<vector<string>> content;
    fname = "uber_condensed.csv";
    fstream file (fname, ios::in);

    if (file.is_open()){
        while(getline(file,line)){
            row.clear();
            stringstream str(line);
            while(getline(str,word,',')){       //CSV are comma deliminated 
                row.push_back(word);
                content.push_back(row);
            }
        }
    } else cout << "Could not open file" << endl;


    int** matrix = new int*[size];
    // for (int i = 0; i < size; ++i) {
    //     matrix[i] = new int[size];
    //     for (int j = 0; j < size; j++) {
    //         matrix[i][j] = adjacencyMatrix[i][j];
    //     }
    // }


    // for (int i=0;(size_t) i<content.size();i++){
    //     int node1 = stoi(content[i][1]);
    //     int node2 = stoi(content[i][2]);

    //     matrix[node1][node2] = stoi(content[i][3]);
    // }

    // printMatrix(matrix,N);
    return matrix;
}  


class Node
{
public:
    vector<pair<int, int>> path;
    int **matrix_reduced;
    int cost;
    int vertex;
    int level;
};

Node* newNode(int **matrix_parent, int size, vector<pair<int, int>> const &path, int level, int i, int j)
{
    Node* node = new Node;
    node->path = path;
    if (level != 0)
        node->path.push_back(make_pair(i, j));

    node->matrix_reduced = new int*[size];
    for (int i = 0; i < size; ++i) {
        node->matrix_reduced[i] = new int[size];
        for (int j = 0; j < size; j++) {
            node->matrix_reduced[i][j] = matrix_parent[i][j];
        }
    }

    for (int k = 0; level != 0 && k < size; k++)
    {
        node->matrix_reduced[i][k] = INF;
        node->matrix_reduced[k][j] = INF;
    }

    node->matrix_reduced[j][0] = INF;
    node->level = level;
    node->vertex = j;
    return node;
}

void reduce(int**matrix_reduced, int *row, int *col, int size, vector<int> nodes)
{
    fill_n(row, size, INF);
    fill_n(col, size, INF);

    for (auto &i : nodes)
        for (auto &j : nodes)
            if (matrix_reduced[i][j] < row[i])
                row[i] = matrix_reduced[i][j];

    for (auto &i : nodes)
        for (auto &j : nodes)
            if (matrix_reduced[i][j] != INF && row[i] != INF)
                matrix_reduced[i][j] -= row[i];

    for (auto &i : nodes)
        for (auto &j : nodes)
            if (matrix_reduced[i][j] < col[j])
                col[j] = matrix_reduced[i][j];

    for (auto &i : nodes)
        for (auto &j : nodes)
            if (matrix_reduced[i][j] != INF && col[j] != INF)
                matrix_reduced[i][j] -= col[j];

    // for (int i = 0; i < size; i++)
    //     for (int j = 0; j < size; j++)
    //         if (matrix_reduced[i][j] < row[i])
    //             row[i] = matrix_reduced[i][j];

    // for (int i = 0; i < size; i++)
    //     for (int j = 0; j < size; j++)
    //         if (matrix_reduced[i][j] != INF && row[i] != INF)
    //             matrix_reduced[i][j] -= row[i];

    // for (int i = 0; i < size; i++)
    //     for (int j = 0; j < size; j++)
    //         if (matrix_reduced[i][j] < col[j])
    //             col[j] = matrix_reduced[i][j];

    // for (int i = 0; i < size; i++)
    //     for (int j = 0; j < size; j++)
    //         if (matrix_reduced[i][j] != INF && col[j] != INF)
    //             matrix_reduced[i][j] -= col[j];
}

int calcCost(int **matrix_reduced, int size, vector<int> nodes)
{
    int row[size];
    int col[size];

    reduce(matrix_reduced, row, col, size, nodes);

    int cost = 0;
    for (auto &i : nodes) {
        cost += (row[i] != INF) ? row[i] : 0;
        cost += (col[i] != INF) ? col[i] : 0;
    }

    return cost;
}

vector<int> printPath(vector<pair<int, int>> const &list) {
    vector<int> res;
    for (int i = 0; (size_t) i < list.size() - 1; i++) {
        //cout << list[i].first << " -> ";
        res.push_back(list[i].first);
    }
    //cout << list[list.size() - 1].first << " -> " << list[list.size() - 1].second << endl;
    res.push_back(list[list.size() - 1].first);
    res.push_back(list[list.size() - 1].second);
    return res;
}

void printRoutes(vector<vector<int>> routes) {
    int numRoutes = (int) routes.size();
    for (int route = 0; route < numRoutes; route++) {
        cout << "Vehicle: " << route << ", ";
        int nodes = (int) routes[route].size();
        for (int node = 0; node < nodes - 1; node++) {
            cout << routes[route][node] << " -> ";
        }
        cout << routes[route][nodes - 1] << endl;
    }
}


class comp {
public:
    bool operator()(const Node* l, const Node* r) const
    {
        return l->cost > r->cost;
    }
};

VRPsolution tspSolve(int **adjacencyMatrix, int size, vector<int> nodes)
{
    int nodeSize = nodes.size();
    priority_queue<Node*, vector<Node*>, comp> pq;
    vector<pair<int,int>> v;
    Node* root = newNode(adjacencyMatrix, size, v, 0, -1, 0);
    root->cost = calcCost(root->matrix_reduced, size, nodes);
    pq.push(root);
    while (!pq.empty())
    {
        Node* min = pq.top();
        pq.pop();
        int i = min->vertex;
        if (min->level == nodeSize - 1)
        {
            min->path.push_back(make_pair(i, 0));
            vector<int> finalPath = printPath(min->path);
            vector<vector<int>> wrapped;
            wrapped.push_back(finalPath);
            VRPsolution done;
            done.routes = wrapped;
            done.cost = min->cost;
            return done;
        }

        for (auto &j : nodes)
        {
            if (min->matrix_reduced[i][j] != INF)
            {
                Node* child = newNode(min->matrix_reduced, size, min->path,
                    min->level + 1, i, j);

                child->cost = min->cost + min->matrix_reduced[i][j]
                            + calcCost(child->matrix_reduced, size, nodes);

                pq.push(child);
            }
        }

        delete min;
    }
}

vector<VRP> genWork(int N, int master, int proc, int pid) {
    vector<int> list;
    list.push_back(master);

    vector<vector<int>> subgraphs;
    subgraphs.push_back(list);

    for (int i = 0; i < N; i++) {
        if (i != master) {
            vector<vector<int>> tempgraphs;
            for (auto &subgraph : subgraphs) {
                tempgraphs.push_back(subgraph);
                subgraph.push_back(i);
                tempgraphs.push_back(subgraph);
            }
            subgraphs = tempgraphs;
        }
    }

    vector<VRP> work;
    for (auto &sg : subgraphs) {
        VRP subprob = {};
        subprob.list = sg;
        subprob.numVehicles = 1;
        //cout << sg.size() << endl;
        if (std::hash<VRP>{}(subprob) % proc == (size_t) pid) {
            work.push_back(subprob);
        }
    }

    return work;
}

int getGranularity(VRP &prob) {
    return pow(2, prob.list.size()) * prob.numVehicles;
}

//recursive function for solving top down
VRPsolution subsetSolve(VRP &prob, unordered_map<VRP, VRPsolution> &solnMap, int **matrix, int size, int master) { //, int proc, int pid) {
    // TURN THIS INTO A LOOP WITH A VECTOR AS A QUEUE
    VRPsolution res;
    unordered_map<VRP, VRPsolution>::const_iterator got = solnMap.find(prob);

    //check for requests
    
    //if solution already in hashtable, find and return VRPsolution
    if (got != solnMap.end()) return got->second;
    else if (prob.numVehicles == 1) {
        res = tspSolve(matrix,size,prob.list);
        solnMap.insert({prob, res});
        return res;
    } 
    // else if (getGranularity(prob) > GRANULARITY && (hash<VRP>{}(prob) % proc != pid)) {
    //     communicate
    //     add to hash table
    //     return value
    // }

    vector<pair<vector<int>, vector<int>>> subsets;
    for (int i = 2; i < pow(2, prob.list.size() - 1); i+=2) { // check - 1
        vector<int> left;
        vector<int> right;

		for (int j = 0; (size_t)j < prob.list.size(); j++) {
            if (j == master) {
                left.push_back(master);
                right.push_back(master);
            } else {
                if ((i & (1 << j)) == 0) {
                    left.push_back(prob.list[j]);
                } else {
                    right.push_back(prob.list[j]);
                }
            }
		}
        subsets.push_back(make_pair(left, right));
	}

    VRPsolution minCost;
    minCost.cost = INF;
    for (auto &p : subsets) {
        VRP next1;
        VRP next2;
        next1.list = p.first;
        next2.list = p.second;

        int maxlen = max(next1.list.size() - 1, next2.list.size() - 1);
        int minlen = min((int) next1.list.size() - 1, (int) next2.list.size() - 1);
        int fulllen = next1.list.size() + next2.list.size() - 2;

        int begin = (prob.numVehicles - maxlen > 1) ? (prob.numVehicles - maxlen) : 1;
        //int iters = min(prob.numVehicles - 1, min((int) next1.list.size() - 1, (int) next2.list.size() - 1));
        int iters = (prob.numVehicles - 1 <= minlen) ? prob.numVehicles - 1 : fulllen - (prob.numVehicles - 1);
        //cout << "Begin and Iteration:  " << begin << "," << iters << endl;
        for (int vehicles = begin; vehicles < iters + begin; vehicles++) {
            next1.numVehicles = (next1.list.size() < next2.list.size()) ? vehicles : prob.numVehicles - vehicles;
            next2.numVehicles = (next1.list.size() < next2.list.size()) ? prob.numVehicles - vehicles : vehicles;

            vector<vector<int>> doub;
            doub.push_back(next1.list);
            doub.push_back(next2.list);

            VRPsolution soln1;
            VRPsolution soln2;

            // got = solnMap.find(next1);
            // //if solution already in hashtable, find and return VRPsolution
            // if (got != solnMap.end()) return got->second;
            // else 
            soln1 = subsetSolve(next1, solnMap, matrix, size, master);
            
            // got = solnMap.find(next2);
            // //if solution already in hashtable, find and return VRPsolution
            // if (got != solnMap.end()) return got->second;
            // else 
            soln2 = subsetSolve(next2, solnMap, matrix, size, master);

            if (minCost.cost > soln1.cost + soln2.cost) {
                soln1.routes.insert(soln1.routes.end(), soln2.routes.begin(), soln2.routes.end());
                minCost.routes = soln1.routes;
                minCost.cost = soln1.cost + soln2.cost;
            }
            if (prob.numVehicles == 4 && (next1.numVehicles > next1.list.size() - 1 || next2.numVehicles > next2.list.size() - 1)) {
                cout << "Next1 Vehicles: " << next1.numVehicles << endl;
                cout << "Next2 Vehicles: " << next2.numVehicles << endl;
                cout << (soln1.cost + soln2.cost) << endl;
                printRoutes(doub);
            }
        }
    }
    return minCost;
}



int main(int argc, char *argv[])
{
    int pid;
    int proc;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    // Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    MPI_Barrier(MPI_COMM_WORLD);
    //Timer totalTime;
    //cout << pid << endl;

    int size = 8;
    int vehicles = 5;
    int master = 0;
    int N = 13;

    int adjacencyMatrix[N][N] =
    // {
    //     {INF, 2},
    //     {2, INF}
    // };
    // {
    //     { INF, 20,  30,  10,  11},
    //     { 15,  INF, 16,  4,   2 },
    //     { 3,   5,   INF, 2,   4 },
    //     { 19,  6,   18,  INF, 3 },
    //     { 16,  4,   7,   16,  INF}
    // };
    {
      {INF, 2451, 713, 1018, 1631, 1374, 2408, 213, 2571, 875, 1420, 2145, 1972},
      {2451, INF, 1745, 1524, 831, 1240, 959, 2596, 403, 1589, 1374, 357, 579},
      {713, 1745, INF, 355, 920, 803, 1737, 851, 1858, 262, 940, 1453, 1260},
      {1018, 1524, 355, INF, 700, 862, 1395, 1123, 1584, 466, 1056, 1280, 987},
      {1631, 831, 920, 700, INF, 663, 1021, 1769, 949, 796, 879, 586, 371},
      {1374, 1240, 803, 862, 663, INF, 1681, 1551, 1765, 547, 225, 887, 999},
      {2408, 959, 1737, 1395, 1021, 1681, INF, 2493, 678, 1724, 1891, 1114, 701},
      {213, 2596, 851, 1123, 1769, 1551, 2493, INF, 2699, 1038, 1605, 2300, 2099},
      {2571, 403, 1858, 1584, 949, 1765, 678, 2699, INF, 1744, 1645, 653, 600},
      {875, 1589, 262, 466, 796, 547, 1724, 1038, 1744, INF, 679, 1272, 1162},
      {1420, 1374, 940, 1056, 879, 225, 1891, 1605, 1645, 679, INF, 1017, 1200},
      {2145, 357, 1453, 1280, 586, 887, 1114, 2300, 653, 1272, 1017, INF, 504},
      {1972, 579, 1260, 987, 371, 999, 701, 2099, 600, 1162, 1200, 504, INF},
    };

    int** matrix = new int*[size];
    for (int i = 0; i < size; ++i) {
        matrix[i] = new int[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = adjacencyMatrix[i][j];
        }
    }


    //read_uber_data(matrix);

    //printMatrix(matrix, size);

    VRP prob;
    vector<int> allPoints;
    for(int i = 0; i < size; i++) {
        allPoints.push_back(i);
    }
    prob.list = allPoints;
    prob.numVehicles = vehicles;
    
    unordered_map<VRP, VRPsolution> routeTable;

    // Pre-compute some values
    // vector<VRP> work = genWork(size, master, proc, pid);
    // for (auto &subProb : work) {
    //     VRPsolution solved = tspSolve(matrix, size, subProb.list);
    //     routeTable.insert({subProb, solved});
    // }

    VRPsolution finished = subsetSolve(prob, routeTable, matrix, size, master);
    cout << "Cost is " << finished.cost << endl;
    printRoutes(finished.routes);
    

    
    // Ring reduce to update values

    // Use hash table and branching algo to find cheapest route
    // The key is (the set of points to visit, num vehicles), and the value is (ordered visits for each vehicle, the cost of the trip)


    for (int i = 0; i < size; ++i)
        delete [] matrix[i];
    delete [] matrix;
    
    MPI_Finalize();
}
