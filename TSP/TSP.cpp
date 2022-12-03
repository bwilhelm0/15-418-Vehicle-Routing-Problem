using namespace std;
#include <queue>
#include <unordered_map>
#include <set>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>

#define INF numeric_limits<int>::max()



struct VRP {
  vector<int> list;
  int numVehicles;

  friend bool operator==(const VRP& a, const VRP& b) {
    return a.list == b.list && a.numVehicles == b.numVehicles;
  }
};

    std::hash<std::vector<int>> vecHash;
// Our custom std::hash specialization for Planet
template <>
struct std::hash<VRP> {
  size_t operator()(const VRP& p) const noexcept {
    return std::hash<int>{}(p.numVehicles) + vecHash(p);
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

void read_uber_data(int** matrix)
{

  int N = 80; //number of nodes 
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

  for (int i=0;i<content.size();i++){
    int node1 = stoi(content[i][1]);
    int node2 = stoi(content[i][2]);

    matrix[node1][node2] = stoi(content[i][3]);
  }

  printMatrix(matrix,N);
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

void reduce(int**matrix_reduced, int *row, int *col, int size)
{
    fill_n(row, size, INF);
    fill_n(col, size, INF);

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (matrix_reduced[i][j] < row[i])
                row[i] = matrix_reduced[i][j];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (matrix_reduced[i][j] != INF && row[i] != INF)
                matrix_reduced[i][j] -= row[i];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (matrix_reduced[i][j] < col[j])
                col[j] = matrix_reduced[i][j];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (matrix_reduced[i][j] != INF && col[j] != INF)
                matrix_reduced[i][j] -= col[j];
}

int calcCost(int **matrix_reduced, int size)
{
    int row[size];
    int col[size];

    reduce(matrix_reduced, row, col, size);

    int cost = 0;
    for (int i = 0; i < size; i++) {
        cost += (row[i] != INF) ? row[i] : 0;
        cost += (col[i] != INF) ? col[i] : 0;
    }

    return cost;
}

void printPath(vector<pair<int, int>> const &list) {
    for (int i = 0; i < list.size() - 1; i++)
        cout << list[i].first << " -> ";
    cout << list[list.size() - 1].first << " -> " << list[list.size() - 1].second << endl;
}



class comp {
public:
    bool operator()(const Node* l, const Node* r) const
    {
        return l->cost > r->cost;
    }
};

int solve(int **adjacencyMatrix, int size)
{
    priority_queue<Node*, vector<Node*>, comp> pq;
    vector<pair<int, int>> v;
    Node* root = newNode(adjacencyMatrix, size, v, 0, -1, 0);
    root->cost = calcCost(root->matrix_reduced, size);
    pq.push(root);
    while (!pq.empty())
    {
        Node* min = pq.top();
        pq.pop();
        int i = min->vertex;
        if (min->level == size - 1)
        {
            min->path.push_back(make_pair(i, 0));
            //printPath(min->path);
            return min->cost;
        }

        for (int j = 0; j < size; j++)
        {
            if (min->matrix_reduced[i][j] != INF)
            {
                Node* child = newNode(min->matrix_reduced, size, min->path,
                    min->level + 1, i, j);

                child->cost = min->cost + min->matrix_reduced[i][j]
                            + calcCost(child->matrix_reduced, size);

                pq.push(child);
            }
        }

        delete min;
    }
}

int retrieve_cost(vector<int> graph, int num_vehicles, int master) {
    // hash table with graph and number of vehicles as key.
    // Solve in branching method where you split a number of vehicles to a certain tsp.

    return 0;
}

vector<vector<vector<int>>> genWork(int N, int master, int proc) {
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

    int i = 0;
    vector<vector<vector<int>>> work;
    work.resize(proc);
    for (auto &sg : subgraphs) {
        work[i % proc].push_back(sg);
        i++;
    }

    return work;
}



int main()
{
    int size = 13;
    int N = 13;
    int proc = 4;

    

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

    read_uber_data(matrix);

    //printMatrix(matrix, size);
    int res = solve(matrix, size);
    //printf("Cost is %d\n", res);

    // Vehicle Routing Algo

    // Creates 2^(N - 1) subsets divided among proc
    vector<vector<vector<int>>> work = genWork(size, 0, proc);

    // Distribute to all processes in ring

    // Process each subset

    // Format vectors in order given but points in order of traversal with cost at end of list

    // Ring reduce to update values

    // Use hash table and branching algo to find cheapest route
    // The key is (the set of points to visit, num vehicles), and the value is (ordered visits for each vehicle, the cost of the trip)
    //unordered_map<VRP, pair<vector<vector<int>>, int>> routeTable;

    for (int i = 0; i < size; ++i)
        delete [] matrix[i];
    delete [] matrix;

    return 0;
}
