using namespace std;
#include "VRP.h"

void printMatrix(int **adjacencyMatrix, int size) {
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size - 1; col++) {
            cout << adjacencyMatrix[row][col] << ", ";
        }
        cout << adjacencyMatrix[row][size - 1] << endl;
    }
}

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

vector<int> printPath(vector<pair<int, int>> const &list, bool printOn) {
    vector<int> res;
    for (int i = 0; (size_t) i < list.size() - 1; i++) {
        if(printOn) cout << list[i].first << " -> ";
        res.push_back(list[i].first);
    }
    if(printOn) cout << list[list.size() - 1].first << " -> " << list[list.size() - 1].second << endl;
    res.push_back(list[list.size() - 1].first);
    res.push_back(list[list.size() - 1].second);
    return res;
}


VRPsolution tspSolve(int **adjacencyMatrix, int size, vector<int> nodes, bool printOn)
{
    int nodeSize = nodes.size();
    priority_queue<Node*, vector<Node*>, comp> pq;
    vector<pair<int,int>> v;
    Node* root = newNode(adjacencyMatrix, size, v, 0, -1, 0);
    root->time = calcCost(root->matrix_reduced, size, nodes);
    pq.push(root);
    while (!pq.empty())
    {
        Node* min = pq.top();
        pq.pop();
        int i = min->vertex;
        if (min->level == nodeSize - 1)
        {
            min->path.push_back(make_pair(i, 0));
            vector<int> finalPath = printPath(min->path, printOn);
            vector<vector<int>> wrapped;
            wrapped.push_back(finalPath);
            VRPsolution done;
            done.routes = wrapped;
            done.time = min->time;
            return done;
        }

        for (auto &j : nodes)
        {
            if (min->matrix_reduced[i][j] != INF)
            {
                Node* child = newNode(min->matrix_reduced, size, min->path,
                    min->level + 1, i, j);

                child->time = min->time + min->matrix_reduced[i][j]
                            + calcCost(child->matrix_reduced, size, nodes);

                pq.push(child);
            }
        }

        delete min;
    }
}
