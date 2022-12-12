// using namespace std;
// #include "VRP.h"

// void printRoutes(vector<vector<int>> routes) {
//     int numRoutes = (int) routes.size();
//     for (int route = 0; route < numRoutes; route++) {
//         cout << "Vehicle: " << route << ", ";
//         int nodes = (int) routes[route].size();
//         for (int node = 0; node < nodes - 1; node++) {
//             cout << routes[route][node] << " -> ";
//         }
//         cout << routes[route][nodes - 1] << endl;
//     }
// }

// int** read_uber_data()
// {

//     int size = 64;//number of nodes
//     string fname;
//     vector <string> row;
//     string line, word, temp;
//     vector<vector<string>> content;
//     fname = "uber_condensed.csv";
//     fstream file (fname, ios::in);

//     if (file.is_open()){
//         while(getline(file,line)){
//             row.clear();
//             stringstream str(line);
//             while(getline(str,word,',')){       //CSV are comma deliminated
//                 row.push_back(word);
//                 content.push_back(row);
//             }
//         }
//     } else cout << "Could not open file" << endl;


//     int** matrix = new int*[size];
//     // for (int i = 0; i < size; ++i) {
//     //     matrix[i] = new int[size];
//     //     for (int j = 0; j < size; j++) {
//     //         matrix[i][j] = adjacencyMatrix[i][j];
//     //     }
//     // }


//     // for (int i=0;(size_t) i<content.size();i++){
//     //     int node1 = stoi(content[i][1]);
//     //     int node2 = stoi(content[i][2]);

//     //     matrix[node1][node2] = stoi(content[i][3]);
//     // }

//     // printMatrix(matrix,N);
//     return matrix;
// }

// vector<vector<int>> genWork(int N, int master) {
//     vector<int> list;
//     list.push_back(master);

//     vector<vector<int>> subgraphs;
//     subgraphs.push_back(list);

//     for (int i = 0; i < N; i++) {
//         if (i != master) {
//             vector<vector<int>> tempgraphs;
//             for (auto &subgraph : subgraphs) {
//                 tempgraphs.push_back(subgraph);
//                 subgraph.push_back(i);
//                 tempgraphs.push_back(subgraph);
//             }
//             subgraphs = tempgraphs;
//         }
//     }

//     return subgraphs;
// }

// int getGranularity(VRP &prob) {
//     return pow(2, prob.list.size()) * prob.numVehicles;
// }


// //checks for incoming data requests and responds
// pair<vector<vector<int>>, vector<MPI_Request>> fillRequests(VRP &prob, unordered_map<VRP, VRPsolution> &solnMap) {
//     //totalPoints + 1 because we send the number of vehicles as well
//     const int requestSize = prob.originalP + 1;

//     //In message there should be a 0 for each vehicle, this is largest possible message size
//     const int answerSize = prob.originalP + prob.originalV;

//     //check for requests
//     int *flags = new int[prob.proc];
//     for (int i = 0; i < prob.proc; i++) {
//         if (i == prob.pid) continue;
//         MPI_Iprobe(i, MPI_ANY_TAG, MPI_COMM_WORLD, &flags[i], MPI_STATUS_IGNORE);
//     }

//     vector<vector<int>> filledReqs;
//     vector<MPI_Request> filledReqStatus;
//     for (int i = 0; i < prob.proc; i++) {
//         if (i == prob.pid) continue;
//         if (flags[i] != 0) {
//             //cout << "communicate to " << i << "from " << prob.pid << endl;
//             VRP req;
//             req.list.resize(requestSize);
//             MPI_Recv((void *) req.list.data(), requestSize, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             //cout << "communicate to " << i << "from " << prob.pid << "passed " << endl;

//             req.numVehicles = req.list[requestSize - 1];
//             for (int i = requestSize - 2; i > 0; i--) {
//                 if (req.list[i] != 0) {
//                     req.list.resize(i + 1);
//                     break;
//                 }
//             }


//             unordered_map<VRP, VRPsolution>::const_iterator foundit = solnMap.find(req);

//             vector<int> answerVal;
//             answerVal.resize(answerSize);
//             MPI_Request answerReq;

//             if (foundit != solnMap.end()) {
//                 cout << "Return found" << endl;
//                 // Found case, send solution
//                 answerVal[0] = foundit->second.time;
//                 int pos = 1;
//                 for (auto &journey : foundit->second.routes) {
//                     copy(journey.begin(), journey.end(), answerVal.begin() + pos);
//                     pos += journey.size() - 1;
//                 }

//                 MPI_Isend((void *) answerVal.data(), answerSize, MPI_INT, i, 0, MPI_COMM_WORLD, &answerReq);
//             } else {
//                 // Not found case, send 0 vector
//                 MPI_Isend((void *) answerVal.data(), answerSize, MPI_INT, i, 0, MPI_COMM_WORLD, &answerReq);
//             }
//             filledReqs.push_back(answerVal);
//             filledReqStatus.push_back(answerReq);
//         }
//     }
//     return make_pair(filledReqs, filledReqStatus);
// }

// // Generates all sub problem pairs for a given VRP (Num vehicles > 1)
// vector<pair<VRP, VRP>> genSubs(VRP &prob) {
//     vector<pair<VRP, VRP>> resSubs;

//     vector<pair<vector<int>, vector<int>>> subsets;
//     for (int i = 2; i < pow(2, prob.list.size() - 1); i+=2) { // check - 1
//         vector<int> left;
//         vector<int> right;

// 		for (int j = 0; (size_t)j < prob.list.size(); j++) {
//             if (j == prob.master) {
//                 left.push_back(prob.master);
//                 right.push_back(prob.master);
//             } else {
//                 if ((i & (1 << j)) == 0) {
//                     left.push_back(prob.list[j]);
//                 } else {
//                     right.push_back(prob.list[j]);
//                 }
//             }
// 		}
//         subsets.push_back(make_pair(left, right));
// 	}

//     for (auto &p : subsets) {
//         VRP next1;
//         VRP next2;
//         next1.list = p.first;
//         next1.master = prob.master;
//         next1.proc = prob.proc;
//         next1.pid = prob.pid;
//         next1.originalP = prob.originalP;
//         next1.originalV = prob.originalV;
//         next1.printPaths = prob.printPaths;

//         next2.list = p.second;
//         next2.master = prob.master;
//         next2.proc = prob.proc;
//         next2.pid = prob.pid;
//         next2.originalP = prob.originalP;
//         next2.originalV = prob.originalV;
//         next2.printPaths = prob.printPaths;


//         int maxlen = max(next1.list.size() - 1, next2.list.size() - 1);
//         int minlen = min((int) next1.list.size() - 1, (int) next2.list.size() - 1);

//         int smallVehicles = 1;
//         int largeVehicles = prob.numVehicles - 1;

//         if (largeVehicles > maxlen) {
//             smallVehicles += largeVehicles - maxlen;
//             largeVehicles = maxlen;
//         }

//         while (smallVehicles <= minlen && largeVehicles >= 1) {
//             next1.numVehicles = (next1.list.size() < next2.list.size()) ? smallVehicles : largeVehicles;
//             next2.numVehicles = (next1.list.size() < next2.list.size()) ? largeVehicles : smallVehicles;

//             resSubs.push_back(make_pair(next1, next2));

//             smallVehicles += 1;
//             largeVehicles -= 1;
//         }
//     }

//     return resSubs;
// }


// //recursive function for solving top down
// VRPsolution subsetSolve(VRP &prob, unordered_map<VRP, VRPsolution> &solnMap, int **matrix, int size) {
//     // TURN THIS INTO A LOOP WITH A VECTOR AS A QUEUE
//     const int requestSize = prob.originalP + 1;        //totalPoints + 1 because we send the number of vehicles as well
//     const int answerSize = prob.originalP + prob.originalV;   //In message there should be a 0 for each vehicle, this is largest possible message size


//     pair<vector<vector<int>>, vector<MPI_Request>> currentlyFilling = fillRequests(prob, solnMap);
//     vector<MPI_Request> filledReqStatus = currentlyFilling.second;
//     vector<vector<int>> filledReqs = currentlyFilling.first;


//     VRPsolution res;
//     unordered_map<VRP, VRPsolution>::const_iterator got = solnMap.find(prob);

//     // if solution already in hashtable, find and return VRPsolution
//     if (got != solnMap.end()) {
//         MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(), MPI_STATUSES_IGNORE);
//         return got->second;
//     }
//     // if problem complexity is above granularity limit, and it is not this pid's, request it
//     if (getGranularity(prob) > GRANULARITY && ((int) (hash<VRP>{}(prob) % prob.proc) != prob.pid)) {
//         int hashedProc = hash<VRP>{}(prob) % prob.proc;

//         //cout << "request communication from " << hashedProc << " to " << prob.pid << endl;
//         vector<int> request;
//         request.resize(requestSize);
//         vector<int> answer;
//         answer.resize(answerSize);

//         copy(prob.list.begin(), prob.list.end(), request.begin());
//         request[requestSize-1] = prob.numVehicles;

//         // Request solution to problem
//         MPI_Send((void *) request.data(), requestSize, MPI_INT, hashedProc, 0, MPI_COMM_WORLD);

//         // Asynch Receive solution if it has been solved, keep checking for sends while waiting to avoid deadlock
//         MPI_Request recRequest;
//         MPI_Irecv((void *) answer.data(), answerSize, MPI_INT, hashedProc, 0, MPI_COMM_WORLD, &recRequest);

//         //manually check if receive request was completed
//         int recFlag = 0;
//         MPI_Request_get_status(recRequest,&recFlag,MPI_STATUS_IGNORE);

//         // Fill requests until the packet has been received
//         while (!recFlag) {
//             pair<vector<vector<int>>, vector<MPI_Request>> newFills = fillRequests(prob, solnMap);

//             for (auto &path : newFills.first) {
//                 filledReqs.push_back(path);
//             }

//             for (auto &status : newFills.second) {
//                 filledReqStatus.push_back(status);
//             }

//             MPI_Request_get_status(recRequest,&recFlag,MPI_STATUS_IGNORE);
//         }

//         // Checks if this was solved or not (first element should be cost, if -1, not solved)
//         if (answer[0] > 0) {
//             cout << "Received solved" << endl;
//             VRPsolution currSolution;
//             currSolution.time = answer[0];

//             vector<int> path;
//             int masters = 0;
//             for (int i = 1; i < answerSize; i++) {
//                 if (answer[i] == prob.master) {
//                     if (masters != 0) {
//                         path.push_back(answer[i]);
//                         currSolution.routes.push_back(path);
//                         path.resize(0);         //resize 0 because we insert
//                     }
//                     masters += 1;

//                     if (masters > prob.numVehicles) {
//                         break;
//                     }
//                 }
//                 path.push_back(answer[i]);
//             }
//             // Add solution to hash table
//             solnMap.insert({prob, currSolution});
//             MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(), MPI_STATUSES_IGNORE);
//             return currSolution;
//         }
//     }
//     // if the problem only has 1 vehicle, then calculate tsp
//     if (prob.numVehicles == 1) {
//         res = tspSolve(matrix, size, prob.list, prob.printPaths);
//         solnMap.insert({prob, res});
//         MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(), MPI_STATUSES_IGNORE);
//         return res;
//     }

//     vector<pair<VRP, VRP>> pairSubs = genSubs(prob);
//     VRPsolution minTime;
//     minTime.time = INF;
//     for (auto &pairSub : pairSubs) {
//         VRPsolution soln1;
//         VRPsolution soln2;

//         soln1 = subsetSolve(pairSub.first, solnMap, matrix, size);
//         soln2 = subsetSolve(pairSub.second, solnMap, matrix, size);

//         if (minTime.time > max(soln1.time, soln2.time)) {
//             soln1.routes.insert(soln1.routes.end(), soln2.routes.begin(), soln2.routes.end());
//             minTime.routes = soln1.routes;
//             minTime.time = max(soln1.time, soln2.time);
//         }
//     }

//     solnMap.insert({prob, minTime});
//     MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(), MPI_STATUSES_IGNORE);
//     return minTime;
// }

// // int main(int argc, char *argv[])
// // {
// //     int pid;
// //     int proc;

// //     // Initialize MPI
// //     MPI_Init(&argc, &argv);
// //     // Get process rank
// //     MPI_Comm_rank(MPI_COMM_WORLD, &pid);
// //     // Get total number of processes specificed at start of run
// //     MPI_Comm_size(MPI_COMM_WORLD, &proc);

// //     MPI_Barrier(MPI_COMM_WORLD);
// //     Timer totalTime;

// //     int size = 10;
// //     int vehicles = 4;
// //     int master = 0; // Keep at 0 until debugged
// //     int N = 13;
// //     bool printPaths = false;

// //     int adjacencyMatrix[N][N] =
// //     // {
// //     //     {INF, 2},
// //     //     {2, INF}
// //     // };
// //     // {
// //     //     { INF, 20,  30,  10,  11, 16},
// //     //     { 15,  INF, 16,  10,  10, 17},
// //     //     { 10,   10,   INF, 10, 10, 18},
// //     //     { 19,  10,   18,  INF, 10, 19},
// //     //     { 16, 14,   17,   16,  INF, 20},
// //     //     { 23, 18,   21,   15,  24, INF}
// //     // };
// //     {
// //       {INF, 2451, 713, 1018, 1631, 1374, 2408, 213, 2571, 875, 1420, 2145, 1972},
// //       {2451, INF, 1745, 1524, 831, 1240, 959, 2596, 403, 1589, 1374, 357, 579},
// //       {713, 1745, INF, 355, 920, 803, 1737, 851, 1858, 262, 940, 1453, 1260},
// //       {1018, 1524, 355, INF, 700, 862, 1395, 1123, 1584, 466, 1056, 1280, 987},
// //       {1631, 831, 920, 700, INF, 663, 1021, 1769, 949, 796, 879, 586, 371},
// //       {1374, 1240, 803, 862, 663, INF, 1681, 1551, 1765, 547, 225, 887, 999},
// //       {2408, 959, 1737, 1395, 1021, 1681, INF, 2493, 678, 1724, 1891, 1114, 701},
// //       {213, 2596, 851, 1123, 1769, 1551, 2493, INF, 2699, 1038, 1605, 2300, 2099},
// //       {2571, 403, 1858, 1584, 949, 1765, 678, 2699, INF, 1744, 1645, 653, 600},
// //       {875, 1589, 262, 466, 796, 547, 1724, 1038, 1744, INF, 679, 1272, 1162},
// //       {1420, 1374, 940, 1056, 879, 225, 1891, 1605, 1645, 679, INF, 1017, 1200},
// //       {2145, 357, 1453, 1280, 586, 887, 1114, 2300, 653, 1272, 1017, INF, 504},
// //       {1972, 579, 1260, 987, 371, 999, 701, 2099, 600, 1162, 1200, 504, INF},
// //     };

// //     int** matrix = new int*[size];
// //     for (int i = 0; i < size; ++i) {
// //         matrix[i] = new int[size];
// //         for (int j = 0; j < size; j++) {
// //             matrix[i][j] = adjacencyMatrix[i][j];
// //         }
// //     }

// //     VRP prob;
// //     vector<int> allPoints;
// //     for(int i = 0; i < size; i++) {
// //         allPoints.push_back(i);
// //     }

// //     prob.list = allPoints;
// //     prob.numVehicles = vehicles;
// //     prob.master = master;
// //     prob.proc = proc;
// //     prob.pid = pid;
// //     prob.originalP = allPoints.size();
// //     prob.originalV = vehicles;
// //     prob.printPaths = printPaths;

// //     vector<VRP> problems;
// //     unordered_map<VRP, VRPsolution> routeTable;

// //     // generate all subsets of top level

// //     vector<pair<VRP, VRP>> allSubs = genSubs(prob);

// //     for (auto &sp : allSubs) {
// //         if ((int) (hash<VRP>{}(prob) % proc) == pid) break;
// //         if ((int) (hash<VRP>{}(sp.first) % proc) == pid) {
// //             subsetSolve(sp.first, routeTable, matrix, size);
// //         }
// //         if ((int) (hash<VRP>{}(sp.second) % proc) == pid) {
// //             subsetSolve(sp.second, routeTable, matrix, size);
// //         }
// //     }

// //     // Good place for ring reduce

// //     if ((int) (hash<VRP>{}(prob) % proc) == pid) {
// //         cout << "GOING IN " << pid << endl;
// //         cout << "Table size " << routeTable.size() << endl;
// //         prob.printPaths = true;
// //         VRPsolution finished = subsetSolve(prob, routeTable, matrix, size);
// //         cout << "Table size after " << routeTable.size() << endl;
// //         cout << "Time Cost is " << finished.time << endl;
// //         printRoutes(finished.routes);
// //     } else {
// //         cout << "FINISHED " << pid << endl;
// //         while (1) {
// //             fillRequests(prob, routeTable);
// //         }
// //     }

// //     for (int i = 0; i < size; i++)
// //         delete [] matrix[i];
// //     delete [] matrix;

// //     MPI_Finalize();
// // }

// // int main(int argc, char *argv[])
// // {
// //     int pid;
// //     int proc;

// //     // Initialize MPI
// //     MPI_Init(&argc, &argv);
// //     // Get process rank
// //     MPI_Comm_rank(MPI_COMM_WORLD, &pid);
// //     // Get total number of processes specificed at start of run
// //     MPI_Comm_size(MPI_COMM_WORLD, &proc);

// //     MPI_Barrier(MPI_COMM_WORLD);
// //     Timer totalTime;

// //     int size = 10;
// //     int vehicles = 4;
// //     int master = 0; // Keep at 0 until debugged
// //     int N = 13;
// //     bool printPaths = false;

// //     int adjacencyMatrix[N][N] =
// //     // {
// //     //     {INF, 2},
// //     //     {2, INF}
// //     // };
// //     // {
// //     //     { INF, 20,  30,  10,  11, 16},
// //     //     { 15,  INF, 16,  10,  10, 17},
// //     //     { 10,   10,   INF, 10, 10, 18},
// //     //     { 19,  10,   18,  INF, 10, 19},
// //     //     { 16, 14,   17,   16,  INF, 20},
// //     //     { 23, 18,   21,   15,  24, INF}
// //     // };
// //     {
// //       {INF, 2451, 713, 1018, 1631, 1374, 2408, 213, 2571, 875, 1420, 2145, 1972},
// //       {2451, INF, 1745, 1524, 831, 1240, 959, 2596, 403, 1589, 1374, 357, 579},
// //       {713, 1745, INF, 355, 920, 803, 1737, 851, 1858, 262, 940, 1453, 1260},
// //       {1018, 1524, 355, INF, 700, 862, 1395, 1123, 1584, 466, 1056, 1280, 987},
// //       {1631, 831, 920, 700, INF, 663, 1021, 1769, 949, 796, 879, 586, 371},
// //       {1374, 1240, 803, 862, 663, INF, 1681, 1551, 1765, 547, 225, 887, 999},
// //       {2408, 959, 1737, 1395, 1021, 1681, INF, 2493, 678, 1724, 1891, 1114, 701},
// //       {213, 2596, 851, 1123, 1769, 1551, 2493, INF, 2699, 1038, 1605, 2300, 2099},
// //       {2571, 403, 1858, 1584, 949, 1765, 678, 2699, INF, 1744, 1645, 653, 600},
// //       {875, 1589, 262, 466, 796, 547, 1724, 1038, 1744, INF, 679, 1272, 1162},
// //       {1420, 1374, 940, 1056, 879, 225, 1891, 1605, 1645, 679, INF, 1017, 1200},
// //       {2145, 357, 1453, 1280, 586, 887, 1114, 2300, 653, 1272, 1017, INF, 504},
// //       {1972, 579, 1260, 987, 371, 999, 701, 2099, 600, 1162, 1200, 504, INF},
// //     };

// //     int** matrix = new int*[size];
// //     for (int i = 0; i < size; ++i) {
// //         matrix[i] = new int[size];
// //         for (int j = 0; j < size; j++) {
// //             matrix[i][j] = adjacencyMatrix[i][j];
// //         }
// //     }


// //     //verify no empty set and no full set

// //     vector<vector<int>> allSubs = genWork(size, master);
// //     unordered_map<VRP, VRPsolution> routeTable;

// //     for (int nVehicles = 1; nVehicles < vehicles; nVehicles++) {
// //         vector<vector<int>> validSubs;
// //         for (auto &sp : allSubs) {
// //             VRP currSP = {};
// //             currSP.list = sp;
// //             currSP.numVehicles = nVehicles;
// //             currSP.master = master;
// //             currSP.proc = proc;
// //             currSP.pid = pid;
// //             currSP.originalP = size;
// //             currSP.originalV = vehicles;
// //             currSP.printPaths = printPaths;
// //             if (nVehicles <= (int) sp.size()) {
// //                 if ((int) (hash<VRP>{}(currSP) % proc) == pid) {
// //                     subsetSolve(currSP, routeTable, matrix, size);
// //                 }
// //                 validSubs.push_back(sp);
// //             }
// //         }
// //         allSubs = validSubs;
// //         //sync
// //     }


// //     VRP prob;
// //     vector<int> allPoints;
// //     for(int i = 0; i < size; i++) {
// //         allPoints.push_back(i);
// //     }

// //     prob.list = allPoints;
// //     prob.numVehicles = vehicles;
// //     prob.master = master;
// //     prob.proc = proc;
// //     prob.pid = pid;
// //     prob.originalP = allPoints.size();
// //     prob.originalV = vehicles;
// //     prob.printPaths = printPaths;

// //     if ((int) (hash<VRP>{}(prob) % proc) == pid) {
// //         subsetSolve(prob, routeTable, matrix, size);
// //     } else {
// //         while (1) {
// //             fillRequests(prob, routeTable);
// //         }
// //     }

// //     for (int i = 0; i < size; i++)
// //         delete [] matrix[i];
// //     delete [] matrix;

// //     MPI_Finalize();
// // }
