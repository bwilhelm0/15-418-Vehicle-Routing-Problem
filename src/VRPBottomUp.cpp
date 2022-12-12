using namespace std;
#include "VRP.h"

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

vector<vector<int>> genWork(int N, int master) {
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

    vector<vector<int>> res;
    if (subgraphs.size() >= 2) {
        res.resize(subgraphs.size() - 2);
        copy(subgraphs.begin() + 1, subgraphs.end() - 1, res.begin());
    }

    return res;
}

// Generates all sub problem pairs for a given VRP (Num vehicles > 1)
vector<pair<VRP, VRP>> genSubs(VRP &prob) {
    vector<pair<VRP, VRP>> resSubs;

    vector<pair<vector<int>, vector<int>>> subsets;
    for (int i = 2; i < pow(2, prob.list.size() - 1); i+=2) { // check - 1
        vector<int> left;
        vector<int> right;

		for (int j = 0; (size_t)j < prob.list.size(); j++) {
            if (j == prob.master) {
                left.push_back(prob.master);
                right.push_back(prob.master);
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

    for (auto &p : subsets) {
        VRP next1;
        VRP next2;
        next1.list = p.first;
        next1.master = prob.master;
        next1.proc = prob.proc;
        next1.pid = prob.pid;
        next1.originalP = prob.originalP;
        next1.originalV = prob.originalV;
        next1.printPaths = prob.printPaths;
        next1.reduceComm = prob.reduceComm;

        next2.list = p.second;
        next2.master = prob.master;
        next2.proc = prob.proc;
        next2.pid = prob.pid;
        next2.originalP = prob.originalP;
        next2.originalV = prob.originalV;
        next2.printPaths = prob.printPaths;
        next2.reduceComm = prob.reduceComm;


        int maxlen = max(next1.list.size() - 1, next2.list.size() - 1);
        int minlen = min((int) next1.list.size() - 1, (int) next2.list.size() - 1);

        int smallVehicles = 1;
        int largeVehicles = prob.numVehicles - 1;

        if (largeVehicles > maxlen) {
            smallVehicles += largeVehicles - maxlen;
            largeVehicles = maxlen;
        }

        while (smallVehicles <= minlen && largeVehicles >= 1) {
            next1.numVehicles = (next1.list.size() < next2.list.size()) ? smallVehicles : largeVehicles;
            next2.numVehicles = (next1.list.size() < next2.list.size()) ? largeVehicles : smallVehicles;

            resSubs.push_back(make_pair(next1, next2));

            smallVehicles += 1;
            largeVehicles -= 1;
        }
    }

    return resSubs;
}

//checks for incoming data requests and responds
pair<vector<vector<int>>, vector<MPI_Request>> fillRequests(VRP &prob, unordered_map<VRP, VRPsolution> &solnMap) {
    //totalPoints + 1 because we send the number of vehicles as well
    const int requestSize = prob.originalP + 1;

    //In message there should be a 0 for each vehicle, this is largest possible message size
    const int answerSize = prob.originalP + prob.originalV + 1;
    const int routeAnswerSize = prob.originalP + prob.originalV;
    const int costAnswerSize = 1;

    //check for requests
    int *flags = new int[prob.proc];
    int *costFlags = new int[prob.proc];
    int *routeFlags = new int[prob.proc];
    for (int i = 0; i < prob.proc; i++) {
        if (i == prob.pid) continue;
        MPI_Iprobe(i, REQUEST, MPI_COMM_WORLD, &flags[i], MPI_STATUS_IGNORE);
        MPI_Iprobe(i, REQUEST_COST, MPI_COMM_WORLD, &costFlags[i], MPI_STATUS_IGNORE);
        MPI_Iprobe(i, REQUEST_ROUTES, MPI_COMM_WORLD, &routeFlags[i], MPI_STATUS_IGNORE);
    }

    vector<vector<int>> filledReqs;
    vector<MPI_Request> filledReqStatus;
    VRP req;

    for (int i = 0; i < prob.proc; i++) {
        if (i == prob.pid) continue;
        if (flags[i] != 0 || costFlags[i] != 0 || routeFlags[i] != 0) {
            MSG_TAG rTag = flags[i] ? REQUEST : (costFlags[i] ? REQUEST_COST : REQUEST_ROUTES);
            MSG_TAG aTag = flags[i] ? ANSWER : (costFlags[i] ? ANSWER_COST : ANSWER_ROUTES);
            req.list.resize(requestSize);
            MPI_Recv((void *) req.list.data(), requestSize, MPI_INT, i, rTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            req.numVehicles = req.list[requestSize - 1];
            for (int i = requestSize - 2; i > 0; i--) {
                if (req.list[i] != 0) { 
                    req.list.resize(i + 1);
                    break;
                }
            }


            unordered_map<VRP, VRPsolution>::const_iterator foundit = solnMap.find(req);

            vector<int> answerVal;
            int size = flags[i] ? answerSize : (costFlags[i] ? costAnswerSize : routeAnswerSize);
            answerVal.resize(size); 
            MPI_Request answerReq;

            if (foundit != solnMap.end()) { 
                // Found case, send solution
                if (flags[i] || costFlags[i]) answerVal[0] = foundit->second.time;

                if (flags[i] || routeFlags[i]) {
                    int pos = flags[i] ? 1 : 0;
                    for (auto &journey : foundit->second.routes) {
                        copy(journey.begin(), journey.end(), answerVal.begin() + pos);
                        pos += journey.size() - 1;
                    }
                }

                MPI_Isend((void *) answerVal.data(), answerSize, MPI_INT, i, aTag, MPI_COMM_WORLD, &answerReq);
            } else { 
                MPI_Isend((void *) answerVal.data(), answerSize, MPI_INT, i, aTag, MPI_COMM_WORLD, &answerReq);
            }
            filledReqs.push_back(answerVal);
            filledReqStatus.push_back(answerReq);
        }
    }
    return make_pair(filledReqs, filledReqStatus);
}


VRPsolution reqSoln(VRP &prob, unordered_map<VRP, VRPsolution> &solnMap, MSG_TAG tag) {
    const int requestSize = prob.originalP + 1;  //totalPoints + 1 because we send the number of vehicles as well
    const int answerSize = tag == REQUEST ? prob.originalP + prob.originalV + 1 : (tag == REQUEST_COST ? 1 : prob.originalP + prob.originalV);   //In message there should be a 0 for each vehicle, this is largest possible message size

    MSG_TAG answerTag = (tag == REQUEST ? ANSWER : (tag == REQUEST_COST ? ANSWER_COST : ANSWER_ROUTES));

    vector<MPI_Request> filledReqStatus;
    vector<vector<int>> filledReqs;

    int hashedProc = hash<VRP>{}(prob) % prob.proc;

    vector<int> request;
    request.resize(requestSize);
    vector<int> answer;
    answer.resize(answerSize);

    copy(prob.list.begin(), prob.list.end(), request.begin());
    request[requestSize-1] = prob.numVehicles;

    // Request solution to problem
    // TODO turn to isend
    MPI_Send((void *) request.data(), requestSize, MPI_INT, hashedProc, tag, MPI_COMM_WORLD);

    // Asynch Receive solution if it has been solved, keep checking for sends while waiting to avoid deadlock
    MPI_Request recRequest;
    MPI_Irecv((void *) answer.data(), answerSize, MPI_INT, hashedProc, answerTag, MPI_COMM_WORLD, &recRequest);

    //manually check if receive request was completed
    int recFlag = 0;
    MPI_Request_get_status(recRequest, &recFlag, MPI_STATUS_IGNORE);

    // Fill requests until the packet has been received
    while (!recFlag) {
        pair<vector<vector<int>>, vector<MPI_Request>> newFills = fillRequests(prob, solnMap);

        for (auto &path : newFills.first) {
            filledReqs.push_back(path);
        }

        for (auto &status : newFills.second) {
            filledReqStatus.push_back(status);
        }

        MPI_Request_get_status(recRequest, &recFlag, MPI_STATUS_IGNORE);
    }

        
    VRPsolution currSolution;
    if (tag == REQUEST || tag == REQUEST_COST) currSolution.time = answer[0];
    else {
        unordered_map<VRP, VRPsolution>::const_iterator got = solnMap.find(prob);

        // if solution already in hashtable, find and return VRPsolution
        if (got != solnMap.end()) {
            currSolution.time = got->second.time;
        } 
        // Otherwise request it
        else {
            cout << "Error, should be found in this case" << endl;
        }
    }
    
    vector<int> path;
    int masters = 0;
    // Wont enter this loop if the request is just for cost
    for (int i = (tag == REQUEST_ROUTES ? 0 : 1); i < answerSize; i++) {
        if (answer[i] == prob.master) {
            if (masters != 0) {
                path.push_back(answer[i]);
                currSolution.routes.push_back(path);
                path.resize(0);         //resize 0 because we insert
            }
            masters += 1;

            if (masters > prob.numVehicles) {
                break;
            }
        }
        path.push_back(answer[i]);
    }

 

    // Add solution to hash table
    solnMap[prob] = currSolution;
    MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(), MPI_STATUSES_IGNORE);
    return currSolution;
}


pair<vector<vector<int>>, vector<MPI_Request>> syncLevel(VRP &prob, unordered_map<VRP, VRPsolution> &solnMap) {
    vector<vector<int>> filledReqs;
    vector<MPI_Request> filledReqStatus;
    vector<int> sVec;
    sVec.resize(1);
    filledReqs.push_back(sVec);

    MPI_Request syncReq;
    for (int i = 0; i < prob.proc; i++) { 
        if (i == prob.pid) continue;
        MPI_Isend((void *) sVec.data(), 1, MPI_INT, i, LEVEL_SYNC, MPI_COMM_WORLD, &syncReq);
        filledReqStatus.push_back(syncReq);
    }

    vector<int> readySync;
    readySync.resize(prob.proc);
    int syncReady = 0;
    while (!syncReady) {
        //check for requests
        int *syncFlags = new int[prob.proc];
        for (int i = 0; i < prob.proc; i++) {
            if (i == prob.pid) continue;
            MPI_Iprobe(i, LEVEL_SYNC, MPI_COMM_WORLD, &syncFlags[i], MPI_STATUS_IGNORE);
        }

        VRP req;

        syncReady = 1;
        for (int i = 0; i < prob.proc; i++) {
            if (i == prob.pid) continue;

            if (syncFlags[i] != 0) {
                req.list.resize(1);
                MPI_Recv((void *) req.list.data(), 1, MPI_INT, i, LEVEL_SYNC, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                readySync[i] = 1;
            } else if (readySync[i] != 1) {
                syncReady = 0;
            }
        }
        pair<vector<vector<int>>, vector<MPI_Request>> newFills = fillRequests(prob, solnMap);
        filledReqs.insert(filledReqs.end(), newFills.first.begin(), newFills.first.end());
        filledReqStatus.insert(filledReqStatus.end(), newFills.second.begin(), newFills.second.end());
    }

    return make_pair(filledReqs, filledReqStatus);
}



VRPsolution vrpSolve(VRP &prob, unordered_map<VRP, VRPsolution> &solnMap, int **matrix, int size) {
    MSG_TAG reqTag = prob.reduceComm ? REQUEST_COST : REQUEST;

    VRPsolution res;
    pair<vector<vector<int>>, vector<MPI_Request>> currentlyFilling = fillRequests(prob, solnMap);
    vector<MPI_Request> filledReqStatus = currentlyFilling.second;
    vector<vector<int>> filledReqs = currentlyFilling.first;

    // if the problem only has 1 vehicle, then calculate tsp
    if (prob.numVehicles == 1) {
        res = tspSolve(matrix, size, prob.list, prob.printPaths);
        solnMap.insert({prob, res});
        MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(), MPI_STATUSES_IGNORE);
        return res;
    }

    vector<pair<VRP, VRP>> pairSubs = genSubs(prob);
    pair<VRP, VRP> finalPair;
    bool reqFirst = false; 
    bool reqSecond = false;
    res.time = INF;
    while (pairSubs.size() > 0) {
        vector<pair<VRP, VRP>> unfinished;
        for (auto &pairSub : pairSubs) {
            VRPsolution soln1;
            VRPsolution soln2;

            unordered_map<VRP, VRPsolution>::const_iterator gotFirst = solnMap.find(pairSub.first);
            unordered_map<VRP, VRPsolution>::const_iterator gotSecond = solnMap.find(pairSub.second);

            // if solution already in hashtable, find and return VRPsolution
            if (gotFirst != solnMap.end()) {
                soln1 = gotFirst->second;
            } 
            // Otherwise request it
            else {
                soln1 = reqSoln(pairSub.first, solnMap, reqTag);
            }
            
            if (gotSecond != solnMap.end()) {
                soln2 = gotSecond->second;
            } 
            else {
                soln2 = reqSoln(pairSub.second, solnMap, reqTag);
            }

            if (soln1.time == 0 || soln2.time == 0) {
                unfinished.push_back(pairSub);
                cout << "FAILED to communicate" << endl;
                continue;
            }

            if (res.time > max(soln1.time, soln2.time)) {
                finalPair = pairSub;
                reqFirst = soln1.routes.size() == 0;
                reqSecond = soln2.routes.size() == 0;
                soln1.routes.insert(soln1.routes.end(), soln2.routes.begin(), soln2.routes.end());
                res.routes = soln1.routes;
                res.time = max(soln1.time, soln2.time);
            }
        }
        pairSubs = unfinished;
    }

    if (reqFirst) {
        VRPsolution first = reqSoln(finalPair.first, solnMap, REQUEST_ROUTES);
        res.routes.insert(res.routes.end(), first.routes.begin(), first.routes.end());
    }
    
    if (reqSecond) {
        VRPsolution second = reqSoln(finalPair.second, solnMap, REQUEST_ROUTES);
        res.routes.insert(res.routes.end(), second.routes.begin(), second.routes.end());
    }

    solnMap[prob] = res;
    MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(), MPI_STATUSES_IGNORE);
    return res;
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

    StartupOptions config = parseOptions(argc, argv);

    int size = config.nodes;
    int vehicles = config.vehicles;
    int master = 0; // Keep at 0 until debugged
    bool printPaths = config.printPaths;
    bool timeLevels = config.timeLevels;

    // Seed random matrix generation
    if (config.seed == -1) {
        srand(time(NULL));
    } else {
        srand(config.seed);
    }

    vector<int> mVec;
    mVec.resize(size * size);
    if (pid == 0) {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; j++) {
                mVec[i * size + j] = (i == j) ? INF : (rand() % 20) + 1;
            }
        }
    }
    MPI_Bcast(mVec.data(), size * size, MPI_INT, 0, MPI_COMM_WORLD);

    int** matrix = new int*[size];
    for (int i = 0; i < size; ++i) {
        matrix[i] = new int[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = mVec[i * size + j];
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    Timer totalTime;

    //verify no empty set and no full set

    vector<vector<int>> allSubs = genWork(size, master);
    unordered_map<VRP, VRPsolution> routeTable;

    VRP prob;
    vector<int> allPoints;
    for(int i = 0; i < size; i++) {
        allPoints.push_back(i);
    }

    prob.list = allPoints;
    prob.numVehicles = vehicles;
    prob.master = master;
    prob.proc = proc;
    prob.pid = pid;
    prob.originalP = allPoints.size();
    prob.originalV = vehicles;
    prob.printPaths = printPaths;
    prob.reduceComm = config.reduceComm;

    double *levelTimes = new double[vehicles - 1];
    for (int nVehicles = 1; nVehicles < vehicles; nVehicles++) {
        vector<vector<int>> validSubs;
        Timer levelTime;
        
        for (auto &sp : allSubs) {
            //cout << pid << " " << sp.size() << endl;
            if (nVehicles <= (int) sp.size()) {
                VRP currSP; 
                currSP.master = master;
                currSP.numVehicles = nVehicles;
                currSP.proc = proc;
                currSP.pid = pid;
                currSP.originalP = size;
                currSP.originalV = vehicles;
                currSP.printPaths = printPaths;
                currSP.list = sp;
                currSP.reduceComm = config.reduceComm;
                if ((int) (hash<VRP>{}(currSP) % proc) == pid) {
                    vrpSolve(currSP, routeTable, matrix, size);
                }
                validSubs.push_back(sp);
            }
        }
        allSubs = validSubs;

        // Visualize workload imbalance
        if (timeLevels) levelTimes[nVehicles - 1] = levelTime.elapsed();

        pair<vector<vector<int>>, vector<MPI_Request>> newFills = syncLevel(prob, routeTable);
        MPI_Waitall(newFills.second.size(), newFills.second.data(), MPI_STATUSES_IGNORE);

    }

    if ((int) (hash<VRP>{}(prob) % proc) == pid) {
        VRPsolution res = vrpSolve(prob, routeTable, matrix, size);
        cout << "Time Cost is " << res.time << endl;
        printRoutes(res.routes);
        double elapsedTime = totalTime.elapsed();
        cout << "Took " << elapsedTime << " seconds" << endl;
        if (timeLevels) cout << "Level timing info: " << endl;
    } 
    pair<vector<vector<int>>, vector<MPI_Request>> newFills = syncLevel(prob, routeTable);
    MPI_Waitall(newFills.second.size(), newFills.second.data(), MPI_STATUSES_IGNORE);

    if (timeLevels) {
        for (int i = 0; i < vehicles - 1; i++) {
            if (pid == 0) cout << "\nLevel " << i + 1 << endl;
            MPI_Barrier(MPI_COMM_WORLD);

            cout << pid << ": " << levelTimes[i] << endl;
        }
    }

    delete levelTimes;

    MPI_Finalize();
}