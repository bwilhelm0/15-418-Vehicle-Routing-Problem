using namespace std;
#include "VRP.h"

// Function to print the final routes taken to solve the final problem
void printRoutes(vector<vector<point>> routes) {
  int numRoutes = (int)routes.size();
  for (int route = 0; route < numRoutes; route++) {
    cout << "Vehicle: " << route << ", ";
    int nodes = (int)routes[route].size();
    for (int node = 0; node < nodes - 1; node++) {
      cout << (int)routes[route][node] << " -> ";
    }
    cout << (int)routes[route][nodes - 1] << endl;
  }
}

void printList(vector<point> list) {
  for (auto &n : list) {
    cout << n << " ";
  }
  cout << endl;
}

// Function to generate the starting subsets based on the N nodes
vector<vector<point>> genWork(point N, point master) {
  vector<point> list;
  list.push_back(
      master); // The master node acts as a depot and is always on the list

  vector<vector<point>> subgraphs;
  subgraphs.push_back(list);

  // iterate through creating the 2^(N - 1) subsets
  for (point i = 0; i < N; i++) {
    if (i != master) {
      vector<vector<point>> tempgraphs;
      for (auto &subgraph : subgraphs) {
        tempgraphs.push_back(subgraph);
        subgraph.push_back(i);
        tempgraphs.push_back(subgraph);
      }
      subgraphs = tempgraphs;
    }
  }

  // Remove the empty set and the full set, not necessary for bottom up approach
  vector<vector<point>> res;
  if (subgraphs.size() >= 2) {
    res.resize(subgraphs.size() - 2);
    copy(subgraphs.begin() + 1, subgraphs.end() - 1, res.begin());
  }

  return res;
}

// Generates all sub problem pairs for a given VRP (Num vehicles > 1)
vector<pair<VRP, VRP>> genSubs(VRP &prob, VRPspecs &info) {
  vector<pair<VRP, VRP>> resSubs;

  // Create every possible bit pattern of N subproblems, split into two sets
  vector<pair<vector<point>, vector<point>>> subsets;
  for (int i = 2; i < pow(2, prob.list.size() - 1); i += 2) {
    vector<point> left;
    vector<point> right;

    for (int j = 0; (size_t)j < prob.list.size(); j++) {
      if (j == info.master) {
        left.push_back(info.master);
        right.push_back(info.master);
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

  // Create VRP structs and assign all possible vehicle combos
  for (auto &p : subsets) {
    VRP next1;
    VRP next2;
    next1.list = p.first;
    next2.list = p.second;

    int maxlen = max(next1.list.size() - 1, next2.list.size() - 1);
    int minlen = min((int)next1.list.size() - 1, (int)next2.list.size() - 1);

    point smallVehicles = 1;
    point largeVehicles = prob.numVehicles - 1;

    if (largeVehicles > maxlen) {
      smallVehicles += largeVehicles - maxlen;
      largeVehicles = maxlen;
    }

    // Assign all possible vehicle combos
    while (smallVehicles <= minlen && largeVehicles >= 1) {
      next1.numVehicles = (next1.list.size() < next2.list.size())
                              ? smallVehicles
                              : largeVehicles;
      next2.numVehicles = (next1.list.size() < next2.list.size())
                              ? largeVehicles
                              : smallVehicles;

      resSubs.push_back(make_pair(next1, next2));

      smallVehicles += 1;
      largeVehicles -= 1;
    }
  }

  return resSubs;
}

// checks for incoming data requests and responds
pair<vector<vector<point>>, vector<MPI_Request>>
fillRequests(VRPspecs &info, unordered_map<VRP, VRPsolution> &solnMap) {
  // totalPoints + 1 because we send the number of vehicles as well
  const int requestSize = info.originalP + 1;

  // In message there should be a 0 for each vehicle, this is largest possible
  // message size
  const int routeAnswerSize = info.originalP + info.originalV;
  const int costAnswerSize = INT2P;
  const int answerSize = routeAnswerSize + costAnswerSize;

  // check for requests
  int *flags = new int[info.proc];
  int *costFlags = new int[info.proc];
  int *routeFlags = new int[info.proc];
  for (int i = 0; i < info.proc; i++) {
    if (i == info.pid)
      continue;
    // Look for all message request types
    MPI_Iprobe(i, REQUEST, MPI_COMM_WORLD, &flags[i], MPI_STATUS_IGNORE);
    MPI_Iprobe(i, REQUEST_COST, MPI_COMM_WORLD, &costFlags[i],
               MPI_STATUS_IGNORE);
    MPI_Iprobe(i, REQUEST_ROUTES, MPI_COMM_WORLD, &routeFlags[i],
               MPI_STATUS_IGNORE);
  }

  vector<vector<point>> filledReqs;
  vector<MPI_Request> filledReqStatus;
  VRP req;

  // Check all processor connections, and answer requests accordingly
  for (int i = 0; i < info.proc; i++) {
    if (i == info.pid)
      continue;
    if (flags[i] != 0 || costFlags[i] != 0 || routeFlags[i] != 0) {
      MSG_TAG rTag =
          flags[i] ? REQUEST : (costFlags[i] ? REQUEST_COST : REQUEST_ROUTES);
      MSG_TAG aTag =
          flags[i] ? ANSWER : (costFlags[i] ? ANSWER_COST : ANSWER_ROUTES);

      // Receive request from specific message type
      req.list.resize(requestSize);
      MPI_Recv((void *)req.list.data(), requestSize * sizeof(point), MPI_CHAR,
               i, rTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Resize the request for lookup
      req.numVehicles = req.list[requestSize - 1];
      for (int i = requestSize - 2; i > 0; i--) {
        if (req.list[i] != 0) {
          req.list.resize(i + 1);
          break;
        }
      }

      unordered_map<VRP, VRPsolution>::const_iterator foundit =
          solnMap.find(req);

      // Resize answer vector based on message type
      vector<point> answerVal;
      int size = flags[i] ? answerSize
                          : (costFlags[i] ? costAnswerSize : routeAnswerSize);
      answerVal.resize(size);
      MPI_Request answerReq;

      if (foundit != solnMap.end()) {
        // Found case, send solution
        if (flags[i] || costFlags[i])
          ((int *)answerVal.data())[0] = foundit->second.time;

        if (flags[i] || routeFlags[i]) {
          int pos = flags[i] ? INT2P : 0;
          for (auto &journey : foundit->second.routes) {
            copy(journey.begin(), journey.end(), answerVal.begin() + pos);
            pos += journey.size() - 1;
          }
        }

        // Should change this to MPI_CHAR and multiply size by sizeof point
        MPI_Isend((void *)answerVal.data(), size * sizeof(point), MPI_CHAR, i,
                  aTag, MPI_COMM_WORLD, &answerReq);
      } else {
        // Send empty array if not found. This is an error
        cout << "Error, shouldn't be requesting here" << endl;
        MPI_Isend((void *)answerVal.data(), size * sizeof(point), MPI_CHAR, i,
                  aTag, MPI_COMM_WORLD, &answerReq);
      }
      filledReqs.push_back(answerVal);
      filledReqStatus.push_back(answerReq);
    }
  }

  // Return answer info so it doesn't go out of scope
  return make_pair(filledReqs, filledReqStatus);
}

// Function for formatting a request and responding while waiting
VRPsolution reqSoln(VRP &prob, VRPspecs &info,
                    unordered_map<VRP, VRPsolution> &solnMap, MSG_TAG tag) {
  info.requests++;
  const int requestSize =
      info.originalP +
      1; // totalPoints + 1 because we send the number of vehicles as well
  const int answerSize =
      tag == REQUEST
          ? info.originalP + info.originalV + INT2P
          : (tag == REQUEST_COST
                 ? INT2P
                 : info.originalP +
                       info.originalV); // In message there should be a 0 for
                                        // each vehicle, this is largest
                                        // possible message size

  // Set the answer tag based on the request
  MSG_TAG answerTag =
      (tag == REQUEST ? ANSWER
                      : (tag == REQUEST_COST ? ANSWER_COST : ANSWER_ROUTES));

  vector<MPI_Request> filledReqStatus;
  vector<vector<point>> filledReqs;

  // Find who owns the data
  int hashedProc = hash<VRP>{}(prob) % info.proc;

  vector<point> request;
  request.resize(requestSize);
  vector<point> answer;
  answer.resize(answerSize);

  // Format the request with number of vehicles at the end
  copy(prob.list.begin(), prob.list.end(), request.begin());
  request[requestSize - 1] = prob.numVehicles;

  // Request solution to problem
  MPI_Request myReq; // Unused because it has to be filled to leave function
  MPI_Isend((void *)request.data(), requestSize * sizeof(point), MPI_CHAR,
            hashedProc, tag, MPI_COMM_WORLD, &myReq);
  filledReqs.push_back(request);
  filledReqStatus.push_back(myReq);

  // Asynch Receive solution if it has been solved, keep checking for sends
  // while waiting to avoid deadlock
  MPI_Request recRequest;
  MPI_Irecv((void *)answer.data(), answerSize * sizeof(point), MPI_CHAR,
            hashedProc, answerTag, MPI_COMM_WORLD, &recRequest);
  filledReqs.push_back(answer);

  // manually check if receive request was completed
  int recFlag = 0;
  MPI_Request_get_status(recRequest, &recFlag, MPI_STATUS_IGNORE);

  // Fill requests until the packet has been received
  while (!recFlag) {
    pair<vector<vector<point>>, vector<MPI_Request>> newFills =
        fillRequests(info, solnMap);

    filledReqs.insert(filledReqs.end(), newFills.first.begin(),
                      newFills.first.end());
    filledReqStatus.insert(filledReqStatus.end(), newFills.second.begin(),
                           newFills.second.end());
    MPI_Request_get_status(recRequest, &recFlag, MPI_STATUS_IGNORE);
  }

  VRPsolution currSolution;
  // Set min time, or look it up if you know you have it
  if (tag == REQUEST || tag == REQUEST_COST)
    currSolution.time = ((int *)answer.data())[0];
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

  vector<point> path;
  point masters = 0;
  // Wont enter this loop if the request is just for cost
  // Format routes
  for (int i = (tag == REQUEST_ROUTES ? 0 : INT2P); i < answerSize; i++) {
    if (answer[i] == info.master) {
      if (masters != 0) {
        path.push_back(answer[i]);
        currSolution.routes.push_back(path);
        path.resize(0); // resize 0 because we insert
      }
      masters += 1;

      if (masters > prob.numVehicles) {
        break;
      }
    }
    path.push_back(answer[i]);
  }

  // Add solution to hash table, possibly overriding existing routes
  if (info.saveRequests)
    solnMap[prob] = currSolution;
  MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(),
              MPI_STATUSES_IGNORE);
  return currSolution;
}

// Function for telling everyone a level is done, and filling requests while
// waiting
// TODO revise to sync in a ring, 2 * N communications instead of N^2
// TODO Modify into ring reduce
void syncLevel(VRPspecs &info, unordered_map<VRP, VRPsolution> &solnMap) {
  vector<vector<point>> filledReqs;
  vector<MPI_Request> filledReqStatus;
  vector<point> sVec;

  // Message doesn't matter
  sVec.resize(1);
  filledReqs.push_back(sVec);

  // Send to everyone that pid is finished with level
  for (int i = 0; i < info.proc; i++) {
    if (i == info.pid)
      continue;
    MPI_Request syncReq;
    MPI_Isend((void *)sVec.data(), 0, MPI_CHAR, i, LEVEL_SYNC, MPI_COMM_WORLD,
              &syncReq);
    filledReqStatus.push_back(syncReq);
  }

  vector<int> readySync;
  readySync.resize(info.proc);
  int syncReady = 0;
  // Wait until everyone is ready for sync
  while (!syncReady) {
    // check for syncs
    vector<int> syncFlags;
    syncFlags.resize(info.proc);
    for (int i = 0; i < info.proc; i++) {
      if (i == info.pid)
        continue;
      MPI_Iprobe(i, LEVEL_SYNC, MPI_COMM_WORLD, &syncFlags[i],
                 MPI_STATUS_IGNORE);
    }

    VRP req;

    syncReady = 1;
    for (int i = 0; i < info.proc; i++) {
      if (i == info.pid)
        continue;
      if (syncFlags[i] != 0) {
        req.list.resize(1);
        MPI_Recv((void *)req.list.data(), 0, MPI_CHAR, i, LEVEL_SYNC,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        readySync[i] = 1;
      } else if (readySync[i] != 1) {
        syncReady = 0;
      }
    }
    pair<vector<vector<point>>, vector<MPI_Request>> newFills =
        fillRequests(info, solnMap);
    filledReqs.insert(filledReqs.end(), newFills.first.begin(),
                      newFills.first.end());
    filledReqStatus.insert(filledReqStatus.end(), newFills.second.begin(),
                           newFills.second.end());
  }

  MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(),
              MPI_STATUSES_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
}

// Function for reducing the mincosts in the hash tables to every processor
pair<vector<vector<point>>, vector<MPI_Request>>
ringReduce(VRPspecs &info, unordered_map<VRP, VRPsolution> &solnMap,
           unordered_map<VRP, VRPsolution> &newSolns) {
  vector<vector<point>> filledReqs;
  vector<MPI_Request> filledReqStatus;

  const int requestSize = info.originalP + 1;

  int destProc = (info.proc + info.pid + 1) % info.proc;
  int srcProc = (info.proc + info.pid - 1) % info.proc;

  // Message will contain size new hash table entries
  int myValSize = newSolns.size();

  MPI_Request sizeReq;
  // Send to PID 1 above with the size of hash table entries to be communicated
  MPI_Isend((void *)&myValSize, 1, MPI_INT, destProc, REDUCE, MPI_COMM_WORLD,
            &sizeReq);
  filledReqStatus.push_back(sizeReq);

  vector<point> myVals;
  myVals.resize(myValSize * (requestSize + INT2P));

  // Format data as vector..., numVehicles, time
  int pos = 0;
  for (auto it = newSolns.begin(); it != newSolns.end(); ++it) {
    copy(it->first.list.begin(), it->first.list.end(), myVals.begin() + pos);
    *((int *)(&myVals[pos + requestSize])) = it->second.time;
    myVals[pos + requestSize - 1] = it->first.numVehicles;
    pos += requestSize + INT2P;

    // Add newSolns into existing solutions, overwriting because these are
    // complete and they shouldnt be present
    solnMap[it->first] = it->second;
  }

  MPI_Request sendRequest;
  // Send to PID 1 above, these are the hash table entries
  MPI_Isend((void *)myVals.data(),
            myValSize * (requestSize + INT2P) * sizeof(point), MPI_CHAR,
            destProc, REDUCE_DATA, MPI_COMM_WORLD, &sendRequest);

  int newData = 0;
  vector<point> req;
  int reqSize;
  int sendFlag = 0;
  int levels = info.steps == -1 ? info.proc - 1 : info.steps;

  // Iterate through levels of ring reduce
  for (int level = 0; level < levels; level++) {
    while (!newData) {
      pair<vector<vector<point>>, vector<MPI_Request>> newFills =
          fillRequests(info, solnMap);
      filledReqs.insert(filledReqs.end(), newFills.first.begin(),
                        newFills.first.end());
      filledReqStatus.insert(filledReqStatus.end(), newFills.second.begin(),
                             newFills.second.end());

      // Check for incoming new data
      MPI_Iprobe(srcProc, REDUCE + 2 * level, MPI_COMM_WORLD, &newData,
                 MPI_STATUS_IGNORE);
    }
    newData = 0;

    // Receive size of data and resize receive buffer
    if (level % 2 == 0) {
      MPI_Recv((void *)&reqSize, 1, MPI_INT, srcProc, REDUCE + 2 * level,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      req.resize(reqSize * (requestSize + INT2P));
    } else {
      MPI_Recv((void *)&myValSize, 1, MPI_INT, srcProc, REDUCE + 2 * level,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      myVals.resize(myValSize * (requestSize + INT2P));
    }

    // totalAdded += reqSize;
    while (!newData) {
      pair<vector<vector<point>>, vector<MPI_Request>> newFills =
          fillRequests(info, solnMap);
      filledReqs.insert(filledReqs.end(), newFills.first.begin(),
                        newFills.first.end());
      filledReqStatus.insert(filledReqStatus.end(), newFills.second.begin(),
                             newFills.second.end());

      // Check for incoming new data
      MPI_Iprobe(srcProc, REDUCE_DATA + 2 * level, MPI_COMM_WORLD, &newData,
                 MPI_STATUS_IGNORE);
    }
    newData = 0;

    if (level % 2 == 0) {
      // Receive actual data
      MPI_Recv((void *)req.data(),
               reqSize * (requestSize + INT2P) * sizeof(point), MPI_CHAR,
               srcProc, REDUCE_DATA + 2 * level, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    } else {
      // Receive size of data and resize receive buffer
      MPI_Recv((void *)myVals.data(),
               myValSize * (requestSize + INT2P) * sizeof(point), MPI_CHAR,
               srcProc, REDUCE_DATA + 2 * level, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }

    // Fill requests until the send has completed and the level is finished
    while (!sendFlag) {
      pair<vector<vector<point>>, vector<MPI_Request>> newFills =
          fillRequests(info, solnMap);
      filledReqs.insert(filledReqs.end(), newFills.first.begin(),
                        newFills.first.end());
      filledReqStatus.insert(filledReqStatus.end(), newFills.second.begin(),
                             newFills.second.end());

      MPI_Request_get_status(sendRequest, &sendFlag, MPI_STATUS_IGNORE);
    }
    sendFlag = 0;

    // Copy data into hash table
    int maxSize = level % 2 == 0 ? reqSize : myValSize;
    vector<point> temp = level % 2 == 0 ? req : myVals;
    for (int ind = 0; ind < (int)(maxSize * (requestSize + INT2P));
         ind += (requestSize + INT2P)) {
      VRP incoming;
      VRPsolution inSoln;

      inSoln.time = *((int *)(&temp[ind + requestSize]));
      incoming.numVehicles = temp[ind + requestSize - 1];

      // This needs to be -2 because this is the last index of the route
      int last = ind + requestSize - 2;
      for (; last > ind; last--) {
        if (temp[last] != 0)
          break;
      }

      incoming.list.insert(incoming.list.end(), temp.begin() + ind,
                           temp.begin() + last + 1);

      // We use insert here because we want it to fail if there is a solution
      // with routes already in the table
      solnMap.insert({incoming, inSoln});
    }

    if (level != levels - 1) {
      MPI_Request a;
      if (level % 2 == 0) {
        MPI_Isend((void *)&reqSize, 1, MPI_INT, destProc,
                  REDUCE + 2 * (level + 1), MPI_COMM_WORLD, &a);
        // Send to PID 1 above with the size of hash table entries to be
        // communicated
        MPI_Isend((void *)req.data(),
                  reqSize * (requestSize + INT2P) * sizeof(point), MPI_CHAR,
                  destProc, REDUCE_DATA + 2 * (level + 1), MPI_COMM_WORLD,
                  &sendRequest);
      } else {
        MPI_Isend((void *)&myValSize, 1, MPI_INT, destProc,
                  REDUCE + 2 * (level + 1), MPI_COMM_WORLD, &a);
        // Send to PID 1 above with the size of hash table entries to be
        // communicated
        MPI_Isend((void *)myVals.data(),
                  myValSize * (requestSize + INT2P) * sizeof(point), MPI_CHAR,
                  destProc, REDUCE_DATA + 2 * (level + 1), MPI_COMM_WORLD,
                  &sendRequest);
      }
      filledReqStatus.push_back(a);
    }
  }

  syncLevel(info, solnMap);

  return make_pair(filledReqs, filledReqStatus);
}

// Function generates all subsets and retrieves necessary info to solve them
VRPsolution vrpSolve(VRP &prob, VRPspecs &info,
                     unordered_map<VRP, VRPsolution> &solnMap,
                     unordered_map<VRP, VRPsolution> &newSolns, int **matrix,
                     int size) {
  MSG_TAG reqTag = info.reduceComm ? REQUEST_COST : REQUEST;

  VRPsolution res;
  // Fill some requests before we begin
  pair<vector<vector<point>>, vector<MPI_Request>> currentlyFilling =
      fillRequests(info, solnMap);
  vector<MPI_Request> filledReqStatus = currentlyFilling.second;
  vector<vector<point>> filledReqs = currentlyFilling.first;

  // if the problem only has 1 vehicle, then calculate tsp
  if (prob.numVehicles == 1) {
    res = tspSolve(matrix, size, prob.list, info.printPaths);

    if (info.ringReduce)
      newSolns[prob] = res;
    else
      solnMap[prob] = res;
    MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(),
                MPI_STATUSES_IGNORE);
    return res;
  }

  // Generate all subset pairs
  vector<pair<VRP, VRP>> pairSubs = genSubs(prob, info);

  // Remember the final pair in case we need to request routes
  pair<VRP, VRP> finalPair;
  bool reqFirst = false;
  bool reqSecond = false;

  // Set time to always be overridden
  res.time = INF;

  while (pairSubs.size() > 0) {
    vector<pair<VRP, VRP>> unfinished;
    for (auto &pairSub : pairSubs) {
      VRPsolution soln1;
      VRPsolution soln2;

      // Dont have to look up in newSolns because that is the current level of
      // DP, and every level is ind of same level
      unordered_map<VRP, VRPsolution>::const_iterator gotFirst =
          solnMap.find(pairSub.first);
      unordered_map<VRP, VRPsolution>::const_iterator gotSecond =
          solnMap.find(pairSub.second);

      // if solution already in hashtable, find VRPsolution
      if (gotFirst != solnMap.end()) {
        soln1 = gotFirst->second;
      }
      // Otherwise request it
      else {
        if (info.ringReduce && info.steps == -1)
          cout << "Error, should never be here" << endl;
        soln1 = reqSoln(pairSub.first, info, solnMap, reqTag);
      }

      if (gotSecond != solnMap.end()) {
        soln2 = gotSecond->second;
      } else {
        if (info.ringReduce && info.steps == -1)
          cout << "Error, should never be here" << endl;
        soln2 = reqSoln(pairSub.second, info, solnMap, reqTag);
      }

      if (soln1.time == 0 || soln2.time == 0) {
        unfinished.push_back(pairSub);
        cout << "FAILED to communicate" << endl;
        continue;
      }

      // Set the current lowest time and route information
      if (res.time > max(soln1.time, soln2.time)) {
        finalPair = pairSub;
        reqFirst = soln1.routes.size() == 0;
        reqSecond = soln2.routes.size() == 0;

        soln1.routes.insert(soln1.routes.end(), soln2.routes.begin(),
                            soln2.routes.end());
        res.routes = soln1.routes;
        res.time = max(soln1.time, soln2.time);
      }
    }
    pairSubs = unfinished;
  }

  // Request routes if they were not in the lookup table
  if (reqFirst) {
    VRPsolution first = reqSoln(finalPair.first, info, solnMap, REQUEST_ROUTES);
    res.routes.insert(res.routes.end(), first.routes.begin(),
                      first.routes.end());
  }
  if (reqSecond) {
    VRPsolution second =
        reqSoln(finalPair.second, info, solnMap, REQUEST_ROUTES);
    res.routes.insert(res.routes.end(), second.routes.begin(),
                      second.routes.end());
  }

  // Set result with mincost and routes
  if (info.ringReduce)
    newSolns[prob] = res;
  else
    solnMap[prob] = res;
  MPI_Waitall(filledReqStatus.size(), filledReqStatus.data(),
              MPI_STATUSES_IGNORE);
  return res;
}

int main(int argc, char *argv[]) {
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
  point vehicles = config.vehicles;
  point master = 0; // Keep at 0 until debugged
  bool printPaths = config.printPaths;
  bool timeLevels = config.timeLevels;

  if ((vehicles >= size || vehicles <= 0 || size <= 1) ||
      (config.steps != -1 && config.steps >= proc)) {
    if (pid == 0)
      cout << "error, invalid input" << endl;
    return 0;
  }

  // Seed random matrix generation
  if (config.seed == -1) {
    srand(time(NULL));
  } else {
    srand(config.seed);
  }

  // Generate matrix
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

  int **matrix = new int *[size];
  for (int i = 0; i < size; ++i) {
    matrix[i] = new int[size];
    for (int j = 0; j < size; j++) {
      matrix[i][j] = mVec[i * size + j];
    }
  }

  if (pid == 0 && config.printMatrix)
    printMatrix(matrix, size);

  MPI_Barrier(MPI_COMM_WORLD);
  Timer totalTime;

  vector<vector<point>> allSubs = genWork(size, master);
  unordered_map<VRP, VRPsolution> routeTable;

  VRP prob;
  VRPspecs info;
  vector<point> allPoints;
  for (point i = 0; i < size; i++) {
    allPoints.push_back(i);
  }

  prob.list = allPoints;
  prob.numVehicles = vehicles;

  info.master = master;
  info.proc = proc;
  info.pid = pid;
  info.originalP = allPoints.size();
  info.originalV = vehicles;
  info.printPaths = printPaths;
  info.reduceComm = config.reduceComm;
  info.ringReduce = config.ringReduce;
  info.steps = config.steps;
  info.saveRequests = config.saveRequests;

  double *levelTimes = new double[vehicles - 1];

  // Iterate through vehicle levels of DP
  for (point nVehicles = 1; nVehicles < vehicles; nVehicles++) {
    vector<vector<vector<point>>> procSubs;
    vector<vector<point>> validSubs;
    Timer levelTime;
    unordered_map<VRP, VRPsolution> newSolns;

    procSubs.resize(proc);
    for (auto &sp : allSubs) {
      if (nVehicles <= (point)sp.size()) {
        VRP currSP;
        currSP.numVehicles = nVehicles;
        currSP.list = sp;

        procSubs[(int)(hash<VRP>{}(currSP) % proc)].push_back(sp);
        validSubs.push_back(sp);
      }
    }

    // Iterate through relevant VRPs
    for (int i = 0; i < (int)procSubs[pid].size(); i++) {
      VRP currSP;
      currSP.numVehicles = nVehicles;
      currSP.list = procSubs[pid][i];

      vrpSolve(currSP, info, routeTable, newSolns, matrix, size);
      // We would've liked to implement work stealing here
    }

    allSubs = validSubs;

    // Visualize workload imbalance
    if (timeLevels)
      levelTimes[nVehicles - 1] = levelTime.elapsed();

    pair<vector<vector<point>>, vector<MPI_Request>> newFills;
    if (config.ringReduce)
      newFills = ringReduce(info, routeTable, newSolns);
    else
      syncLevel(info, routeTable);
    MPI_Waitall(newFills.second.size(), newFills.second.data(),
                MPI_STATUSES_IGNORE);
  }

  if ((int)(hash<VRP>{}(prob) % proc) == pid) {
    unordered_map<VRP, VRPsolution> newSolns;
    VRPsolution res = vrpSolve(prob, info, routeTable, newSolns, matrix, size);
    cout << "Time Cost is " << res.time << endl;
    printRoutes(res.routes);
    double elapsedTime = totalTime.elapsed();
    cout << "Took " << elapsedTime << " seconds" << endl;
    if (timeLevels)
      cout << "Level timing info: " << endl;
  }
  syncLevel(info, routeTable);
  MPI_Barrier(MPI_COMM_WORLD);

  if (timeLevels) {
    for (point i = 0; i < vehicles - 1; i++) {
      if (pid == 0)
        cout << "\nLevel " << i + 1 << endl;
      MPI_Barrier(MPI_COMM_WORLD);

      cout << pid << ": " << levelTimes[i] << endl;
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  if (config.printRequests) {
    if (pid == 0)
      cout << "\nRequest summary " << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    cout << pid << " Number of requests made " << info.requests << endl;
    MPI_Barrier(MPI_COMM_WORLD);
  }

  delete levelTimes;

  MPI_Finalize();
}