**TITLE: Parallel Vehicle Routing Problem in the United States **
Kobe Zhang
Brendan Wilhelm
**URL**: https://bwilhelm0.github.io/15-418-Vehicle-Routing-Problem/

**SUMMARY**: We are going to implement and optimize a solver for the vehicle routing problem and other forms of the traveling salesman problem using the PSC. We will measure cost using large datasets of real-world Uber data, and minimize the round-trip cost of a vehicle fleet that covers every node.

**BACKGROUND**: Given a large map with various nodes and weights on the edges, the vehicle routing problem must guarantee that each node is visited by at least one vehicle within the vehicle fleet. Our approach to this project would be to first optimize a single instance of the traveling salesman problem, requiring 1 vehicle to visit every single node while minimizing the cost. Then we would parallelize the problem by partitioning the vertices of the map into subsets and assigning each vehicle to travel through a unique subset.
  One way that this application can benefit from parallelism, especially on a large computing cluster, is that a work efficient sequential solution to the traveling salesman problem can be employed on every compute node. Then this allows for a brute force approach to be utilized due to the computing power, guaranteeing the optimal solution to the problem. The parallelism can also be useful when partitioning the map by vehicle, as there is a lot of processing involved in creating all of the different vehicle-map assignment possibilities.
  
**THE CHALLENGE**: One of the largest challenges that we can see for this problem is the preprocessing necessary to ensure workload balance. This could get very complex because there are many ways to divide the work, and finding the optimal workload balance will be tricky with a high processor count. For example, with a work stealing implementation, this would require a large amount of communication between processors to ensure a high load balance. If we then want to reduce communication and increase locality, we might try to reuse data so that in an MPI framework less particles need to be communicated. 
  Another factor that we will need to consider is that there is a data dependency between generating the tasks for each processor, and if not executed carefully could result in detrimental bottlenecks. One thought we had, since there is a large amount of processing to be done to create tasks for processors, is that we could employ some form of pipelining to reduce latency. A constraint that we will need to keep in mind is the large size of the dataset compared to the individual node’s memory size.

RESOURCES: For our project, we would need to access the Pittsburgh Supercomputer Cluster (PSC). For our approach to the vehicle routing problem, we will begin by referencing the following papers, although a single instance of the traveling salesman problem is not our focus as we will focus on parallelizing the Vehicle Routing Problem. We will reference the algorithms from these papers to implement our own solution.

**Braun, Heinrich. "On solving travelling salesman problems by genetic algorithms." International Conference on Parallel Problem Solving from Nature. Springer, Berlin, Heidelberg, 1990.
Ouaarab, Aziz, Belaïd Ahiod, and Xin-She Yang. "Discrete cuckoo search algorithm for the traveling salesman problem." Neural Computing and Applications 24.7 (2014): 1659-1669.**

**GOALS AND DELIVERABLES**
PLAN TO ACHIEVE:
Implement a parallelized version of the Vehicle Routing Problem in C++ . 
Process and utilize large Uber datasets to create a graph of the United States for the program to use.
Create a visual graph of the inputs and outputs of the program in python showing the results of the program and the changes in the solution over time as it works.
For the poster session, create speedup graphs after each optimization and also run the aforementioned python graphing as a live demo
Implement a fine-grained dynamic workload balancing scheme to improve speedup (i.e. work stealing, distributed work queues…etc.)
HOPE TO ACHIEVE:
Have 0.7x linear speedup over a single core parallel version, because we think this is reasonably difficult given the communication to computation ratio.
Utilize real-time Uber data using the Uber APIs
Explore a parallel implementation of the traveling merchant problem (Only possible if we are doing VERY well)

**PLATFORM CHOICE**: We chose to utilize the PSC because our parallel implementation would benefit from a large amount of CPUs. Our program does not seem to fit with the GPU architecture due to low amount of instruction level parallelism and variance in workload. We’ll be using C++ for the parallel implementation due to familiarity, speed of execution, and compatibility with Open MPI. We would also be using some python for data manipulation and displaying results.




SCHEDULE:

Week
Nov 13-19
Genetic algorithm research
Uber data processing 
Nov 20-26
Sequential traveling salesman problem and naive workload generation implementation
Nov 27-Dec 3
Analysis of initial implementation and improvements
Python visualization 
Dec 4-10
Continued optimization with improved workload balancing schemes
Dec 11-17
Final program iteration and analysis
Poster creation, live demo verification, finishing touches









