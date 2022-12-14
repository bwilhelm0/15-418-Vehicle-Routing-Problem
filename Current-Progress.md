Over the past few weeks we processed and reformatted uber data sets to be used with our eventual Vehicle Routing Problem (VRP) and we also implemented a sequential Traveling Salesman Problem that will be instrumental in our parallel implementation of the VRP. For the uber data processing we condensed the data and averaged out the distance values for trips that had the same starting and ending points, we then proceeded to reduce the precision of longitude and latitude values to lower the amount of nodes we would need to process within our graph. We also isolated the unique edges/trips and named each location with a node number so that we can create an adjacency matrix that we’ll need for our algorithms. We implemented a sequential traveling salesman problem solver that can be found within our github at TSP.cpp. We’ll be using this sequential implementation as a basis to create our parallel VRP solver.

Additionally, we’ve started a rough outline of our approach to workload balancing and we’re planning on implementing our approach in the upcoming week.

Schedule
Nov 30th-Dec 3rd: Kobe-Integrate Uber Data with TSP Algorithm, Brendan-Implement OpenMPI code for naive workload balance
Dec 4th-7th: Brendan-Analyze performance and identify bottlenecks, implement work stealing
Dec 8th-10th: Kobe-Create distributed data structure and ring reduce code
Dec 11th-14th: Brendan-Parallelize workload generation and consumption with lookup table
Dec 15th-17th: Brendan & Kobe-Create performance comparison charts for different implementations, create poster with algorithm description
