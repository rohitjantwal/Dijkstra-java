# Dijkstra-java
Dijkstra algorithm implementation in java

Name: Rohit Jantwal
UFID: 1351-4976
E-mail: rohitjantwal@ufl.edu

Compiler used: Java SE 7
How to compile: javac dijkstra.java

How to run:
For random mode:
java dijkstra –r n d x where ‘n’, ‘d’, and ‘x’ are integer values representing number of nodes, density percentage and source node respectively. For user input mode:
Simple scheme:
java dijkstra -s filename
where –s represents simple scheme mode and filename is the name of the file containing the input.

Fibonacci heap scheme:
java dijkstra –f filename where –f represents Fibonacci heap scheme and filename is the name of the file containing the input.

Assuming that that the vertices are from 0 to n-1. An example input can be:
0
3 3
0 1 3
1 2 8
0 2 1
This graph consists of three vertices {0,1,2] and three edges (0,1), (1,2) and (0,2) with costs 3,8 and 1 respectively.

Structure of the program:
 main(): Main method takes input from the command line and stores values in the corresponding variables. It also detects the type of mode to be run i.e. random mode or simple scheme mode. It also calculates the time for random mode scheme.
 randommode(): This defines the random mode function which accepts 3 variables: n, d, x. Initializes the number of vertices and density given by the user. This checks whether the graph is connected or not using is_connect function. Generates the required ((n*(n-1)/2)*density) number of edges. And finally prints the percentage density.
 simplescheme(): This defines simple scheme mode function which accept a single parameter Rnode() which inserts nodes adjacent to each other. Initializes a HashSet of thise vertices whose shortest path is not yet found. Finds the shortest path from a node to all undetermined nodes. Also checks if the distance between the nodes and the new added nodes is changed and returns shortest distance in the end.
 fiboscheme(): This defines Fibonacci scheme which accepts a single parameter Rnode() which inserts nodes adjacent to each other. Constructs a Fibonacci heap of the nodes. Checks for distance of the adjacent node that has just been removed. Checks if node has already been inserted, otherwise inserts it. If it is already in the heap then finds the corresponding node in the heap and lowers the distance. Returns the shortest distance graph in at the end.
 load_file(): This function accepts a file in the given format and passes the variables to get the source node, number of nodes and number of edges. And return the edges at the end.
 get_number_of_nodes(): This function accepts the edges as parameter and using the start and end of nodes calculates the number of nodes. Returns the number of nodes.
 create_graph(): This function accepts edges and integer as parameter and creates a graph accordingly using class Rnode. is_connect(): This function accepts graph as the parameter and checks the distance between two nodes. If the distance is larger than the maximum possible distance then the graph is not connected. This returns true if the nodes are found to be connected.
 Class Vertex: Constructor in this class lets us input the start node, end node and the distance. Also checks for the start node being the end node itself.
 Class Adjacent_node: The constructor of this class inputs the position of the node and the distance of the node from the source.
 Class Rnode: This class maintains the list of adjacent nodes. Sets the neighbor node of a node and sets the index of the nodes.
 Class Fibo_heap: The constructor of this class invokes the input and output function for index, degree, distance, parent, child, left sibling, right sibling and child cut of a node.

link(): This function links the defines the relation of child and parent between nodes.
merge(): This function is used to merge two Fibonacci heaps and set root and child nodes accordingly.
insert(): This function accepts a integer and a double value as parameters. This is used to insert a new node into a heap or merge it into the given heap.
placenode():This is used in remove min function and it links the nodes having same degree in root list.
removemin(): This function returns the node with minimum value and adjusts the root list accordingly using placenode() to adjust the position of min node.
lowerkey(): This function checks if the current key is less than the new key.
cut(): This function accepts two parameters which are nodes themselves. This is used to divide two nodes and that are currently joined together in the Fibonacci heap and add them to the root list.
setcut(): This is a complimentary function of cut child and arranges the other node whose child was cut in the cut() function and then add it to the Fibonacci heap.
delete(): This is used to delete node from the heap.
searchnode(): This is used to search a node by its root value. Returns the node that has been found.
print(): This accepts the root node as parameter and returns the heap whose root node is the given node.
Result comparison (Expected):

Simple scheme: O( n2 ).
Fibonacci heap: O(nlogn + e).
This shows that for small values of ‘n’ simple scheme should give better results but as n increases, Fibonacci scheme will start performing better. Fibonacci scheme works efficiently only for sparse graph and for dense graphs, simple scheme takes the lead.
Number of nodes(n)
Density(d)
Simple scheme

F heap scheme
1000
0.1 %
-
-

1%
75
67

20%
82
71

50%
89
110

75%
96
159

100%
98
184

Number of nodes(n)
Density(d)
Simple scheme
F heap scheme
3000

0.1%
1384
298

1%
1636
862

20%
1906
1894

50%
2021
2381

75%
2294
2778

100%
2564
3845
Number of nodes(n)
Density(d)
Simple scheme
F heap scheme
5000

0.1%
5204
2431

1%
6531
5803

20%
8986
8538

50%
10316
11901

75%
11164
12877

100%
12100
14611

Result Observations:
The above collected data almost meet the expectations according to the complexity of the algorithms.
As n is increased, n2 greatly increases but there is not much effect on nlogn+e where e is almost negligible therefore Fibonacci heap performs better with large n values.
As density increases so does the number of edges and hence e starts playing important part in the nlogn+e formula and for enough dense graphs, simple scheme starts performing better than Fibonacci due to the nlogn in f heap.
