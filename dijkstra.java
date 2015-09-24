/**
*
* @author rohit
*/

import java.io.File;
import java.util.Scanner;
import java.util.Random;
import java.util.Set;
import java.util.ArrayList;
import java.util.HashSet;
import java.io.FileNotFoundException;


public class dijkstra 
{ 
	public static final double INFINITY = Double.POSITIVE_INFINITY;
	public static int fibs = 0;
	
	public static void main(String[] args) throws FileNotFoundException
	{
		int n, i, j;
		Vertex[] edges = null;
		long time;
		Rnode[] r;
			
		if(args[0].equals("-r"))
                {
                    System.out.println("Random mode:");
                    int x = Integer.parseInt(args[1]);
                    int y = Integer.parseInt(args[2]);

                    time = System.currentTimeMillis();
                    randommode(x,y);								//calls to random mode function and passing parameters x and y.
                    time = System.currentTimeMillis() - time;
                    //System.out.println(": " + time + " ms");
                }
		else
                {
                    if(args.length == 2)							//the graph info is read from the file
                    {
                        String filename = args[1];
                        edges = load_file(filename);
                    }
                    
                    n = get_number_of_nodes(edges);
                    r = create_graph(edges, n);
                    double[][] disvector = new double[n][n];
                    
                    switch (args[0]) 
                    {
                        case "-s":
                            disvector = simplescheme(r);				//calls Simple Scheme
                            break;
                        case "-f":
                            disvector = fiboscheme(r);					//calls Fibonacci-Heap Scheme
                            break;
                        default:
                        	System.out.println("A wrong argument has been passed !");
                        	System.exit(0);
                    }

                    	for(j=0; j<n; j++)
                        {
                            if(disvector[fibs][j] == Double.POSITIVE_INFINITY)
                                System.out.print("Invalid\n");
                            else
                                System.out.print((int)disvector[fibs][j] + "\n");
                        }
                }
	}
	
		//Random Mode definition
	
	public static void randommode(int a, int b)
	{
		int n = a; 								//number of nodes in graph
		double density = (b * 0.01); 			//density of graph given by user
		Vertex[] edges;
		Vertex tempEdge;
		Rnode[] r = null;
		int numOfEdges, currentNum;
		int start, end;
		Random gen = new Random();
		boolean is_connect = false;
		long start1=0, end1, start2=0, end2, stime, ftime;
		int percentage;
		boolean sameEdge;
		System.out.println("Number of nodes\t " + "Density\t" + "Simple scheme\t " + "F-heap scheme\t ");
				is_connect = false;
				numOfEdges = (int) (n * (n-1) * density * 0.5);
				edges = new Vertex[numOfEdges];
				while(!is_connect)			//checks connection of graph, execute the three schemes only when the graph is a connected
				{
					currentNum = 0;
					while(currentNum < numOfEdges)		//generate (n*(n-1)/2)*d) edges
					{
						start = gen.nextInt(n);
						end = gen.nextInt(n);
						{
							tempEdge = new Vertex(start, end, gen.nextInt(1000) + 1);
							sameEdge = false;
							for(Vertex e : edges)
							{
								if(tempEdge.equals(e))	//checks if the new generated edge already exists iin the graph
								{	
									sameEdge = true;
									break;
								}
							}
							if(sameEdge == false)
							{
								edges[currentNum] = tempEdge;
								currentNum ++;
							}
						}
					}
					r = create_graph(edges, n);
					is_connect = is_connect(r);		//this function checks for connection in the graph
				}
				
				start1 = System.currentTimeMillis();
				simplescheme(r);
				end1 = System.currentTimeMillis();
				stime = end1 - start1;
				
				start2 = System.currentTimeMillis();
				fiboscheme(r);
				end2 = System.currentTimeMillis();
				ftime = end2 - start2;
				
				percentage = (int) (density * 100.0);		//prints the density percentage
                                
				System.out.println(n + "\t\t\t" + percentage + "%\t\t" + stime + "\t\t" + ftime);
	}
        
		//Simple Scheme definition
	
	public static double[][] simplescheme(Rnode[] r)
	{															//graph r is passed as parameter
		int n = r.length;
		double[][] disvector = new double[n][n];
		for(int i=0; i<n; i++)
			for(int j=0; j<n; j++)
				disvector[i][j] = INFINITY;
		
		int minIndex = 0;
		double minDist;
		Set<Rnode> V = new HashSet<Rnode>();	//makes a hashset of set of nodes, the shortest paths to those nodes are not yet found

		for(int source=0; source<n; source++)
		{
			disvector[source][source] = 0;
			
			V.clear();
			for(Rnode gvertex : r)
			{
				V.add(gvertex);
			}
			
			while(!V.isEmpty())		//finds shortest path to nodes
			{
				minDist = Double.POSITIVE_INFINITY;
				for(Rnode gvertex : V)	//find the node that has the shortest dist among all undetermined nodes
				{

					if(minDist > disvector[source][gvertex.getIndex()])
					{
						minDist = disvector[source][gvertex.getIndex()];
						minIndex = gvertex.getIndex();
					}
				}
				
				V.remove(r[minIndex]);

				//the dist between the source and a neighbor of the new added node may be reduced
				for(Adjacent_node an : r[minIndex].getAdjnodes())
				{
					if(disvector[source][an.getIndex()] > disvector[source][minIndex] + an.getDist())
						disvector[source][an.getIndex()] = disvector[source][minIndex] + an.getDist();
				}	
			}
		}
		return disvector;			//return shortest dist graph
	}

		//Fibonacci heap scheme
	public static double[][] fiboscheme(Rnode[] r)
	{															//graph r as parameter
		Fibo_heap fheap = new Fibo_heap();
		int n = r.length;
		double[][] disvector = new double[n][n];
		for(int i=0; i<n; i++)
			for(int j=0; j<n; j++)
				disvector[i][j] = INFINITY;
		int minIndex = 0;
		Set<Rnode> V = new HashSet<Rnode>();
	
		for(int source=0; source<n; source++)
		{
			disvector[source][source] = 0;
			V.clear();
			for(Rnode gvertex : r)
			{
				V.add(gvertex);
			}
			
			fheap.insert(source, 0);
			
			V.remove(r[source]);
			
			while(fheap.getMin() != null)
			{
				minIndex = fheap.removemin().getIndex();
				//dist to the adjacent nodes of the recently removed node may be reduced
				for(Adjacent_node an : r[minIndex].getAdjnodes())
				{
					if(disvector[source][an.getIndex()] > disvector[source][minIndex] + an.getDist())
					{
						disvector[source][an.getIndex()] = disvector[source][minIndex] + an.getDist();
						//if the node has not been inserted to the heap, then insert it to the heap
						if(V.contains(r[an.getIndex()]))
						{
							fheap.insert(an.getIndex(), disvector[source][an.getIndex()]);
							V.remove(r[an.getIndex()]);
						}
						
						//if the node is already in the heap, find the corresponding heap node, and decrease the dist
						else
							fheap.lowerKey(fheap.searchNode(an.getIndex(), fheap.getMin()), disvector[source][an.getIndex()]);
					}
				}	
			}
		}
		return disvector;			//return shortest dist graph
	}
	
		//read informatin from file
	public static Vertex[] load_file(String filename) throws FileNotFoundException
	{														//name of the file passed as parameter
		Scanner in = new Scanner(new File(filename));
		String s;
		String[] var1, var2, var3;
		StringBuffer buffer = new StringBuffer();
		String fibst;
		int x,y;
		fibst = in.nextLine();
		fibs = Integer.parseInt(fibst);			//reads the source node
		fibst = in.nextLine();

		while(in.hasNextLine())
		{
			s = in.nextLine();
			buffer.append(s);
			buffer.append("/");
		}
		
		String combination = buffer.toString();
		var1 = combination.split("/");
		y = var1.length;	
		x = (2*var1.length);
		Vertex[] edges = new Vertex[x];

		int j=0;
		for(int i=0; i<x; i =i+2)
		{
			var2 = var1[j].split(" ");
			var3 = var2;
			edges[i] = new Vertex(Integer.parseInt(var2[0]), Integer.parseInt(var2[1]), Integer.parseInt(var2[2]));
			edges[i+1] = new Vertex(Integer.parseInt(var3[1]), Integer.parseInt(var3[0]), Integer.parseInt(var3[2]));
            j=j+1;
		}
		
		in.close();
		return edges;			//returns edges
	}
	

	//compute number of nodes
	public static int get_number_of_nodes(Vertex[] edges)
	{												//edge information in class Vertex as parameter
		int n = 0;					//number of nodes;
		for(Vertex e : edges)
		{
			n = Math.max(n, e.getStart());
			n = Math.max(n, e.getEnd());
		}
		n = n + 1;					//number of nodes is equal to the maximum index of nodes + 1
		return n;
	}
	

	//create a graph using the edge information, and number of nodes
	public static Rnode[] create_graph(Vertex[] edges, int n)
	{													//edge info and number of nodes as parameters
		Rnode[] rvertex = new Rnode[n];
		for(int i=0; i<n; i++)
		{
			rvertex[i] = new Rnode(i);
		}
		
		for(Vertex e: edges)
		{
			rvertex[e.getStart()].setNeighbour(e.getEnd(), e.getDist());
		}
			
		return rvertex;				//return the created graph
	}
	

	//test for the connectivity of nodes
	public static boolean is_connect(Rnode[] r)
	{												//graph r as parameter
		double[][] disvector = simplescheme(r);		
		//execute simple scheme to check the dist between any pair of nodes

		for(int i=0; i<r.length; i++)
			for(int j=0; j<r.length; j++)
			{												//if the dist is larger than the possible maximum dist,
				if(disvector[i][i] > 1000 * r.length + 1)		//then the graph is not connected
					return false;
			}
		return true;		//return true if the graph is connected
	}
}


class Vertex 
{
  	private int start;
	private int end;
	private double dist;

	//constructor for class Vertex
	public Vertex(int start, int end, double dist)
	{										//start point, end point, dist in between as parameters
		this.start = start;
		this.end = end;
		this.dist = dist;
	}
	
	public int getStart(){
		return start;
	}
	
	public int getEnd(){
		return end;
	}
	
	public double getDist(){
		return dist;
	}
	
	public boolean equals(Vertex e){
		if(e != null && start == e.getStart() && end == e.getEnd())
			return true;
		else return false;
	}

}

class Adjacent_node 
{
  private int index;
	private double dist;
		//constructor for class edge
	public Adjacent_node(int i, double d)
	{							//node index and dist to adjacent node as parameters
		index = i;
		dist = d;
	}
	
	public int getIndex(){
		return index;
	}

	public double getDist(){
		return dist;
	}
}


	//class for fibonacci heap scheme
class Fibo_heap 
{
  private Fnode min;
	private int size;

 	public static class Fnode
	{
		private int index;
		private int degree;
		private double dist;
		private Fnode parent;
		private Fnode child;
		private Fnode leftSibling;
		private Fnode rightSibling;
		private boolean childCut;
		
		public Fnode(int index, double dist)
		{
			this.index = index;
			this.dist = dist;
			this.degree = 0;
			this.childCut = false;
		}
		
		public void setIndex(int i){
			index = i;
		}
		
		public int getIndex(){
			return index;
		}
		
		public void setDist(double d){
			dist = d;
		}
		
		public double getDistcance(){
			return dist;
		}
		
		public void setParent(Fnode p){
			parent = p;
		}
		
		public Fnode getParent(){
			return parent;
		}
		
		public void setChild(Fnode c){
			child = c;
		}
		
		public Fnode getChild(){
			return child;
		}
		
		public void setLeftSibling(Fnode ls){
			leftSibling = ls;
		}
		
		public Fnode getLeftSibling(){
			return leftSibling;
		}
		
		public void setRightSibling(Fnode rs){
			rightSibling = rs;
		}
		
		public Fnode getRightSibling(){
			return rightSibling;
		}
		
		public void increaseDgree(){
			degree ++;
		}
		
		public void decreaseDegree(){
			degree --;
		}
		
		public int getDegree(){
			return degree;
		}
		
		public void setChildCut(boolean b){
			childCut = b;
		}
		
		public boolean getChildCut(){
			return childCut;
		}
		
	}

	//constructor to create an empty Fibonacci heap
	public Fibo_heap()
	{
		min = null;
		size = 0;
	}
	
	//constructor to create a Fibonacci heap with one node
	public Fibo_heap(Fnode x){
		min = x;					//parameter x becomes the min of new heap
		x.setChild(null);
		x.setParent(null);
		x.setLeftSibling(x);
		x.setRightSibling(x);
		size = 1;
	}
	
	public Fnode getMin(){
		return min;
	}
	
	//to get size of heap
	public int getSize(){
		return size;		//return number of verices
	}
	

	//make y, child of x
	public void link(Fnode y, Fnode x) 		
	{									//x becomes y's parent, y becomes x's child				
		Fnode l = y.getLeftSibling();
		Fnode r = y.getRightSibling();
		
		l.setRightSibling(r);
		r.setLeftSibling(l);
		y.setParent(x);
		y.setChildCut(false);
		
		if(x.getDegree() == 0)				//if y becomes the only child
		{
			y.setLeftSibling(y);
			y.setRightSibling(y);
		}
		
		else								//add y to the child list of x
		{
			y.setRightSibling(x.getChild());
			y.setLeftSibling(x.getChild().getLeftSibling());
			x.getChild().getLeftSibling().setRightSibling(y);
			x.getChild().setLeftSibling(y);
		}	
		x.setChild(y);
		x.increaseDgree();
			
	}
	

	//merge of two F-heaps
	public void merge(Fibo_heap h)
	{									//another F-heap as parameter
		Fnode min2 = h.getMin();
		size += h.getSize();
		
		if(min == null && min2 !=null)
		{
			min = min2;
		}

		if(min != null && min2 != null)		//add h's root list to this heap's root list, adjust min if necessary
		{
			Fnode l = min.getLeftSibling();
			
			min2.getLeftSibling().setRightSibling(min);
			l.setRightSibling(min2);
			min.setLeftSibling(min2.getLeftSibling());
			min2.setLeftSibling(l);
			
			if(min.getDistcance() > min2.getDistcance())
				min = min2;
		}
	}
	

	//insert new node to heap
	public void insert(int index, double dist)
	{											//index of new node, value of dist of new node as parameters
		Fnode fvertex = new Fnode(index, dist);
		Fibo_heap h = new Fibo_heap(fvertex);
		this.merge(h);
	}
	

	//supporting function for removemin). links nodes having same degree in root list
	public void placenode()
	{
		int bSize = (int) (Math.log(size) / Math.log(2)) + 2; 	
		//degMap[i] means the node whose degree is i
		Fnode[] degMap = new Fnode[bSize];					
		for(int i=0; i<bSize; i++)
			degMap[i] = null;
		
		int d;
		Fnode x = min;
		Fnode start = min;
		Fnode y, next;
		do
		{
			d = x.getDegree();
			next = x.getRightSibling();
			while(degMap[d] != null)
			{
				y = degMap[d];
				if(x.getDistcance() > y.getDistcance())		//exchange x, y if necessary (refer to link())
				{
					Fnode var1 = y;
					y = x;
					x = var1;
				}
				if(y == start)								//maintain the ending mark for the loop
					start = start.getRightSibling();
				
				if(y == next)								//maintain the next pointer for x (should be in the root list)
					next = next.getRightSibling();
				this.link(y, x);
				degMap[d] = null;
				d ++;
			}
			degMap[d] = x;
			x = next;
		}while(x != start);
		
		
		min = null;
		for(int i=0; i<bSize; i++)		//adjust min 
		{
			if(min == null || (degMap[i] != null && min.getDistcance() > degMap[i].getDistcance()))
				min = degMap[i];
		}
	}
	

	//remove min from root list
	public Fnode removemin()
	{
		Fnode z = min;
		if(z != null)
		{
			Fnode x = z.getChild();		//set children's parent to null
			if(x != null)
			{
				x.setParent(null);
				for(x = x.getRightSibling(); x != z.getChild(); x = x.getRightSibling())
					x.setParent(null);
			
			
				Fnode l = z.getLeftSibling();
				x = z.getChild();
			
				l.setRightSibling(x);	//add z's children to root list
				x.getLeftSibling().setRightSibling(z);
				z.setLeftSibling(x.getLeftSibling());
				x.setLeftSibling(l);
			}
			
			z.getLeftSibling().setRightSibling(z.getRightSibling());	//remove z from root list
			z.getRightSibling().setLeftSibling(z.getLeftSibling());
					
			if(z == z.getRightSibling())		
				min = null;
				
			else
			{
				min = z.getRightSibling();
				this.placenode();		//adjustment of min is done in placenode()
			}
			size--;
		}
		return z;		//returns min node
	}


	//reduce dist in node x to 'k'
	public void lowerKey(Fnode x, double k)
	{										//a node in the heap and the new value of dist as parameters
		if(k > x.getDistcance())
		{
			System.out.println("New key is more then current key !");
			System.exit(1);
		}
		
		x.setDist(k);
		Fnode y = x.getParent();
		if(y != null && x.getDistcance() < y.getDistcance())
		{
			cut(x, y);			//move x to root list
			setCut(y);	//cascading cut if necessary
		}
		
		if(x.getDistcance() < min.getDistcance())
			min = x;			//adjust min if necessary
	}
	

	//cut node x from y, add x to root list
	public void cut(Fnode x, Fnode y)
	{										//y is x's parent
		Fnode z = y.getChild();							//remove x from the child list of y
		if (z == x && x.getRightSibling() == x)			//x is y's only child
			y.setChild(null);
		
		if (z == x && x.getRightSibling() != x)			//x is not the only child, but the pointer points to x
		{
			Fnode l = x.getLeftSibling();
			Fnode r = x.getRightSibling();
			l.setRightSibling(r);
			r.setLeftSibling(l);
			y.setChild(r);
		}
		
		if(z != x)				//x is not the only child, and the pointer does not point to x
		{
			Fnode l = x.getLeftSibling();
			Fnode r = x.getRightSibling();
			
			l.setRightSibling(r);
			r.setLeftSibling(l);
		}
		
		y.decreaseDegree();
		
		x.setRightSibling(min);						//add x to root list, as y exists, min is not null
		min.getLeftSibling().setRightSibling(x);
		x.setLeftSibling(min.getLeftSibling());
		min.setLeftSibling(x);
		
		x.setParent(null);
		x.setChildCut(false);
	}
	

	//arrange y if it is set true
	public void setCut(Fnode y)
	{								//node y as parameter which has just lost one child node
		Fnode z = y.getParent();
		if(z != null)
		{
			if(y.getChildCut() == false)
				y.setChildCut(true);
			
			else
			{
				this.cut(y, z);
				this.setCut(z);
			}
		}
	}
	

	//delete node from heap
	public void delete(Fnode x)
	{										//the node to be deleted passed as parameter
		this.lowerKey(x, Double.NEGATIVE_INFINITY);
		this.removemin();
	}
	

	//search a node in the heap whose root node is x
	public Fnode searchNode(int index, Fnode x) 
	{									//index value of node, root node of heap as parameters
		Fnode wantedNode = null;
		if(x != null)
		{
			Fnode y  = x;
			do
			{
				if(y.getIndex() == index)
					return y;
				if((wantedNode = searchNode(index, y.getChild())) != null)
					return wantedNode;
				y = y.getRightSibling();
				
			}while (y != x);
		}
		return wantedNode;				//return found node
	}
	

	//print the heap whose root node is x
	public void print(Fnode x)
	{					//root node passed as parameter
		Fnode y = x;
		if(x != null)
		{
			do
			{
				System.out.print(y.getIndex() + "\t");
				print(y.getChild());
				
				System.out.println();
				
				y = y.getRightSibling();
			}
			while(y != x);
		}
	}
	
}


class Rnode 
{
	private int index;
	private ArrayList<Adjacent_node> adjnodes;		//list of adjacent nodes
	private int neighborNumber = 0;
	
	public Rnode(int index)
	{
		this.index = index;
		adjnodes = new ArrayList<Adjacent_node>();
	}
	
	//insert an adjacent node
	public void setNeighbour(int i, double d)
	{						//index and dist of adjacent node as parameters
		adjnodes.add(new Adjacent_node(i, d));
	}
	

	//get the adjacent node
	public ArrayList<Adjacent_node> getAdjnodes()
	{
		return adjnodes;			//adjacent node list
	}
	
	public int getIndex()
	{
		return index;
	}
	
	public int getNeighbourNumber()
	{
		return neighborNumber;
	}
}