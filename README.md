GR-APL
======
A global routing application for VLSI design	

					
To compile:
-----------
There are two versions of the routing code, each located within its own directory.  The multi-threaded version is located in the MT/ directory, 
while the path-saving, single-threaded version is located in the MEM/ directory.  To compile either, cd to the directory and type:

	make
	
After compilation, there will be four executable files located in the MT/build/ and five executable files in the MEM/build directory:

	MT/route_MT : approximately solves LP relaxation to the ILP global routing formulation (multi-threaded, no path-saving)
	MEM/route   : approximately solves LP relaxation to the ILP global routing formulation (single-threaded, path-saving)
	MEM/route_MT: approximately solves LP relaxation to the ILP global routing formulation (multi-threaded, path-saving)
	frac_round_MRD : performs randomized rounding to solution produced by route_MT
	prune_yal : format YAL benchmark files
	prune_ibm : format IBM benchmark files
	
To execute:
-----------
Once compiled, the program should be executed as follows:

	route_MT bench.txt args.txt

where bench.txt is the text file containing the instance of the problem and args.txt is the argument file. 
For a sample argument file, see ARGS_EX.txt.

To round the solution, type:

	frac_round_MRD trees.txt stats.txt sol.txt [time]
	
where trees.txt are the trees computed by route_MT, stats.txt is the statistics file produced by route_MT, 
sol.txt is the integer solution file to be written, and an optional integer representing the number of 
minutes spent in rounding.   If no time is specified, rounding will continue until the program is halted 
manually.

The files prune_yal and prune_ibm are used to format the input files to route_MT.  These programs remove
unnecessary nets from the instance that do not need to be dealt with in the global routing phase (repeated
terminal nets, 1 terminal nets).  Example:

	prune_ibm in.txt out.txt [tol]
	
where in.txt is the original benchmark, out.txt is the formatted benchmark and the optional integer tol is
the maximum number of terminals allowed in a net.   (i.e./ tol of 9 will only write nets with less than 10
terminals to the output file)

Visualization MATLAB Extension
------------------------------
A visualization MATLAB extension has been provided to produce animations of the algorithm in action.  The code
will produce files with the name Congestion_xx.txt in the build/output/ directory.  To produce the animations, 
copy these files to the VIZ/ directory, open MATLAB, cd to the VIZ/ directory and type:

	grid_movie('Congestion', N, 'movie.avi', gw, gh, fps);

where 'Congestion' is the prefix of the files in the output/ directory, N is the number of iterations performed by
route or route_MT, gw and gh represent the width and height of the grid graph respectively, and fps is the number of
frames per second for the file movie.avi. 

List of Files:
--------------

	binheap.h 	 : header file for a binary heap min-priority queue
	binheap.c 	 : implementation for a binary heap min-priority queue
	fatal.h		 : error handling
	frac_round_MRD.c : main file for rounding program
	frac_round_aux.c : auxilliary functions for frac_round_MRD
	globals.MT.c	 : initialization routines for global variables used in main.MT.c
	globals.c	 : initialization routines for global variables use in main.c
	graph.h		 : header file for graphs
	graph.c 	 : implementation of graphs
	graphtype.h	 : header file for structured graphs (virtual layer, grid)
	graphtype.c	 : intialization routines for structured graphs
	hash.h		 : header file for hash tables (unused)
	hash.c 		 : implementation of hash tables (unused)
	helper.h	 : auxilliary routines
	helper.c	 : auxilliary routines
	linesearch.h	 : header file for line search and potential function evaluations
	linesearch.c	 : implementation of line search and potential function evaluations
	list.h		 : header file for linked lists
	list.c		 : implementation of linked lists
	main.MT.c	 : implementation of MEM/route_MT
	main.c		 : implementation of route
	main.ecmm.c	 : implementation of edge congestion minimization model route_ECMM (unused)
	mst.h		 : header file for computing minimum spanning trees
	mst.c		 : implementation of prims algorithm for computing MST's
	parseargs.h	 : header file for parsing command line arguments (symbol table)
	parseargs.c	 : implementation of parser
	prune_ibm.c	 : main file for prune_ibm
	prune_yal.c	 : main file for prune_yal
	rounding.h	 : header file for randomized rounding
	rounding.c	 : implementation of randomized rounding
	routing.h	 : header file for routing specific code (nets, parsing benchmarks, etc)
	routing.c	 : implementation of routing specific code
	set.h		 : header file for ADT-Set
	set.c		 : implementation of ADT-Set
	shortest.MT.h	 : header file for computing shortest paths in MEM/route_MT
	shortest.MT.c	 : implementation of Dijkstra's algorithm for computing shortest paths in MEM/route_MT
	shortest.h	 : header file for computing shortest paths in route
	shortest.c	 : implementation of Dijkstra's algorithm for computing shortests paths in route
	steiner.MT.h	 : header file for computing Steiner trees in route_MT
	steiner.MT.c	 : implementation for computing Steiner trees in route_MT
	steiner.h	 : header file for computing Steiner trees in route
	steiner.c	 : implementation for computing Steiner trees in route
 

References:
-----------

C. Dickson, Global Routing in VLSI: Algorithms, Theory, and Computation, Masters Thesis, Department of Computing
and Software, McMaster University, 2007.  

-CD 06/2007
