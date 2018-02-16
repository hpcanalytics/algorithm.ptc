/*
 *	Parallel TC
 *  Assumption: in input file, nodes are listed in ascending order (galib graph format)
 *  Output: Total number of triangles
 *  Author: Arif (sm10.vbi.vt.edu); created: June 2012
 */

#include "TrianglePar.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>



using namespace std;


// pass parameters to actual triangle counting object

void ComputeCC(const char *gphfile_name, int numProc, int me)
{
	
	ClusterCo cc;   // Triangle counting object
	
	cc.initPar(numProc, me);
	
	clock_t start, absstart;
	
	absstart=start= clock();
				
	
	cc.ReadGraph(gphfile_name);
	
	time_t end;
	
		
	cc.AllNodeCC();
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (me==0)
		start=clock();
	
	long int tcount=cc.getTcount(), sumt=0;
	MPI_Reduce( &tcount, &sumt, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD );
	
	if (me==0){
		
		cout << endl<< "Total triangles in the graph: "<< sumt << endl;
		//cout<<"Time taken for reduce operation: "<< (clock()-start)/(double) CLOCKS_PER_SEC<<" Sec."<<endl;
		cout<<"Total time elapsed: "<< (clock()-absstart)/(double) CLOCKS_PER_SEC<<" Sec."<<endl;
					
	}
	 
	
}

// main function of the program

int main(int argc, char **argv)
{
	int numProc, me;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	
	if (argc != 2) {
		cout << "Input file name required. Usage: pccnod <graph_file_name>" << endl;
		return 0;
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);	
	ComputeCC(argv[1], numProc, me);
	
	
	MPI_Finalize();
	
	return 0;
}
