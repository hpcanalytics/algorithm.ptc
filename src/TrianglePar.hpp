#ifndef TRIANGLEPAR_H_
#define TRIANGLEPAR_H_

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <time.h>
#include <mpi.h>

//#include "SGraph.hpp"
#include "utility.hpp"
#include "Set.hpp"
#define LINESIZE 4000
#define ERROR -1

using namespace std;

typedef  int ULI;
typedef  int  Vertex;
typedef  int  Vsize;



class ClusterCo
{
private:
				
	ULI     N;				// number of nodes
	double  avgcc;

	Set     *nlist, *nlist2;
		
	int *count, *count2, *selfcount;
	int *deg, *deg2;
	int numProc, me;
	Set nArray, coreArray;                              // alternative of STL map
	ULI start, end;
	ULI csize, esize;
	long int tcount;
	
public:
	
	ClusterCo();
	~ClusterCo();
	void initPar(int npro, int self);		
	void NodeCC(Vertex v);
	void AllNodeCC();	
	double WriteCC();
	double AvgCC(); 
	void ReadGraph(const char *fname);
	void OrderAndStore();
	void catOut(const char *oFileName);
	void CommData1();
	void CommData2(); 
	ULI GetN(); 
	double ReduceFiles(string fname[], int FileCount, string oFileName);
	
	long int getTcount();
};	


	
ClusterCo::ClusterCo()	
{ 	
	N=0;
	start= 0;
	end=0;
						
	avgcc = 0;
	nlist=NULL;
	nlist2=NULL;
	count=NULL;
	count2=NULL;
	deg=NULL;
	deg2= NULL;
	selfcount=NULL;
	tcount=0;
		
}


// Destructor

ClusterCo::~ClusterCo()
{ 
	
	if (deg) delete [] deg;
	if (deg2) delete [] deg2;
	nArray.destroy();
	coreArray.destroy();
	
}

void ClusterCo::initPar(int npro, int self)
{
	numProc=npro;
	me=self;

}

double ClusterCo::AvgCC() {
	 
	return avgcc;
	
}

ULI ClusterCo::GetN()
{
	return N;	
}

long int ClusterCo::getTcount() {
	 
	return tcount;
	
}

// compute cluster coefficient of node v in graph g

void ClusterCo::NodeCC(Vertex v)
{
 	long int j = nArray.member_index(v)-1;;
 	
 	Vsize  degree = nlist2[j].Size();
	
	if (degree <= 1) 
		return ;
	
	Vertex u;
	long int pos;
	int k;
	Set A(degree);
	
	for (Vsize i=0; i<degree; i++) {
		
		u = nlist2[j][i];
					
		pos = nArray.member_index(u)-1;   // binary search
		
		A.intersect2(0, nlist2[j], nlist2[pos]);  
	
		k= A.Size();
		
		tcount = tcount + k ;
	
	}
	
	A.destroy();
		
}


void ClusterCo::AllNodeCC()
{
	if(start>=N)
		return;
	

	Vertex v;
	
	//for (Vsize i=0; i < csize; i++)
		//count[i]=0;
	
	///for (Vsize i=0; i < esize; i++)
		///count2[i]=0;
	
	/*	
	for (v = start; v < end; v++) {
		
		NodeCC(v);
			
	}
	
	*/
	
	int i;
	for (i = 0; i < csize; i++) {
		
		NodeCC(coreArray[i]);
			
	}
				
	if(nlist)
		for (ULI i=0; i < csize;i++)
			nlist[i].destroy();
			
	if(nlist2)
		for (ULI i=0; i < esize;i++)
			nlist2[i].destroy();
					
	if (nlist) delete []nlist;
	if (nlist2) delete []nlist2;
	//cout << "P"<<me<<" :testing: after allnodecc"<<endl;
	
}

void ClusterCo::CommData1()    // MPI built-in functinality
{
	int quota, *recCount;
	quota= (int) ceil((double) N/numProc);
	int i;
	
	int qProc= (int) floor((double) N/quota);
	
	// prepare recCount Array
	
	recCount= new int[N];        // won't work if number of element per quota exceeds 2 billion 
	
	for (i= 0; i < qProc; i++)
	{
		recCount[i]= quota;	
	}
	
	if (numProc > qProc)
		recCount[qProc] = N - (qProc*quota);
	
	for (i= qProc+1; i < numProc; i++)
	{
		recCount[i]= 0;	
	}
	
	// Prepare receive buffer
	
	if (recCount[me] > 0)
		selfcount = new int[recCount[me]]; 
	
	else
		selfcount= new int [1];
	
	MPI_Reduce_scatter( count, selfcount, recCount, MPI_INT, MPI_SUM, MPI_COMM_WORLD );	
	
	if (count) delete [] count;
	
	
}

void ClusterCo::CommData2()    // explicit comm
{
	char filename[512];
	ofstream *occ;
	occ=new ofstream[numProc];
	int NumChar;
		
	for (int j=0;j<numProc;j++){
		
		NumChar=sprintf(filename,"./OUT/out_%d_%d.txt",me,j);
	
		occ[j].open(filename);

		if(!occ[j].is_open()){
			cout << "intermediate file creation failed." << endl;
			exit(1);
		}
	}
	
	int quota = (int) ceil((double) N/numProc);
	
		
	for (long int i=0; i < esize; i++){
		
		int p = nArray[i] % numProc; 
		occ[p] << nArray[i] << "\t" << count2[i] << "\t" << deg2[i] << endl;	

	}
	
	
	for (int j=0;j<numProc;j++)
			occ[j].close();
	if (count) delete [] count;
	if (count2) delete [] count2;	
	
	//cout << "P"<<me << " :: testing: after commdata2"<<endl;	
		
}


double ClusterCo::ReduceFiles(string fname[], int FileCount, string oFileName)
{
	
	struct record { ULI vertex; int count; int degree;}; 
	typedef struct record Record;
	
	ifstream *filePointer;  
	filePointer=new ifstream[FileCount]; 

	ofstream oFile;
	oFile.open(oFileName.c_str());
	
	if(!oFile.is_open()){
		cout<<" Intermediate output file creation failed."<<endl;
		exit (1);

	}

	ULI j=0;
	
	bool *flagArray=new bool[FileCount];
	int i, recdeg, sumFlag=FileCount;
	Record *recArray;
	recArray= new Record[FileCount];  
 
	ULI vertexU, tempCount;
	
	int min, first=1;

	for ( i=0;i<FileCount;i++)
	{
	
		filePointer[i].open(fname[i].c_str());

		if(!filePointer[i].is_open()){
			cout<<fname[i]<<" not found: "<<i<<endl;
			
			exit(1);
		}
		
		if (filePointer[i]>>vertexU>>tempCount>>recdeg){
			flagArray[i]=true;   //to keep track number of <vertex, count> pair still remaining for file i
			recArray[i].vertex=vertexU;
			recArray[i].count=tempCount;
			recArray[i].degree=recdeg;
			
			if (first) {min = vertexU; first=0; }
			
			else if( vertexU < min) min= vertexU; 
		}
		else
		{
			flagArray[i]= false;
			sumFlag--;
		}
				

	}


	//ULI itr = start;
	
	while(sumFlag){

		ULI newCount=0, minvDeg;
		
					
		for(i=0;i<FileCount;i++)
		{
		
			if(flagArray[i]&&(recArray[i].vertex == min)){
				
				newCount=newCount+recArray[i].count;
				minvDeg= recArray[i].degree;             // can be improved later on; redundant assignment
				
				if (filePointer[i]>>vertexU>>tempCount>>recdeg){
					recArray[i].vertex=vertexU;
					recArray[i].count=tempCount;
					recArray[i].degree=recdeg;
					
				}
				else{
					flagArray[i]= false;	
					sumFlag--;
				}
			}

	
		}
		
		double cc=0;
		
		//cout<<"count: "<<itr<< "\t" <<newCount<<"\t"<<deg[itr-start]<<endl;
		if(!newCount){
			oFile<< min <<"\t"<< newCount <<endl;
			
		}	
		else{
			
			cc= (double) (newCount * 2) / (double)(minvDeg * (minvDeg-1));
			oFile << min <<"\t"<< cc << endl;
	    
		}   
	
		avgcc += cc; 
			
		// calculate new min
		
		first=1;

		for ( i=0;i<FileCount;i++)
		{
					
			if (flagArray[i]){
				
				vertexU = recArray[i].vertex;
				
				if (first) {min = vertexU; first=0; }
				
				else if( vertexU < min) min= vertexU; 
			}
							
	
		}
	
		
					

	}  //while


	for (i=0; i<FileCount;i++)
		if(filePointer[i]!=NULL)
			filePointer[i].close();

	
	oFile.close();

	delete []flagArray;
	delete []recArray;
	delete []filePointer;
	//cout << "P"<<me << " :: testing: after reduction"<<endl;
	return avgcc;
}


double ClusterCo::WriteCC()
{
	
	int quota;
	ULI start, end;
	quota= (int) ceil((double) N/numProc);
	start= me*quota;
	if(start>=N)
		return 0;
		
	if (start+ quota < N)
		end = start +quota;
	else 
		end = N;
	
	char filename[512];
	int NumChar=sprintf(filename,"./OUT/out%d",me);
	
	ofstream outfile(filename, ios::out);
	
	double cc=0;
	ULI j=0;
	
	for(ULI i=start; i< end; i++){
	
		if(!selfcount[j]){
			cc= selfcount[j++];
			outfile<< i <<"\t\t"<< cc <<endl;
		}	
		else{
		    cc= (double) selfcount[j++]*2 / (deg[i] * (deg[i]-1));
		    outfile<< i <<"\t\t"<< cc <<endl;
		    
		}   
		
		avgcc += cc; 
	}
	
	//avgcc = avgcc / N;   //temp commented
		
	outfile.close();
	if (selfcount) delete [] selfcount;
	
	//cout<<"Average CC: "<< avgcc <<endl;  //temp commented
	
	return avgcc;
}



void ClusterCo::catOut(const char *oFileName)
{
	time_t startMerge=time(NULL);
	cout<<"Concatenating files..."<<endl;
	
	ofstream oFile;
	oFile.open(oFileName);
	if(!oFile.is_open())
	{
		cout<<"Output file creation failed!"<<endl;
		exit(1);
	}
	
	
	char CatCommand[1024];
	
	for (int i=0; i< numProc; i++)
	{
		sprintf(CatCommand, "cat ./OUT/merge%d.txt >> %s", i, oFileName);
		system(CatCommand);	
		
	}
	
	oFile.close();
	cout<<"\n Final output file created."<<endl;
	time_t mergeTime=time(NULL)-startMerge;
	cout<<" MergeTime: "<<mergeTime<<" Seconds"<<endl;

}



void ClusterCo::ReadGraph(const char *fname) 
{
	ULI     u, v, temp;					// vertices
	Vsize   i, j, k, tmpdeg, dummy;
						
	double  wt;
	
	ifstream ifp(fname);
	
	if(!ifp.is_open()) {
		cout << "\% Cannot open " << fname << endl;
		exit(1);	
	}

	//cout << "\% Reading graph file: " << fname << " ... "<<endl; //cout.flush();
	
	
	ifp >> N;									// read first line -- the number of vertex
	
	ULI quota;
	
	quota= (int) ceil((double) N/numProc);
	start= me*quota;
	if(start>=N)
	{	
		csize=esize=0;
		ifp.close();
		return;
	}
	
		
	
	
	
	if (start+ quota < N)
		end = start +quota;
	else 
		end = N;
	
	
	//deg= new int[N];
			
	// Skip 'start' number of entries
	
	for (i=0;i<start;i++){
		ifp >> u >> tmpdeg;
		//deg[i]=tmpdeg;
		for(k=0; k<tmpdeg; k++){
      		ifp >> v >> wt >> temp;
      		//ifp >> v ;
		}
	}
	
	csize= end - start; 
    
    nlist = new Set[csize];
    deg= new int[csize];
    count= new int[csize];
	
	nArray.init(csize);
	coreArray.init(csize);
	
	// enlist neighbors that need to further explore
			
	for (i = start; i < end; i++) {
		
		j = i-start;
		
    	ifp >> u >> tmpdeg;						// read nodes and its degree
		
		nArray.Dinsert(u);
		coreArray.Dinsert(u);
		
		//--deg[j] = tmpdeg;
		//--nlist[j].init(tmpdeg);
						
		for(k=0; k<tmpdeg; k++) {
      		
      		ifp >> v >> wt >> temp;					// read the adjacent nodes
			//ifp >> v ;
			
			//--nlist[j].insert(v);	
			
			//--if (!(v >= start && v < end))	
			nArray.Dinsert(v);
			
			
   		}
   		/*
   		if(nlist[j].Size()){
				nlist[j].finalize2();
				nlist[j].sort2();
		}
		*/
	}
	
	
	if (nArray.Size()){
		//nArray.finalize2();  //nArray finalize
		nArray.sort2();      //sort ascending
	}
	
	long int ii=0,kk=0;
	long int asz= nArray.Size();
	
	//while (ii < asz && nArray[ii] < end ) ii++;
	
	while(ii < asz)  // duplicate elimination
	{
		
		nArray[kk]=nArray[ii++];
		
		while (ii < asz && nArray[ii]==nArray[kk])	ii++;
		
		kk++;
	}   
	
	nArray.size = kk;
	esize = kk;
	
	nArray.finalize2();
	
	if (esize){
		
		nlist2 = new Set[esize];
	    count2 = new int[esize];
	    deg2= new int [esize];
		
	}
	
	
	
	// exploring further for elements in neighbor list	
	//2nd pass to the file
	
	//--if(me){
	
		ifp.clear();              // forget we hit the end of file
		ifp.seekg(0, ios::beg);   // move to the start of the file
		ifp >> N ;   // again scan N ?
	//--}
	
		
	if (esize){
						
		// now scan for neighbors and degree info
		
		ULI ia=0;
		
		while (ifp.good() && ia < esize)
		{
			ifp >> u >> tmpdeg;	
			
			
			if (nArray[ia]==u){                     // so the assumption is, strictly increasing order
				
					
				nlist2[ia].init(tmpdeg);
				deg2[ia] = tmpdeg;
		 		
				for(k=0; k<tmpdeg; k++) {
		      		
		      		ifp >> v >> wt >> temp;					// read the adjacent nodes
								
					//if (v > u)
					nlist2[ia].insert(v);		// add the edge (u, v) to the graph
		   		}
				
			
				
				ia++;	
			}
			
			else 
			{
				for(k=0; k<tmpdeg; k++) {
		      		
		      		ifp >> v >> wt >> temp;					// read the adjacent nodes
								
				}
				
			}
				
		}
	}	
	
	//cout << "P"<<me << " :: testing: after read"<<endl;
	
	ifp.close();
	
	OrderAndStore();
	
	
}



void ClusterCo::OrderAndStore()
{
	//clock_t start= clock();
	
	int dg; 
		
	int slot=0, idx;
	
	for (ULI i=0; i < esize; i++) {
	
    	dg= deg2[i];
		
		slot=0;
		
		for(int k=0; k < dg; k++) {
      		      		
			idx = nArray.member_index(nlist2[i][k])-1;
			if (idx >= 0 )        // then not in outside world
				if ((deg2[idx] > dg) || (deg2[idx]==dg && nlist2[i][k] > nArray[i]))
					nlist2[i][slot++]= nlist2[i][k];
					      												
   		}
   		
   		nlist2[i].size=slot;
   		if (nlist2[i].Size()){
			nlist2[i].finalize2();    // to release excess memory
			nlist2[i].sort2();
		}
		
	}	
	
	//if (me ==0)
		//cout << " P0: Time needed for reorganizing memory: "<< (clock()-start) / (double) CLOCKS_PER_SEC << " Sec" << endl;
	//else if (me == numProc-1)
		//cout << " P_last: Time needed for reorganizing memory: "<< (clock()-start) / (double) CLOCKS_PER_SEC << " Sec" << endl;	
	
}



#endif /* #ifndef TRIANGLEPAR_H_ */


