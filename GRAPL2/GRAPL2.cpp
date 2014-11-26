// GRAPL2.cpp : Defines the entry point for the console application.
//
#include "ParseIBMBenchmark.h"
#include "ParseMCNCBenchmark.h"
#include "VLSIProblem.h"
#include "VirtualLayerGraph.h"
#include "GridGraph.h"
#include "Timer.h"
#include <iostream>
//#include <windows.h>

using namespace std;

void MemoryGrind()
{
	long int Itr = 0;
	while(true)
	{
		cout << "Iteration " << Itr << endl;
		Graph* G = new Graph();
		for (int i = 0; i < 100; i++)
		{
			char NodeNameString[64];
			sprintf( NodeNameString, "%7d", i );
			Node* n1 = G->InsertNode(NodeNameString);

			sprintf( NodeNameString, "%7d", i );
			Node *n2 = G->InsertNode(NodeNameString);

			G->InsertUndirectedEdge(n1,n2,1.f);
		}
		G->DeleteReferencedNodes();
		delete G;
		Itr++;
	}
}



#define NUM_THREADS 4

#define DEBUG_PRINTS 1

struct SteinerThreadData {
	VLSIProblem* TheProblem;		// The data about the problem
	unsigned int StartNet;			// The first net this thread will generate
	unsigned int EndNet;			// The last net this thread will generate
	unsigned int ThreadIndex;		
};

/*
DWORD WINAPI SteinerThread( LPVOID lpParam )
{
	SteinerThreadData *LocalThreadData = (SteinerThreadData*)lpParam;
	VLSIProblem Problem = *(LocalThreadData->TheProblem);

	// Create VL-Graph
#if DEBUG_PRINTS
	cout << "Thread " << LocalThreadData->ThreadIndex << " initializing graph..." << endl;
#endif
	VirtualLayerGraph *G = new VirtualLayerGraph(Problem.Width, Problem.Height, Problem.HorizonalEdgeLength, Problem.VerticalEdgeLength, 1.f);

	for (unsigned int i = LocalThreadData->StartNet; i <= LocalThreadData->EndNet; i++)
	{
#if DEBUG_PRINTS
		cout << LocalThreadData->ThreadIndex << ": " << "Generating net " << i << endl;
#endif

		vector<Node*> Terminals;
		multiset<TreeNode, CompareTreeEdge> STEdges;
		set< UINTVector2, CompareUINTVector2 >::iterator TerminalIt;
		set<UINTVector2, CompareUINTVector2>& TerminalSet = Problem.Nets[i].Terminals;

		for (TerminalIt = TerminalSet.begin(); TerminalIt != TerminalSet.end(); TerminalIt++)
		{
			Node* TerminalNode = G->GetNode( (*TerminalIt).x, (*TerminalIt).y, 0 );

			Terminals.push_back(TerminalNode);
		}

		G->SteinerTree(Terminals, STEdges);
		STEdges.empty();
	}

	return 0;
}

void ThreadedMain()
{
	ParseMCNCBenchmark IBMDummy;
	VLSIProblem Problem;
	SteinerThreadData *ThreadData = new SteinerThreadData[NUM_THREADS];
	DWORD   dwThreadIdArray[NUM_THREADS];
	HANDLE  hThreadArray[NUM_THREADS]; 

	GRAPLTimer Timer;
	Time TotalTime;


	IBMDummy.ParseFile("C:\\Documents and Settings\\Chris Dickson\\Desktop\\MCNC\\ind3.yal.dat", Problem);

	// Split out the nets evenly amongst all the threads.   Give the last thread the leftovers
	int NetsPerThread = Problem.Nets.size() / NUM_THREADS;
	int StartNet = 0;
	for (int i = 0; i < NUM_THREADS; i++)
	{
		ThreadData[i].ThreadIndex = (unsigned int) i;
		ThreadData[i].TheProblem = &Problem;
		ThreadData[i].StartNet = StartNet;
		ThreadData[i].EndNet = StartNet + NetsPerThread - 1;
		StartNet += NetsPerThread;
	}
	ThreadData[ NUM_THREADS - 1 ].EndNet += Problem.Nets.size() % NUM_THREADS;


	// Create the threads
	Timer.Start();
	for( int i=0; i<NUM_THREADS; i++ )
	{
		// Create the thread to begin execution on its own.

		hThreadArray[i] = CreateThread( 
			NULL,                   // default security attributes
			0,                      // use default stack size  
			SteinerThread,       // thread function name
			&ThreadData[i],          // argument to thread function 
			0,                      // use default creation flags 
			&dwThreadIdArray[i]);   // returns the thread identifier 

	} // End of main thread creation loop.

	// Wait until all threads have terminated.

	WaitForMultipleObjects(NUM_THREADS, hThreadArray, TRUE, INFINITE);

	Timer.Stop();
	TotalTime = Timer.GetElapsed();
	cout << "Total Time: " << TotalTime.TotalSeconds() << endl;

}
*/
void SteinerGrind()
{
	ParseIBMBenchmark IBM01;
	VLSIProblem Problem;

	long int itr = 0;

	IBM01.ParseFile("ibmDummy.txt", Problem);
	
	while (1)
	{
		cout << "Iteration " << itr << endl;
		// Create VL-Graph
		VirtualLayerGraph *G = new VirtualLayerGraph(Problem.Width, Problem.Height, Problem.HorizonalEdgeLength, Problem.VerticalEdgeLength, 1.f);

		for (unsigned int i = 0; i < Problem.Nets.size(); i++)
		{
			vector<Node*> Terminals;
			multiset<TreeNode, CompareTreeEdge> STEdges;
			set< UINTVector2, CompareUINTVector2 >::iterator TerminalIt;
			set<UINTVector2, CompareUINTVector2>& TerminalSet = Problem.Nets[i].Terminals;

			for (TerminalIt = TerminalSet.begin(); TerminalIt != TerminalSet.end(); TerminalIt++)
			{
				Node* TerminalNode = G->GetNode( (*TerminalIt).x, (*TerminalIt).y, 0 );

				Terminals.push_back(TerminalNode);
			}

			G->SteinerTree(Terminals, STEdges);
		}
		
		G->DeleteReferencedNodes();
		delete G;		
		
		itr++;
	}

}


int main(int argc, char* argv[])
{
	SteinerGrind();
	return 0;
	


	ParseIBMBenchmark IBM01;
	ParseMCNCBenchmark Bio;
	VLSIProblem Problem;
	GRAPLTimer MyTimer;
	Time FileParsingTime, GridGenTime, SteinerTreeTime;


	MyTimer.Start();
	IBM01.ParseFile("C:\\Documents and Settings\\Chris Dickson\\Desktop\\FromOptlab\\benchmarks\\ibmDummy.txt", Problem);
	//Bio.ParseFile("C:\\Documents and Settings\\Chris Dickson\\Desktop\\MCNC\\bio.yal.dat", Problem);
	MyTimer.Stop();
	FileParsingTime = MyTimer.GetElapsed();
	MyTimer.Reset();

	// Create VL-Graph
	MyTimer.Start();
 	VirtualLayerGraph *G = new VirtualLayerGraph(Problem.Width, Problem.Height, Problem.HorizonalEdgeLength, Problem.VerticalEdgeLength, 1.f);
	MyTimer.Stop();
	GridGenTime = MyTimer.GetElapsed();
	MyTimer.Reset();


	cout << "Generating Nets:" << endl;
	MyTimer.Start();
	for (unsigned int i = 0; i < Problem.Nets.size(); i++)
	{
		cout << "\t" << i << " of " << Problem.Nets.size() << endl;
		vector<Node*> Terminals;
		multiset<TreeNode, CompareTreeEdge> STEdges;
		set< UINTVector2, CompareUINTVector2 >::iterator TerminalIt;
		set<UINTVector2, CompareUINTVector2>& TerminalSet = Problem.Nets[i].Terminals;

		for (TerminalIt = TerminalSet.begin(); TerminalIt != TerminalSet.end(); TerminalIt++)
		{
			Node* TerminalNode = G->GetNode( (*TerminalIt).x, (*TerminalIt).y, 0 );

			Terminals.push_back(TerminalNode);
		}

		G->SteinerTree(Terminals, STEdges);
		STEdges.empty();
	}
	MyTimer.Stop();
	SteinerTreeTime = MyTimer.GetElapsed();
	MyTimer.Reset();


	cout << "FileParsingTime: " << FileParsingTime.TotalSeconds() << endl;
	cout << "GridGenTime: " << GridGenTime.TotalSeconds() << endl;
	cout << "SteinerTreeTime: " << SteinerTreeTime.TotalSeconds() << endl;

	return 0;
}
