#pragma once
#include "Graph.h"
#include "DistanceNetwork.h"
#include "StringUtils.h"
#include <string>
#include <math.h>

class VirtualLayerGraph : public Graph
{
public:
	VirtualLayerGraph(unsigned int M, unsigned int N, float HorizontalEdgeWeight=1.f, float VerticalEdgeWeight=1.f, float ViaEdgeWeight=1.f);
	~VirtualLayerGraph()
	{
		// Delete the node lookups, but not the nodes themselves
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N; j++)
			{
				delete[] NodeGrid[i][j];
			}
			delete [] NodeGrid[i];
		}
		delete [] NodeGrid;
	}

	Node* GetNode(unsigned int _M, unsigned int _N, unsigned int _Z)
	{
		assert( _M >= 0 && _M < M );
		assert( _N >= 0 && _N < N );
		assert(_Z == 0 || _Z == 1);

		return NodeGrid[_M][_N][_Z];
	}

	Node* GetNode(const char* _Name)
	{
		unsigned int i,j,k;
		GetNodeCoordinates(_Name,i,j,k);
		return NodeGrid[i][j][k];
	}
protected:
	virtual void PostProcessSteinerTreeGraph(Graph* ST);

	void GetNodeCoordinates(const char* Name, unsigned int& i, unsigned int &j, unsigned int &k);
	void GetNodeCoordinates(Node* TheNode, unsigned int& i, unsigned int& j, unsigned int &k);


private:
	Node ****NodeGrid;
	unsigned int M,N;
	
};


VirtualLayerGraph::VirtualLayerGraph(unsigned int _M, unsigned int _N, float HorizontalEdgeWeight, float VerticalEdgeWeight, float ViaEdgeWeight)
{
	assert(_M > 0 && _N > 0);
	M = _M; 
	N = _N;

	// Insert all the nodes
	NodeGrid = new Node***[M];
	for (unsigned int i = 0; i < M; i++)
	{
		NodeGrid[i] = new Node **[N];

		for (unsigned int j = 0; j < N; j++)
		{
			NodeGrid[i][j] = new Node *[2];

			for (unsigned int k = 0; k < 2; k++)
			{
				unsigned int NumDigits_i = (unsigned int)(floor( log10( float(abs( (int)i )) ) )) + 1; 
				unsigned int NumDigits_j = (unsigned int)(floor( log10( float(abs( (int)j )) ) )) + 1;  
				unsigned int NumDigits_k = 1;
				char* NodeName = new char[NumDigits_i + NumDigits_j + NumDigits_k + 3];
				sprintf(NodeName,"%d,%d,%d",i,j,k);
				Node* n = new Node(NodeName);
				NodeGrid[i][j][k] = InsertNode(n);
				delete[] NodeName;
			}
		}
	}

	// Insert all horizontal edges
	for (unsigned int i = 0; i < M; i++)
	{
		for(unsigned int j = 0; j < N; j++)
		{
			// left
			if ( i > 0 )
			{
				InsertDirectedEdge( NodeGrid[i][j][0], NodeGrid[i-1][j][0], HorizontalEdgeWeight );
			}

			// right
			if ( i < (M - 1) )
			{
				InsertDirectedEdge( NodeGrid[i][j][0], NodeGrid[i+1][j][0], HorizontalEdgeWeight );
			}
		}
	}

	// Insert all vertical edges
	for (unsigned int i = 0; i < M; i++)
	{
		for(unsigned int j = 0; j < N; j++)
		{
			// up
			if ( j < (N - 1) )
			{
				InsertDirectedEdge( NodeGrid[i][j][1], NodeGrid[i][j+1][1], VerticalEdgeWeight);
			}

			// down
			if ( j > 0 )
			{
				InsertDirectedEdge( NodeGrid[i][j][1], NodeGrid[i][j-1][1], VerticalEdgeWeight);
			}
		}
	}

	// Insert all the vias
	for (unsigned int i = 0; i < M; i++)
	{
		for (unsigned int j = 0; j < N; j++)
		{
			InsertUndirectedEdge( NodeGrid[i][j][0], NodeGrid[i][j][1], ViaEdgeWeight );
		}
	}
}

void VirtualLayerGraph::GetNodeCoordinates(const char* Name, unsigned int& i, unsigned int &j, unsigned int &k)
{
	string TheNodeSTLString(Name);
	string CommaString(",");
	vector<string> Results;

	Results = SplitString(TheNodeSTLString, ',');

	assert(Results.size() == 3);
	i = (unsigned int)(atoi(Results[0].data()));
	j = (unsigned int)(atoi(Results[1].data()));
	k = (unsigned int)(atoi(Results[2].data()));
}

void VirtualLayerGraph::GetNodeCoordinates(Node* TheNode, unsigned int& i, unsigned int& j, unsigned int &k)
{
	GetNodeCoordinates(TheNode->Name,i,j,k);
}

void VirtualLayerGraph::PostProcessSteinerTreeGraph(Graph* ST)
{
	// check all the leafs (terminals by this point) and see if the have an edge connecting to them that is a via
	// If so, remove it from the graph
	vector<Node*> AllTheLeafs;
	vector<Node*>::iterator LeafIt;
	ST->AllLeafs(AllTheLeafs,true);
	
	for (LeafIt = AllTheLeafs.begin(); LeafIt != AllTheLeafs.end(); LeafIt++)
	{
		char ToNodeName[256];
		unsigned int i,j,k;
		// Check and see if this leaf has a via
		Node* From = *LeafIt;
		Node* To;

		GetNodeCoordinates(From, i,j,k);
		
		if ( k == 0 )
		{
			sprintf(ToNodeName,"%d,%d,1",i,j);
		}
		else
		{
			sprintf(ToNodeName,"%d,%d,0",i,j);
		}
		To = ST->GetNode(ToNodeName);


		if ( To != NULL && ST->HasEdge(From,To) )
		{
			ST->RemoveDirectedEdge(From,To);
			ST->RemoveDirectedEdge(To,From);
		}
	}
}
