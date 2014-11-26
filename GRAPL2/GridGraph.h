#pragma once
#include "Graph.h"
#include "DistanceNetwork.h"

#include <math.h>

class GridGraph : public Graph
{
public:
	GridGraph(unsigned int M, unsigned int N, float HoritzontalEdgeLength=1.f, float VerticalEdgeLength=1.f);
	~GridGraph(){}

	Node *GetNode(unsigned int _M, unsigned int _N)
	{
		assert( _M >= 0 && _M < M );
		assert( _N >= 0 && _N < N );
		return NodeGrid[_M][_N];
	}

	unsigned int Width() const { return M; }
	unsigned int Height() const { return N; }


private:
	Node ***NodeGrid;
	unsigned int M, N;
};


GridGraph::GridGraph(unsigned int _M, unsigned int _N, float HorizontalEdgeLength, float VerticalEdgeLength)
{
	assert(_M > 0 && _N > 0);
	M = _M; 
	N = _N;

	NodeGrid = new Node**[M];

	// Insert all the verticies of the grid
	for (unsigned int i = 0; i < M; i++)
	{
		NodeGrid[i] = new Node*[N];
		for (unsigned int j = 0; j < N; j++)
		{
			unsigned int NumDigits_i = (unsigned int)(floor( log10( float(abs( (int)i )) ) )) + 1; 
			unsigned int NumDigits_j = (unsigned int)(floor( log10( float(abs( (int)j )) ) )) + 1;  
			char* NodeName = new char[NumDigits_i + NumDigits_j + 2];
			sprintf(NodeName,"%d,%d",i,j);
			Node* n = new Node(NodeName);
			NodeGrid[i][j] = InsertNode(n);
			
		}
	}

	// Insert the edges of the grid
	for (unsigned int i = 0; i < M; i++)
	{
		for (unsigned int j = 0; j < N; j++)
		{
			// up
			if ( j < (N - 1) )
			{
				InsertDirectedEdge( NodeGrid[i][j], NodeGrid[i][j+1], VerticalEdgeLength);
			}

			// down
			if ( j > 0 )
			{
				InsertDirectedEdge( NodeGrid[i][j], NodeGrid[i][j-1], VerticalEdgeLength );
			}

			// left
			if ( i > 0 )
			{
				InsertDirectedEdge( NodeGrid[i][j], NodeGrid[i-1][j], HorizontalEdgeLength );
			}

			// right
			if ( i < (M - 1) )
			{
				InsertDirectedEdge( NodeGrid[i][j], NodeGrid[i+1][j], HorizontalEdgeLength );
			}
		}
	}
}