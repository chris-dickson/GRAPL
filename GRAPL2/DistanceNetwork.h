#pragma once
#include <map>
#include <vector>
#include "Graph.h"
#include "TreeNode.h"

class DistanceNetwork : public Graph
{
public:
	DistanceNetwork( Graph* G, vector<Node*> NodeList );
	Graph* ReplaceWithGraphEdges(Graph* G);
private:
	vector<vector<Node*> > ShortestPathLists;
};



DistanceNetwork::DistanceNetwork( Graph* G, vector<Node*> NodeList )
{
	vector<Node*>::iterator NodeIt;
	// Insert nodes in NodeList
	for (NodeIt = NodeList.begin(); NodeIt != NodeList.end(); NodeIt++)
	{
		InsertNode(*NodeIt);
	}

	vector<Node*> AllTheNodes;
	AllNodes(AllTheNodes);

	// compute shortest path between all sets of points in NodeList
	for (NodeIt = AllTheNodes.begin(); NodeIt != AllTheNodes.end(); NodeIt++)
	{
		vector<Node*>::iterator dst;
		G->SSSP(G->GetNode((*NodeIt)->Name));

		for (dst = AllTheNodes.begin(); dst != AllTheNodes.end(); dst++)
		{
			if ( NodeIt != dst )
			{
				vector<Node*> ShortestPath;
				G->ExtractShortestPathTo(G->GetNode((*NodeIt)->Name), G->GetNode((*dst)->Name), ShortestPath);
				
				ShortestPathLists.insert(ShortestPathLists.begin(), ShortestPath);

				float PathWeight = G->PathLength(ShortestPath);

				InsertUndirectedEdge(*NodeIt, *dst, PathWeight);
			}
		}
	}

}


Graph* DistanceNetwork::ReplaceWithGraphEdges(Graph* G) 
{
	Graph* NewGraph = new Graph();
	
	multiset<TreeNode, CompareTreeEdge> AllTheEdges;
	multiset<TreeNode, CompareTreeEdge>::iterator EdgeIt;
	AllEdges(AllTheEdges);

	for (EdgeIt = AllTheEdges.begin(); EdgeIt != AllTheEdges.end(); EdgeIt++)
	{
		Node *SourceNode, *DestNode;
		SourceNode = EdgeIt->n1;
		DestNode   = EdgeIt->n2;

		// Find the shortest paths with source/dest pair
		vector< vector<Node*> >::iterator PathListIt;
		for (PathListIt = ShortestPathLists.begin(); PathListIt != ShortestPathLists.end(); PathListIt++)
		{
			unsigned int PathListLength = (*PathListIt).size();
			if ( (!strcmp((*PathListIt)[0]->Name,SourceNode->Name)) && (!strcmp((*PathListIt)[PathListLength-1]->Name,DestNode->Name)) ) 
			{
				break;
			}
		}

		assert( PathListIt != ShortestPathLists.end() );

		// Add the path edges to the new graph
		vector<Node*>::const_iterator n = (*PathListIt).begin();
		while( n != (*PathListIt).end() )
		{
			if ( n + 1 != (*PathListIt).end() )
			{
				Node* src = NewGraph->InsertNode(*n);
				Node* dst = NewGraph->InsertNode(*(n + 1));
				NewGraph->InsertUndirectedEdge(src, dst, G->GetEdgeWeight(src, dst));
			}
			n++;
		}

	}

	return NewGraph;
}