#pragma once
#include "Node.h"
#include "TreeNode.h"
#include "AdjacenyListElement.h"
#include "AdjacencyListHeader.h"
#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
using namespace std;

/**********************************************************************************************//**
 * @brief	Adjacency list representation of a graph. 
 *
 * @author	Chris Dickson
 * @date	8/3/2011
 **************************************************************************************************/

class AdjacencyList
{
	public:
		AdjacencyList(){ }

		~AdjacencyList() 
		{
			map<Node*, AdjacencyListHeader*>::iterator HeaderIt;

			for(HeaderIt = HeaderMap.begin(); HeaderIt != HeaderMap.end(); HeaderIt++)
			{
				AdjacencyListHeader* ToDelete = HeaderIt->second;
				delete ToDelete;
			}
			HeaderMap.clear();
		}

		/**********************************************************************************************//**
		 * @brief	Inserts a node described by Name.   Checks to see if it exists already.  
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 *
		 * @param	Name	The name of the node to insert.
		 *
		 * @return	null if it fails, else a pointer to the node.
		 **************************************************************************************************/

		Node* InsertNode(const char* Name)
		{
			map<Node*, AdjacencyListHeader*>::iterator MapIt;
			for (MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				if ( !strcmp(MapIt->first->Name, Name) )
				{
					break;
				}
			}
			if ( MapIt != HeaderMap.end() )
			{
				return MapIt->first;
			}

			AdjacencyListHeader* newHeader = new AdjacencyListHeader(Name);
			HeaderMap[newHeader->TheNode] = newHeader;
			return newHeader->TheNode;
		}

		Node* InsertNode(Node* Other)
		{
			if ( HeaderMap.find(Other) != HeaderMap.end() )
			{
				return Other;
			}

			AdjacencyListHeader* newHeader = new AdjacencyListHeader(Other);
			HeaderMap[Other] = newHeader;
			return Other;
		}

		Node* GetNode(const char* name)
		{
			map<Node*, AdjacencyListHeader*>::iterator MapIt;
			for(MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				if ( !strcmp(MapIt->first->Name, name) )
				{
					return MapIt->first;
				}
			}
			return NULL;
		}

		void RemoveNode(Node* ToRemove)
		{
			// Assumes that all the edges ending at this node have been 
			// removed already!
			
			map<Node*, AdjacencyListHeader* >::iterator ToDelete;
			ToDelete = HeaderMap.find(ToRemove);
			if ( ToDelete != HeaderMap.end() )
			{
				HeaderMap.erase(ToDelete);
			}
			else
			{
				cout << "Warning!   Tried to delete node " << ToRemove->Name << " which doesn't exist!" << endl;
			}
		}

		unsigned int NumNodes() const
		{
			return HeaderMap.size();
		}

		unsigned int NumEdges() const
		{
			unsigned int Count = 0;
			AdjacencyListHeader* AdjHeader;
			map<Node*, AdjacencyListHeader*>::const_iterator MapIt;
			
			for (MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				AdjHeader = (*MapIt).second;
				Count += AdjHeader->NumEdges();
			}
			return Count;
		}

		void AllLeafs(vector<Node*>& Leafs, bool bIsUndirected) 
		{	
			unsigned int EdgeCountForLeaf = bIsUndirected ? 1 : 0;

			map<Node*, AdjacencyListHeader* >::iterator MapIt;
			for (MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				if ((MapIt->second)->NumEdges() == EdgeCountForLeaf)
				{
					Leafs.push_back( (MapIt->first) );
				}
			}
		}

		void AllNodes(vector<Node*>& Nodes)
		{
			map<Node*, AdjacencyListHeader*>::iterator MapIt;
			for (MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				Nodes.push_back( MapIt->first );
			}
		}

		void AllNodesAdjacentTo(Node* TheNode, vector<Node*>& AdjacentNodes)
		{
			map<Node*, AdjacencyListHeader*>::iterator MapIt;
			MapIt = HeaderMap.find(TheNode);

			assert(MapIt != HeaderMap.end());

			AdjacencyListHeader *TheNodeHeader = MapIt->second;
			TheNodeHeader->AllAdjacentNodes(AdjacentNodes);
		}

		void AllEdges(multiset<TreeNode, CompareTreeEdge>& Edges) 
		{
			map<Node*, AdjacencyListHeader*>::iterator MapIt;
			for (MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				Node* from = MapIt->first;
				AdjacencyListHeader *fromHeader = MapIt->second;
				fromHeader->AddToEdgeList(Edges);
			}
		}


		/**********************************************************************************************//**
		 * @brief	Gets an edge weight.
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 *
		 * @param	from	Source node for the edge.
		 * @param	to  	Destination node for the edge
		 *
		 * @return	The edge weight, -1 if it doesn't exist.
		 **************************************************************************************************/

		float GetEdgeWeight( Node* from,  Node* to) const
		{
			map<Node*,AdjacencyListHeader*>::const_iterator MapIt = HeaderMap.find(from);
			assert(MapIt != HeaderMap.end());

			AdjacencyListHeader *FromHeader = MapIt->second;

			return FromHeader->GetEdgeWeight(*to);
		}

		bool HasEdge( Node* from, Node* to) const
		{
			map<Node*,AdjacencyListHeader*>::const_iterator MapIt = HeaderMap.find(from);
			assert(MapIt != HeaderMap.end());

			AdjacencyListHeader *FromHeader = MapIt->second;

			return FromHeader->HasEdge(*to);
		}

		/**********************************************************************************************//**
		 * @brief	Inserts an edge.   Reweights it if it exists already.
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 *
		 * @param	from  	Source node for the.
		 * @param	to	  	to.
		 * @param	weight	The weight.
		 **************************************************************************************************/

		void InsertEdge(  Node* from,  Node* to, float weight )
		{
			map<Node*,AdjacencyListHeader*>::iterator MapIt = HeaderMap.find(from);
			assert(MapIt != HeaderMap.end());

			AdjacencyListHeader *FromHeader = MapIt->second;
			
			FromHeader->InsertEdge( to, weight );
		}

		void RemoveEdge(  Node* from,  Node* to)
		{
			map<Node*,AdjacencyListHeader*>::iterator MapIt = HeaderMap.find(from);
			assert(MapIt != HeaderMap.end());

			AdjacencyListHeader *FromHeader = MapIt->second;
			FromHeader->RemoveEdge( to );
		}

		void RemoveAllEdges()
		{
			map<Node*, AdjacencyListHeader*>::iterator MapIt;
			for ( MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				(MapIt->second)->RemoveAllEdges();
			}

		}

		void Print() const
		{
			map<Node*, AdjacencyListHeader*>::const_iterator MapIt;

			for (MapIt = HeaderMap.begin(); MapIt != HeaderMap.end(); MapIt++)
			{
				cout << (*MapIt).first->Name << " -> ";
				(*MapIt).second->Print();
				cout << endl;
			}
		}

	private:

		map< Node*, AdjacencyListHeader* > HeaderMap;
};