#pragma once
#include <vector>
#include <iostream>
#include "Node.h"
#include "Graph.h"
using namespace::std;

/**********************************************************************************************//**
 * @brief	Adjacency list header.  This stores the node itself, plus a list of adjacent nodes
 * 			with each weight.
 *
 * @author	Chris Dickson
 * @date	8/3/2011
 **************************************************************************************************/

class AdjacencyListHeader
{
	public:

		/**********************************************************************************************//**
		 * @brief	Constructor.   Takes the name of the node this header is for as input
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 *
		 * @param	NameIn	The name of the node
		 **************************************************************************************************/

		AdjacencyListHeader(const char* NameIn)
		{
			TheNode = new Node(NameIn);
		}

		AdjacencyListHeader(Node* Other)
		{
			TheNode = Other;
		}


		/**********************************************************************************************//**
		 * @brief	Finaliser.   
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 **************************************************************************************************/

		~AdjacencyListHeader()
		{
			RemoveAllEdges();
		}

		/**********************************************************************************************//**
		 * @brief	Inserts an edge.   If the edge exists, we reset the weight provided.
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 *
		 * @param	ToNode	to node.
		 * @param	weight	The weight.
		 **************************************************************************************************/

		void InsertEdge( Node *ToNode, float weight)
		{
			// reweight it if it exists already
			for (unsigned int i = 0; i < AdjacentVerticies.size(); i++)
			{
				if ( !strcmp(AdjacentVerticies[i]->GetName(),ToNode->Name) )
				{
					AdjacentVerticies[i]->SetWeight(weight);
					return;
				}
			}

			AdjacentVerticies.insert( AdjacentVerticies.end(), 1, new AdjacencyListElement( ToNode, weight ));
		}

		void RemoveEdge( Node* ToNode )
		{
			vector< AdjacencyListElement* >::iterator AVIt;

			for (AVIt = AdjacentVerticies.begin(); AVIt != AdjacentVerticies.end(); AVIt++)
			{
				if ( !strcmp((*AVIt)->GetName(), ToNode->Name) )
				{
					AdjacentVerticies.erase(AVIt);
					return;
				}
			}
			cout << "WARNING:   Tried to remove an edge from " << TheNode->Name << " to " << ToNode->Name << " that doesn't exist!" << endl;
		}

		void RemoveAllEdges()
		{
			// Delete all the edges (AdjacencyListElements)
			vector<AdjacencyListElement*>::iterator AdjVertIt;
			for ( AdjVertIt = AdjacentVerticies.begin(); AdjVertIt != AdjacentVerticies.end(); AdjVertIt++)
			{
				AdjacencyListElement *ToDelete = *AdjVertIt;
				delete ToDelete;
			}
			AdjacentVerticies.clear();
		}

		/**********************************************************************************************//**
		 * @brief	Gets an edge weight.
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 *
		 * @param	ToNode	to node.
		 *
		 * @return	The edge weight, -1 if it doesn't exist.
		 **************************************************************************************************/

		float GetEdgeWeight( const Node& ToNode )
		{
			for (unsigned int i = 0; i < AdjacentVerticies.size(); i++)
			{
				if ( AdjacentVerticies[i]->GetName() == ToNode.Name )
				{
					return AdjacentVerticies[i]->GetWeight();
				}
			}
			cout << "Warning!   Edge from " << TheNode->Name << " to " << ToNode.Name << " doesn't exist!" << endl;
			return -1;
		}

		bool HasEdge( const Node& ToNode )
		{
			for (unsigned int i = 0; i < AdjacentVerticies.size(); i++)
			{
				if ( AdjacentVerticies[i]->GetName() == ToNode.Name )
				{
					return true;
				}
			}
			return false;
		}

		unsigned int NumEdges() const
		{
			return AdjacentVerticies.size();
		}

		void AddToEdgeList(multiset<TreeNode, CompareTreeEdge>& Edges)
		{
			for (unsigned int i = 0; i < NumEdges(); i++)
			{
				TreeNode newEdge;
				newEdge.n1 = TheNode;
				newEdge.n2 = AdjacentVerticies[i]->GetAdjacentNode();
				newEdge.weight = AdjacentVerticies[i]->GetWeight();
				Edges.insert(newEdge);
			}
			
		}

		void AllAdjacentNodes(vector<Node*>& NodeList)
		{
			for(unsigned int i = 0; i < NumEdges(); i++)
			{
				NodeList.push_back(AdjacentVerticies[i]->GetAdjacentNode());
			}
		}


		/**********************************************************************************************//**
		 * @brief	Prints this object.
		 *
		 * @author	Chris Dickson
		 * @date	8/3/2011
		 **************************************************************************************************/

		void Print() const
		{
			for (unsigned int i = 0; i < AdjacentVerticies.size(); i++)
			{
				cout << "[" << AdjacentVerticies[i]->GetName() << "]:" << AdjacentVerticies[i]->GetWeight();
				if ( i < AdjacentVerticies.size() - 1 )
				{
					cout << " -> ";
				}
			}
		}

		Node* TheNode;
		vector< AdjacencyListElement* > AdjacentVerticies;
};