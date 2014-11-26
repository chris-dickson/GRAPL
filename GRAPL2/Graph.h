#pragma once
#include <float.h>
#include <set>
#include <list>
#include <algorithm>
#include "TreeNode.h"
#include "AdjacencyList.h"
#include "IHeapInsertable.h"
#include "BinaryHeap.h"
#include "SSSP_NodeData.h"
#include "SSSP_DataCache.h"

/**********************************************************************************************//**
 * @brief	A representation of an arbitrary directed graph with floating point edge weights
 *
 * @author	Chris Dickson
 * @date	8/3/2011
 **************************************************************************************************/

class Graph
{
	public:

		Graph()
		{
			Name = NULL;
		}

		Graph(const char* _name)
		{
			SetDebugName(_name);
		}

		~Graph()
		{
			if(Name != NULL)
			{
				delete[] Name;
				Name = NULL;
			}

			// Cleans up any memory allocated while computing shortest paths
			CleanupSSSP();
		}

		void DeleteReferencedNodes()
		{
			vector<Node*> AllTheNodes;
			AllNodes(AllTheNodes);

			for(vector<Node*>::iterator NodeIt = AllTheNodes.begin(); NodeIt != AllTheNodes.end(); NodeIt++)
			{
				delete *NodeIt;
			}
		}

		void SetDebugName(const char* _name)
		{
			if ( Name != NULL )
			{
				delete[] Name;
				Name = NULL;
			}

			Name = new char[strlen(_name)+1];
			strcpy(Name,_name);
		}

		Node* InsertNode(const char* Name)
		{
			return AdjacencyData.InsertNode(Name);
		}

		Node* InsertNode(Node* Other)
		{
			return AdjacencyData.InsertNode(Other);
		}

		void RemoveNode(Node* ToRemove)
		{
			// Remove all edges ending at this node;
			multiset< TreeNode, CompareTreeEdge > AllTheEdges;
			multiset< TreeNode, CompareTreeEdge >::iterator EdgeIt;
			AllEdges(AllTheEdges);

			for (EdgeIt = AllTheEdges.begin(); EdgeIt != AllTheEdges.end(); EdgeIt++)
			{
				if ( !strcmp((*EdgeIt).n2->Name,ToRemove->Name) )
				{
					RemoveDirectedEdge((*EdgeIt).n1, (*EdgeIt).n2 );
				}
			}

			// Remove the header for this node
			AdjacencyData.RemoveNode(ToRemove);
		}

		void RemoveNode(const char* Name)
		{
			Node* ToRemove = GetNode(Name);
			RemoveNode(ToRemove);
		}

		Node* GetNode(const char* Name)
		{
			return AdjacencyData.GetNode(Name);
		}

		unsigned int NumNodes() const
		{
			return AdjacencyData.NumNodes();
		}

		unsigned int NumEdges() const
		{
			return AdjacencyData.NumEdges();
		}

		float GetEdgeWeight(  Node* from,  Node* to ) const
		{
			return AdjacencyData.GetEdgeWeight(from, to);
		}

		bool HasEdge( Node* from, Node* to ) const
		{
			return AdjacencyData.HasEdge(from,to);
		}

		void InsertDirectedEdge( Node* from,  Node* to, float weight)
		{
			assert(from != NULL && to != NULL);
			AdjacencyData.InsertEdge(from, to, weight);
		}

		void InsertUndirectedEdge( Node* n1,  Node* n2, float weight)
		{
			assert( n1 != NULL && n2 != NULL );
			InsertDirectedEdge(n1, n2, weight);
			InsertDirectedEdge(n2, n1, weight);
		}

		void RemoveDirectedEdge( Node* from,  Node* to)
		{
			assert( from != NULL && to != NULL );
			AdjacencyData.RemoveEdge(from,to);
		}

		void RemoveAllEdges()
		{
			AdjacencyData.RemoveAllEdges();
		}

		void AllEdges(multiset<TreeNode, CompareTreeEdge>& Edges);
		void AllNodes(vector<Node*>& Nodes);
		void AllLeafs(vector<Node*>& Leafs, bool bIsUndirected=false);

 		void MST(multiset<TreeNode, CompareTreeEdge>& Edges);
 		void MST();			// converts this graph to an MST
		void SteinerTree(vector<Node*>& Terminals, multiset<TreeNode, CompareTreeEdge>& Edges, bool bCacheSSSPData=false);
 		float GetTreeWeight(multiset<TreeNode, CompareTreeEdge>&Edges);

		void Print() const
		{
			return;
			if ( Name != NULL )
			{
				cout << Name << ": " << endl;
			}
			AdjacencyData.Print();
			cout << endl;
		}
		
 		void SSSP(Node* source);
 		void ExtractShortestPathTo(Node* Source, Node* Destination, vector<Node*>& PathList);
		void InvalidateSSSPCache();
		float PathLength(vector<Node*>& PathList) const;

	protected:
		virtual void PostProcessSteinerTreeGraph(Graph* ST) { }

	private:

 		void InitializeSSSP(const Node* source);
		void CleanupSSSP();
 		void Relax(SSSP_NodeData* u, SSSP_NodeData* v, BinaryHeap<float> *Q) const;

		AdjacencyList AdjacencyData;
		char* Name;
		map<Node*, SSSP_NodeData*>* NodeDataMap;
		map<Node*, map<Node*, SSSP_NodeData*>* > NodeDataCache;
};

void Graph::Relax(SSSP_NodeData* u, SSSP_NodeData* v, BinaryHeap<float> *Q) const
{
	float DMod = u->D() + AdjacencyData.GetEdgeWeight(u->TheNode, v->TheNode);
	if ( v->D() > DMod )
	{
		Q->DecreaseKey(DMod, v->HeapIndex);
		v->p = u->TheNode;
	}
}

void Graph::InitializeSSSP(const Node* source)
{
	vector<Node*> AllTheNodes;
	AllNodes(AllTheNodes);

	// Create a new map for this node
	NodeDataMap = new map<Node*, SSSP_NodeData*>();

	for (unsigned int NodeIdx = 0; NodeIdx < AllTheNodes.size(); NodeIdx++)
	{
		SSSP_NodeData *NewNodeData = new SSSP_NodeData;

		if ( AllTheNodes[NodeIdx] == source )
		{
			*(NewNodeData->key) = 0.f;
		}
		else
		{
			*(NewNodeData->key) = FLT_MAX;
		}

		NewNodeData->p = NULL;
		NewNodeData->TheNode = AllTheNodes[NodeIdx];
		NewNodeData->NodeIndex = NodeIdx;

		(*NodeDataMap)[ AllTheNodes[NodeIdx] ] = NewNodeData;
	}
}

void Graph::CleanupSSSP()
{
	// deletes the cache data
	InvalidateSSSPCache();
}

void Graph::SSSP(Node* source) 
{
	if ( NodeDataCache.find(source) != NodeDataCache.end() )
	{
		return;
	}

	InitializeSSSP(source);

	vector <SSSP_NodeData*> S;
	BinaryHeap<float> *Q = new BinaryHeap<float>();

	map<Node*, SSSP_NodeData*>::iterator NodeDataIt;
	for (NodeDataIt = NodeDataMap->begin(); NodeDataIt != NodeDataMap->end(); NodeDataIt++)
	{
		Q->Insert((IHeapInsertable<float>*)(NodeDataIt->second));
	}

	while (!Q->IsEmpty())
	{
		SSSP_NodeData *u = (SSSP_NodeData*)(Q->RemoveMin());
		S.push_back(u);

		Node* uNode = u->TheNode;
		vector<Node*> AdjacentNodesToU;
		AdjacencyData.AllNodesAdjacentTo(uNode, AdjacentNodesToU);

		SSSP_NodeData *uNodeData = (*NodeDataMap)[uNode];

		for (unsigned int i = 0; i < AdjacentNodesToU.size(); i++)
		{
			map<Node*, SSSP_NodeData*>::iterator AdjNodeIt = NodeDataMap->find(AdjacentNodesToU[i]);
		
			Relax( uNodeData, AdjNodeIt->second, Q);
		}
	}

	NodeDataCache[source] = NodeDataMap;
	delete Q;
}

void Graph::ExtractShortestPathTo(Node* Source, Node* Destination, vector<Node*>& PathList) 
{
	Node* CurrentNode = Destination;
	
	while ( CurrentNode != NULL )
	{
		PathList.insert(PathList.begin(), CurrentNode);
		map<Node*, SSSP_NodeData*>* NodeData = NodeDataCache[Source];
		//map<Node*, SSSP_NodeData*>* NodeData = NodeDataMap;
		map<Node*, SSSP_NodeData*>::iterator NodeIt;
		
		NodeIt = NodeData->find(CurrentNode);
		assert(NodeIt != NodeData->end());

		SSSP_NodeData* SSSPNodeData = NodeIt->second;

		CurrentNode = SSSPNodeData->p;
	}
}

void Graph::InvalidateSSSPCache()
{
	map< Node*, map<Node*, SSSP_NodeData*>* >::iterator CacheIt;
	map<Node*, SSSP_NodeData*>::iterator NodeDataIt;
	
	for (CacheIt = NodeDataCache.begin(); CacheIt != NodeDataCache.end(); CacheIt++)
	{
		map<Node*, SSSP_NodeData*>* NodeDataElement = CacheIt->second;

		for ( NodeDataIt = NodeDataElement->begin(); NodeDataIt != NodeDataElement->end(); NodeDataIt++ )
		{
			SSSP_NodeData* ToDelete = NodeDataIt->second;
			delete ToDelete;
		}
		delete NodeDataElement;
	}
	NodeDataCache.clear();
}

float Graph::PathLength(vector<Node*>& PathList) const
{
	vector<Node*>::iterator src, dst;
	float Length = 0.f;
	for ( src = PathList.begin(); src+1 != PathList.end(); src++)
	{
		dst = src+1;
		Length += GetEdgeWeight(*src, *dst);
	}
	return Length;
}

void Graph::AllEdges(multiset<TreeNode, CompareTreeEdge>& Edges) 
{
	Edges.clear();
	AdjacencyData.AllEdges(Edges);
}

void Graph::AllNodes(vector<Node*>& Nodes)
{
	Nodes.clear();
	AdjacencyData.AllNodes(Nodes);
}

void Graph::AllLeafs(vector<Node*>& Leafs, bool bIsUndirected) 
{
	Leafs.clear();
	AdjacencyData.AllLeafs(Leafs, bIsUndirected);
}

list<multiset<Node*, CompareNodePtr> >::iterator FindNodeInForest( list<multiset<Node*, CompareNodePtr> >& Forest, Node* n ) 
{
	list<multiset<Node*, CompareNodePtr> >::iterator ForestIt;
	multiset<Node*, CompareNodePtr>::iterator NodeIt;

	for (ForestIt = Forest.begin(); ForestIt != Forest.end(); ForestIt++)
	{
		for (NodeIt = (*ForestIt).begin(); NodeIt != (*ForestIt).end(); NodeIt++)
		{		
			if ( (*NodeIt) == n )
			{
				return ForestIt;
			}
		}
	}
	assert(false);			// Couldn't find the node anywhere.....this is bad
	return ForestIt;
}

void Graph::MST(multiset<TreeNode, CompareTreeEdge>& MSTOut) 
{ 
	list< multiset<Node*, CompareNodePtr> > Forest;
	multiset<TreeNode, CompareTreeEdge> Edges;
	multiset<Node*, CompareNodePtr>::iterator TreeIt;
	list<multiset<Node*, CompareNodePtr> >::iterator n1It, n2It;
	vector< TreeNode > HalfMst;

	vector<Node*> AllTheNodes;
	AllNodes(AllTheNodes);

	for (unsigned int i = 0; i < AllTheNodes.size(); i++)
	{
		multiset<Node*, CompareNodePtr> Tree;
		
		Tree.insert( AllTheNodes[i] );
		Forest.insert(Forest.begin(), Tree);
	}

	// Get a set of all the edges in the graph.   They will already be sorted.
	AllEdges(Edges);

	multiset<TreeNode, CompareTreeEdge>::iterator Edge;
	for( Edge = Edges.begin(); Edge != Edges.end(); Edge++ ) 
	{

		// Find the sets for n1 and n2
		n1It = FindNodeInForest(Forest, Edge->n1);
		n2It = FindNodeInForest(Forest, Edge->n2);

		if ( n1It != n2It )
		{
			multiset<Node*, CompareNodePtr>::iterator N1Begin, N1End, N2Begin, N2End;
			HalfMst.insert(HalfMst.begin(),*Edge);

			// Merge the multisets with n1 and n2
			vector<Node*> UnionSet((*n1It).size() + (*n2It).size());
			set_union((*n1It).begin(), (*n1It).end(), (*n2It).begin(), (*n2It).end(), UnionSet.begin());

			// Remove old nodes from forest
			Forest.erase(n1It);
			Forest.erase(n2It);

			// UNION Fores[n1idx] with forest[n2idx]
			multiset<Node*, CompareNodePtr> UnionMultiset;
			for (vector<Node*>::iterator UnionIt = UnionSet.begin(); UnionIt != UnionSet.end(); UnionIt++)
			{
				UnionMultiset.insert(*UnionIt);		
			}
	
			Forest.insert(Forest.begin(), UnionMultiset);
		}
	}

	// Double up the MST edges to assume undirected edges!
	for (vector<TreeNode>::iterator EdgeIt = HalfMst.begin(); EdgeIt != HalfMst.end(); EdgeIt++)
	{
		TreeNode NewEdge;
		NewEdge.n1 = (*EdgeIt).n2;
		NewEdge.n2 = (*EdgeIt).n1;
		NewEdge.weight = (*EdgeIt).weight;

		MSTOut.insert(NewEdge);
		MSTOut.insert(*EdgeIt);
	}
}

void Graph::MST()
{
	multiset<TreeNode, CompareTreeEdge> MSTEdges;
	MST(MSTEdges);

	RemoveAllEdges();

	SetDebugName("Should have no edges");
	Print();

	multiset<TreeNode, CompareTreeEdge>::iterator EdgeIt;
	for (EdgeIt = MSTEdges.begin(); EdgeIt != MSTEdges.end(); EdgeIt++)
	{
		InsertDirectedEdge( EdgeIt->n1, EdgeIt->n2, EdgeIt->weight );
	}
}

#include "DistanceNetwork.h"

void Graph::SteinerTree(vector<Node*>& Terminals, multiset<TreeNode, CompareTreeEdge>& Edges, bool bCacheSSSPData)
{
	DistanceNetwork *DN = new DistanceNetwork(this, Terminals);
	DN->SetDebugName("Distance Network");
	DN->Print();
	DN->MST();
	DN->SetDebugName("MST of Distance Network");
	DN->Print();
	Graph* NDN = DN->ReplaceWithGraphEdges(this);
	NDN->SetDebugName("Replaced with GG edges");
	NDN->Print();
	NDN->MST();
	NDN->SetDebugName("MST of GG[NDN]");
	NDN->Print();

	vector<Node*> Leafs;
	NDN->AllLeafs(Leafs,true);

	// Delete all leaves not in NodeList
	vector<Node*>::iterator LeafIt, TerminalIt;
	for (LeafIt = Leafs.begin(); LeafIt != Leafs.end(); LeafIt++)
	{
		bool bIsNonTerminalLeaf = true;
		for (TerminalIt = Terminals.begin(); TerminalIt != Terminals.end(); TerminalIt++)
		{
			if ( !strcmp( (*TerminalIt)->Name, (*LeafIt)->Name ))
			{
				bIsNonTerminalLeaf = false;
			}
		}

		if ( bIsNonTerminalLeaf )
		{
			NDN->RemoveNode((*LeafIt)->Name);
		}
	}

	if ( !bCacheSSSPData )
	{
		CleanupSSSP();
	}

	PostProcessSteinerTreeGraph(NDN);

	// NDN is now an approximate minimal steiner tree, save it out
	NDN->AllEdges(Edges);

	// destroy the intermediate graphs
	delete NDN;
	delete DN;
}

float Graph::GetTreeWeight(multiset<TreeNode, CompareTreeEdge>&Edges) 
{
	multiset<TreeNode, CompareTreeEdge>::iterator Edge;
	float Sum = 0.f;

	for ( Edge = Edges.begin(); Edge != Edges.end(); Edge++)
	{
		Sum += (*Edge).weight;
	}

	return Sum;
}
