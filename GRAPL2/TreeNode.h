#pragma once
#include "Node.h"
class TreeNode
{
public:
	Node *n1, *n2;
	float weight;

	bool operator<(const TreeNode& other) const
	{
		if ( weight < other.weight )
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

class CompareTreeEdge
{
public:
	bool operator()(TreeNode s1, TreeNode s2)
	{
		if(s1.weight < s2.weight )
			return true;
		else
			return false;
	}
};