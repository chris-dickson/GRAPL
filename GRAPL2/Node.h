#pragma once

#include <string.h>

/**********************************************************************************************//**
 * @brief	A representation for a node used in a graph.   Labeled by a string.
 *
 * @author	Chris Dickson
 * @date	8/3/2011
 **************************************************************************************************/

class Node
{
  public:
	Node() 
	{ 
		Name = ""; 
	}

	Node(const char* NameIn) 
	{ 
		Name = new char[strlen(NameIn)+1];
		strcpy(Name ,NameIn); 
	}

	~Node() 
	{
		delete[] Name;
	}

	char* Name;
};

class CompareNodePtr
{
public:
	bool operator()(Node *n1, Node *n2)
	{
		return n1 == n2;
	}
};
