#pragma once
#include "IHeapInsertable.h"
#include "Node.h"

class SSSP_NodeData : public IHeapInsertable<float>
{


protected:
	virtual void OnIndexChanged(unsigned int NewIndex)
	{
		HeapIndex = NewIndex;
	}

	virtual void ModifyKey(float NewKey)
	{
		*key = NewKey;
	}

public:

	unsigned int HeapIndex;

	float D() const { return *key; }

	Node* p;
	Node* TheNode;
	unsigned int NodeIndex;
};