#pragma  once

#include "IHeapInsertable.h"

template<typename KeyType> class Heap 
{
public:
	Heap<KeyType>() {}
	~Heap<KeyType>() {}

	virtual void Insert(IHeapInsertable<KeyType>* Element) = 0;

	virtual bool IsEmpty() const = 0;

	virtual IHeapInsertable<KeyType>* RemoveMin() = 0;

	virtual void DecreaseKey(KeyType NewKey, unsigned int HeapIndex ) = 0;

	virtual void IncreaseKey(KeyType NewKey, unsigned int HeapIndex) = 0;
};