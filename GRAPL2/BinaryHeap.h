#pragma once
#include "Heap.h"
#include <vector>
#include <assert.h>

using namespace::std;


template<typename KeyType> class BinaryHeap : public Heap<KeyType>
{
  public:
	  BinaryHeap<KeyType>() { }
	  ~BinaryHeap<KeyType>() { }

	  virtual void Insert(IHeapInsertable<KeyType>* Element)
	  {
		  HeapData.push_back(Element);
		  Element->OnIndexChanged(HeapData.size() - 1);
		  SiftUp(HeapData.size() - 1);
	  }

	  virtual bool IsEmpty() const
	  {
		  return HeapData.size() == 0;
	  }

	  virtual IHeapInsertable<KeyType>* RemoveMin()
	  {
		  IHeapInsertable<KeyType>* toReturn = NULL;
		  if ( HeapData.size() == 0 )
		  {
			  assert(false);
		  }
		  else
		  {
			  toReturn = HeapData[0];
			  HeapData[0] = HeapData[ HeapData.size() - 1];
			  HeapData[0]->OnIndexChanged(0);
			  HeapData.pop_back();

			  SiftDown(0);
		  }
		  return toReturn;
	  }

	  virtual void DecreaseKey(KeyType NewKey, unsigned int HeapIndex )
	  {
		  ModifyKey(NewKey, HeapIndex);
		  SiftUp(HeapIndex);
	  }

	  virtual void IncreaseKey(KeyType NewKey, unsigned int HeapIndex)
	  {
		  ModifyKey(NewKey, HeapIndex);
		  SiftDown(HeapIndex);
	  }

private:
	  void SwapNodes(unsigned int parentIndex, unsigned int nodeIndex)
	  {
		  IHeapInsertable<KeyType>* tmp;

		  tmp = HeapData[parentIndex];
		  HeapData[parentIndex] = HeapData[nodeIndex];
		  HeapData[nodeIndex]->OnIndexChanged(parentIndex);

		  HeapData[nodeIndex] = tmp;
		  HeapData[nodeIndex]->OnIndexChanged(nodeIndex);
	  }

	  void SiftUp(unsigned int nodeIndex)
	  {
		  unsigned int parentIndex;
		  if (nodeIndex != 0) 
		  {
			  parentIndex = Parent(nodeIndex);
			  if (*HeapData[parentIndex] > *HeapData[nodeIndex]) 
			  {
				  SwapNodes(parentIndex, nodeIndex);
				  SiftUp(parentIndex);
			  }
		  }
	  }

	  void SiftDown(unsigned int nodeIndex) {

		  unsigned int leftChildIndex, rightChildIndex, minIndex;

		  leftChildIndex = Left(nodeIndex);
		  rightChildIndex = Right(nodeIndex);

		  if (rightChildIndex >= HeapData.size()) 
		  {
			  if (leftChildIndex >= HeapData.size())
			  {
				  return;
			  }
			  else
			  {
				  minIndex = leftChildIndex;
			  }
		  } 
		  else 
		  {
			  if (*HeapData[leftChildIndex] <= *HeapData[rightChildIndex])
			  {
				  minIndex = leftChildIndex;
			  }
			  else
			  {
				  minIndex = rightChildIndex;
			  }

		  }

		  if (*HeapData[nodeIndex] > *HeapData[minIndex]) 
		  {
			  SwapNodes(nodeIndex, minIndex);
			  SiftDown(minIndex);
		  }
	  }

	  void ModifyKey(KeyType NewKey, unsigned int HeapIndex)
	  {
		  assert(HeapIndex >= 0 && HeapIndex < HeapData.size() );
		  HeapData[HeapIndex]->ModifyKey(NewKey);
	  }
	  inline unsigned int Left(unsigned int i) { return 2 * i + 1; }
	  inline unsigned int Right(unsigned int i) { return 2 * i + 2; }
	  inline unsigned int Parent(unsigned int i) { return (i-1) / 2; }

	  vector<IHeapInsertable<KeyType>*> HeapData;
};