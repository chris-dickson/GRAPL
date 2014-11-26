#include "Node.h"

/**********************************************************************************************//**
 * @brief	Adjacency list element.  Stores the adjacent node plus the weight of the node.
 *
 * @author	Chris Dickson
 * @date	8/3/2011
 **************************************************************************************************/

class AdjacencyListElement
{
	public:
		AdjacencyListElement(){}

		AdjacencyListElement( Node* NodeIn, float WeightIn)
		{
			AdjacentNode = NodeIn;
			Weight = WeightIn;
		}

		~AdjacencyListElement(){}

		const char* GetName() const { return AdjacentNode->Name; }
		float GetWeight() const { return Weight; }
		void SetWeight(float newWeight) { Weight = newWeight; }
		Node* GetAdjacentNode() { return AdjacentNode; }
		
	private:
		Node	*AdjacentNode;
		float	Weight;
};