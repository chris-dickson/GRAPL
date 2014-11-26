#pragma once

#include "UINTVector2.h"
#include <vector>
#include <set>
#include "Net.h"

class VLSIProblem
{
public:
	unsigned int Width;
	unsigned int Height;
	unsigned int HorizontalCapacity;
	unsigned int VerticalCapacity;
	float HorizonalEdgeLength;
	float VerticalEdgeLength;

	std::vector< Net > Nets;
};