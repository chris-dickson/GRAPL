#pragma once
#include <string>
#include "UINTVector2.h"
class Net
{
public:
	std::string Name;
	std::set<UINTVector2, CompareUINTVector2> Terminals;
};