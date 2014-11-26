#pragma once
#include "VLSIProblem.h"

class IParseBenchmark
{
public:
	IParseBenchmark(){}

	virtual int ParseFile(const char* FileName, VLSIProblem& ProblemOut) = 0;
};