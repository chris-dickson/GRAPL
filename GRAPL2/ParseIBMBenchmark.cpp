#include "ParseIBMBenchmark.h"
#include "VLSIProblem.h"
#include "IParseBenchmark.h"
#include "StringUtils.h"
#include "Net.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

using namespace::std;

int ParseIBMBenchmark::ParseFile(const char* FileName, VLSIProblem& ProblemOut)
{
	ProblemOut.Nets.clear();

	string line;
	ifstream myfile (FileName);
	if (myfile.is_open())
	{
		vector< string > SplitResults;

		// grid dimensions
		getline(myfile,line);
		SplitResults = SplitString(line,' ');
		ProblemOut.Width = (unsigned int)(atoi(SplitResults[1].data()));
		ProblemOut.Height = (unsigned int)(atoi(SplitResults[2].data()));

		// vertical capacity
		getline(myfile,line);
		SplitResults = SplitString(line,' ');
		ProblemOut.VerticalCapacity = (unsigned int)(atoi(SplitResults[2].data()));

		// horizontal capacity
		getline(myfile,line);
		SplitResults = SplitString(line,' ');
		ProblemOut.HorizontalCapacity = (unsigned int)(atoi(SplitResults[2].data()));

		// Even spacing for IBM benchmarks
		ProblemOut.HorizonalEdgeLength = 1.f;
		ProblemOut.VerticalEdgeLength = 1.f;

		// num nets
		unsigned int NumNets;
		getline(myfile, line);
		SplitResults = SplitString(line, ' ');
		NumNets = (unsigned int)(atoi(SplitResults[2].data()));

		// parse the nets
		for (unsigned int i = 0; i < NumNets; i++)
		{
			Net NewNet;
			unsigned int NumTerminals;
			getline(myfile, line);
			SplitResults = SplitString(line, ' ');
			NumTerminals = (unsigned int)(atoi(SplitResults[2].data()));
			NewNet.Name = SplitResults[0];

			for (unsigned int j = 0; j < NumTerminals; j++)
			{
				UINTVector2 Terminal;

				getline(myfile, line);
				SplitResults = SplitString(line, ' ');
				Terminal.x = (unsigned int)(atoi(SplitResults[2].data()));
				Terminal.y = (unsigned int)(atoi(SplitResults[3].data()));

				NewNet.Terminals.insert(Terminal);
			}

			if ( NewNet.Terminals.size() > 1 )
			{
				ProblemOut.Nets.insert(ProblemOut.Nets.end(), NewNet);
			}
		}

		myfile.close();
	}
	else
	{
		cout << "Failed to open file " << FileName << endl;
		return 1;
	}

	return 0;
}