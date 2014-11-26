#pragma once

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

class ParseMCNCBenchmark : public IParseBenchmark
{
public:

	virtual int ParseFile(const char* FileName, VLSIProblem& ProblemOut);

};

int ParseMCNCBenchmark::ParseFile(const char* FileName, VLSIProblem& ProblemOut)
{
	ProblemOut.Nets.clear();

	string line;
	ifstream myfile (FileName);
	if (myfile.is_open())
	{
		vector< string > SplitResults;

		// num nets
		unsigned int NumNets;
		getline(myfile, line);
		SplitResults = SplitString(line, ' ');
		NumNets = (unsigned int)(atoi(SplitResults[3].data()));

		// num modules (unused)
		getline(myfile,line);

		// grid height
		getline(myfile,line);
		SplitResults = SplitString(line,' ');
		ProblemOut.Height = (unsigned int)(atoi(SplitResults[3].data()));

		// grid width
		getline(myfile,line);
		SplitResults = SplitString(line,' ');
		ProblemOut.Width = (unsigned int)(atoi(SplitResults[3].data()));

		// capacity
		getline(myfile,line);
		SplitResults = SplitString(line,' ');
		ProblemOut.VerticalCapacity = (unsigned int)(atoi(SplitResults[3].data()));
		ProblemOut.HorizontalCapacity = ProblemOut.VerticalCapacity;

		// Horizontal length
		getline(myfile,line);
		SplitResults = SplitString(line, ' ');
		ProblemOut.HorizonalEdgeLength = (float)(atof(SplitResults[3].data()));

		// Vertical length
		getline(myfile,line);
		SplitResults = SplitString(line, ' ');
		ProblemOut.VerticalEdgeLength = (float)(atof(SplitResults[3].data()));

		// Parse the two empty lines
		getline(myfile, line);
		getline(myfile, line);

		// parse the nets
		for (unsigned int i = 0; i < NumNets; i++)
		{
			vector<string> NetInformationResults;
			vector<string> TerminalPartsResults;

			Net NewNet;
			getline(myfile, line);
			SplitResults = SplitString(line, ':');

			// Parse the net information
			string NetInformation = SplitResults[0];
 			NetInformationResults = SplitString(NetInformation,',', true);
 			NewNet.Name = NetInformationResults[0];

			TerminalPartsResults = SplitString(SplitResults[1], ' ');

			for (unsigned int j = 1; j < TerminalPartsResults.size() - 1; j++)
			{
				vector<string> TerminalCoordinates;
				UINTVector2 Terminal;
				TerminalCoordinates = SplitString(TerminalPartsResults[j], ',');

				Terminal.x = (unsigned int)(atoi(TerminalCoordinates[0].data())) - 1;
				Terminal.y = (unsigned int)(atoi(TerminalCoordinates[1].data())) - 1;

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