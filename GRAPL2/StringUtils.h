#pragma once
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <vector>

using namespace::std;

inline std::string trim_right(const std::string &source , const std::string& t = " ")
{
	std::string str = source;
	return str.erase( str.find_last_not_of(t) + 1);
}

inline std::string trim_left( const std::string& source, const std::string& t = " ")
{
	std::string str = source;
	return str.erase(0 , source.find_first_not_of(t) );
}

inline std::string trim(const std::string& source, const std::string& t = " ")
{
	std::string str = source;
	return trim_left( trim_right( str , t) , t );
} 

// TODO: reference additional headers your program requires here
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems, bool bTrimWhitespace) {
	std::stringstream ss(s);
	std::string item;
	while(std::getline(ss, item, delim)) {

		if ( bTrimWhitespace )
		{
			item = trim(item);
		}

		elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> SplitString(const std::string &s, char delim, bool bTrimWhitespace=false) {
	std::vector<std::string> elems;
	return split(s, delim, elems, bTrimWhitespace);
}

