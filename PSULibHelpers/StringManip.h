#pragma once

#include <string>
#include <vector>

using namespace std;

class StringManip
{
public:
	static string trim(string string_to_trim);
	static vector<string> stringToVector(string input, string delim);
		
};

