#include "StringManip.h"
#include <xlocbuf>
#include <codecvt>



using namespace std;

string StringManip::trim(string string_to_trim)
{
	string string_trimmed = string_to_trim;
	size_t startpos = string_trimmed.find_first_not_of(" ");
	size_t endpos = string_trimmed.find_last_not_of(" ");
	if (startpos == string::npos)
		string_trimmed = "";
	else if (endpos == string::npos)
		string_trimmed = string_trimmed.substr(startpos);
	else
		string_trimmed = string_trimmed.substr(startpos, endpos - startpos + 1);
	return string_trimmed;
}

string StringManip::removeChar(string string_to_trim, string char_remove)
{
	string string_trimmed = string_to_trim;
	int startpos = string_trimmed.find(char_remove);
	while (startpos > -1)
	{
		string_trimmed = string_trimmed.replace(startpos, 1, "");
		startpos = string_trimmed.find(char_remove);
	}	
	return string_trimmed;
}

vector<string> StringManip::stringToVector(string input, string delim)
{
	string newin = input;
	vector<string> vals;
	int pos = newin.find(delim);
	while (pos > -1)
	{
		string val = newin.substr(0, pos);
		if (val != "")
			vals.push_back(val);
		newin = newin.substr(pos + 1);
		pos = newin.find(delim);
	}

	vals.push_back(newin);
	return vals;
}

wstring StringManip::utf8ToUtf16(const std::string& utf8Str)
{
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> conv;
	return conv.from_bytes(utf8Str);
}

std::string StringManip::ws2s(const std::wstring& wstr)
{
	using convert_typeX = std::codecvt_utf8<wchar_t>;
	std::wstring_convert<convert_typeX, wchar_t> converterX;

	return converterX.to_bytes(wstr);
}
