#include "LogFile.h"
#include <Windows.h>


LogFile* LogFile::instance = 0;

LogFile::LogFile()
{
	_fileName = "";
	time(&_start);
	_starttime = localtime(&_start);
	stringstream run;
	run << (_starttime->tm_year) - 100 << "_" << (_starttime->tm_mon) + 1 << "_";
	run << _starttime->tm_mday << "_" << _starttime->tm_hour << "_" << _starttime->tm_min << "_";
	run << _starttime->tm_sec;
	_runId = run.str();	
}

LogFile* LogFile::getInstance()
{
	if (!instance)
		instance = new LogFile();
	return instance;
}

void LogFile::writeMessage(string msg)
{
	stringstream message;
	time_t now;
	time(&now);
	tm* nowtime = localtime(&now);
	double diff = difftime(now, _start);

	message << nowtime->tm_hour << ":" << nowtime->tm_min << ":" << nowtime->tm_sec << "(elap:" << diff << "):\t" << msg << "\n";
	cout << message.str();

	ofstream myfile(_fileName, std::ios_base::app);
	if (myfile.is_open())
	{
		myfile << message.str();
	}

}

string LogFile::runId()
{
	return _runId;
}

bool LogFile::setLogFile(string fileName, string path)
{	
	_fileName = fileName;
	//Now need to ensure there is a directory here
	return CreateFolder(path.c_str());
}


bool LogFile::CreateFolder(string path)
{
	wstring wpath = utf8ToUtf16(path);	
	LPCWSTR lp = wpath.c_str(); 
	bool success = true;
	if (CreateDirectoryW(lp, NULL))
	{
		string reppath = path + "\\Reports\\";
		wstring wrpath = utf8ToUtf16(reppath);
		LPCWSTR lpr = wrpath.c_str();
		if (!CreateDirectoryW(lpr, NULL))		
			success = false;		
	}
	else	
		success = true; //TODO need to check if the directory exists
	
	return success;
}

wstring LogFile::utf8ToUtf16(const std::string& utf8Str)
{
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> conv;
	return conv.from_bytes(utf8Str);
}


