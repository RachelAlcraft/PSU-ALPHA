// PHD-RemoveSimilarity.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <CSVFile.h>
//#include <algorithm>
#include <map>
#include <set>
#include <DataFrame.h>
#include <sstream>
#include <iomanip>


using namespace std;


	//////////////////////////////////////////////////////////////////////////
			///////// HELPER FUNCTIONS ///////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

double getResolution(string pdbCode, map<string, double>* resDictionary)
{
	map<string, double>::iterator iter = resDictionary->find(pdbCode); 
	if (iter != resDictionary->end())
		return iter->second;
	else
		return 0;
}


int getHomologyPosition(string pdbCode, vector<set<string>>* similarLists)
{	
	for (unsigned int i = 0; i < similarLists->size(); ++i)
	{
		set<string> thisHomology = (*similarLists)[i];
		set<string>::iterator iter = thisHomology.find(pdbCode);
		if (iter != thisHomology.end())
			return i;
	}
	return -1;
}
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

int main()
{
	// ***  USER INPUT *********************************************************
	string inputPath = "F:\\PSUA\\Code\\PSU-ALPHA\\PHD-RemoveSimilarity\\";

	string inputname = "InOut\\rcsb_pdb_ids_high.txt";
	string outputFilename = "InOut\\high_nonsim90_1_3.csv";
	//string rejectedFileName = "InOut\\high_sim.csv";
	double resolutionlimit = 1.3; // 0 for no limit
	string similarity = "90";
	// *************************************************************************

	//////////////////////////////////////////////////////////////////////////
	cout << "1. Load the files\n";
	CSVFile* similarStructures = new CSVFile(inputPath + "Data\\bc-100.out", " ", true);
	if (similarity == "90")
		similarStructures = new CSVFile(inputPath + "Data\\bc-90.out", " ", true);
	else if (similarity == "95")
		similarStructures = new CSVFile(inputPath + "Data\\bc-95.out", " ", true);			
	CSVFile pdblist(inputPath + inputname, ",", true);	
	CSVFile cmpd_res(inputPath + "Data\\cmpd_res.idx", ";", true);

	//////////////////////////////////////////////////////////////////////////
	cout << "2. Create a dictionary of pdb to resolution\n";		
	map<string, double> pdb_res;
	for (unsigned int i = 4; i < cmpd_res.fileVector.size(); ++i)
	{
		std::cout << "2.  " << i << "\\" << cmpd_res.fileVector.size() << "\n";
		vector <string> entry = cmpd_res.fileVector[i];
		string pdb = entry[0];
		pdb = pdb.replace(pdb.find("\t"), 1, "");
		string strres = entry[1];

		double res = atof(strres.c_str());
		if (res < 0.0000001)//arbitrary small number
		{
			//std::cout << "Not x-ray " << pdb << "\n";
		}
		else
		{
			if (pdb_res.find(pdb) == pdb_res.end())
			{
				pdb_res.insert(pair<string, double>(pdb, res));
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	cout << "3. Turn the similarity list into a list of sets\n";	
	vector<set<string>> similarLists;	
	int totalI = similarStructures->fileVector.size();
	for (unsigned int i = 0; i < similarStructures->fileVector.size(); ++i)
	{
		std::cout << "3.  i=" << i <<  "\\" << totalI << "\n";
		set<string> row;

		vector<string> this_row = similarStructures->fileVector[i];
		for (unsigned int j = 0; j < this_row.size(); ++j)
		{
			string pdb = this_row[j]; // includes _chain letter
			pdb = pdb.substr(0, 4);
			row.insert(pdb);
		}
		similarLists.push_back(row);
	}

	//////////////////////////////////////////////////////////////////////////
	cout << "4. Go through the inputs and add to either the yes list of the homology dictionary\n";
	vector <string> yeslist;
	map<int, string> homologyList;
	for (unsigned int i = 0; i < pdblist.fileVector[0].size(); ++i)
	{		
		string pdb = pdblist.fileVector[0][i];
		std::cout << "4.  " << i << "\\" << pdblist.fileVector[0].size() << " " << pdb << "\n";
		int simPos = getHomologyPosition(pdb, &similarLists);
		if (simPos == -1)
		{
			yeslist.push_back(pdb);
		}
		else
		{
			if (homologyList.find(simPos) == homologyList.end())
			{
				homologyList[simPos] = pdb;
			}
			else
			{
				string currentPdb = homologyList[simPos];
				double currentResolution = getResolution(currentPdb, &pdb_res);
				double candidateResolution = getResolution(pdb, &pdb_res);

				if (candidateResolution < currentResolution)
					homologyList[simPos] = pdb;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	cout << "5. Add the homology and yes together\n";
	vector<string> pdbsAllResolution;
	for (unsigned int i = 0; i < yeslist.size(); ++i)
		pdbsAllResolution.push_back(yeslist[i]);

	for (map<int, string>::iterator iter = homologyList.begin(); iter != homologyList.end(); ++iter)
		pdbsAllResolution.push_back(iter->second);

	

	//////////////////////////////////////////////////////////////////////////
	cout << "6. Only keep those of the desired resolution\n";
	vector<string> pdbsFinalList;
	for (unsigned int i = 0; i < pdbsAllResolution.size(); ++i)
	{
		std::cout << "6.  "  << i << "\\" << pdbsAllResolution.size() << "\n";
		string pdb = pdbsAllResolution[i];
		double resolution = getResolution(pdb, &pdb_res);
		if (resolution > 0 && (resolutionlimit == 0 || resolution <= resolutionlimit))
			pdbsFinalList.push_back(pdb);
	}


	//////////////////////////////////////////////////////////////////////////
	cout << "7. Print the non homologous list out\n";
	DataFrame data_nonsim(inputPath + outputFilename);
	data_nonsim.headerVector.push_back("PDB");
	data_nonsim.headerVector.push_back("RES");
	int pdbcount = pdbsFinalList.size();
	int count = 0;
	for (unsigned int i = 0; i < pdbsFinalList.size(); ++i)
	{		
		cout << "7.  " << i << "/" << pdbcount << "\n";
		string pdb = pdbsFinalList[i];
		double resolution = getResolution(pdb, &pdb_res);				
		//then add it to our list
		vector<string> observation;
		observation.push_back(pdb);
		stringstream ss;
		ss << setprecision(4);
		ss << resolution;
		observation.push_back(ss.str());
		data_nonsim.fileVector.push_back(observation);				
	}

	//finally print the list of pdb files that are high res and not 100% sequence similar
	cout << "8. Now print to: " << inputPath + outputFilename << "\n";
	data_nonsim.print();	
}




