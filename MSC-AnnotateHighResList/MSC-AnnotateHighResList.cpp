// MSC-AnnotateHighResList.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <CSVFile.h>
#include <DataFrame.h>
#include <StringManip.h>
//#include <PDBFile.h>
//#include <LogFile.h>
//#include <ProteinManager.h>
using namespace std;

/*
This annotates the list of pdbs that are unique and high res with info from their downloaded PDB files
From Mark: their resolution, R and Rfree values and whether or not they have experimental data =structure factors
Also: How many chains, whether there are small molecules
*/

int main()
{

	string dir = "F:\\PSUA\\Code\\PSU-ALPHA\\Project\\PDBData\\";
	string pdbdir = "F:\\PSUA\\ProjectData\\HighResFiles\\";
	string highResFile = "highres_and_unique.txt";
	string annHighResFile = "ann_highres_unique.txt";

	CSVFile highRes(dir + highResFile, ",",true);
	DataFrame annHighRes(dir + annHighResFile);
	annHighRes.headerVector.push_back("PDB");
	annHighRes.headerVector.push_back("RES");
	annHighRes.headerVector.push_back("CLASS");
	annHighRes.headerVector.push_back("COMPLEX");
	annHighRes.headerVector.push_back("RVALUE");
	annHighRes.headerVector.push_back("RFREE");
	annHighRes.headerVector.push_back("STRUCFAC");
	annHighRes.headerVector.push_back("CHAINS");
	annHighRes.headerVector.push_back("DATE");
	annHighRes.headerVector.push_back("COMMENTS");

	for (unsigned int i = 1; i < highRes.fileVector.size(); ++i)
	{
		vector<string> observation;
		string pdb = highRes.fileVector[i][0];
		string res = highRes.fileVector[i][1];
		// header problems mneans I cannot actually use my pdb class or the whole project stops building, needs fixing TODO
		CSVFile pdbfile(pdbdir + pdb + ".pdb","@",true); //dummy seperator as I want the whole line
		CSVFile structurefile(pdbdir + pdb + "-sf.cif", "@",false); //dummy seperator only checking file exists
		string sf = structurefile.exists ? "Y" : "N";
		string rval = "NA";
		string rfree = "NA";
		string chains = "NA";
		string name = "NA";
		string date = "NA";
		string complex = "NA";
		/*
		HEADER    OXIDOREDUCTASE                          17-SEP-98   1BVR
		COMPND   3 CHAIN: A, B, C, D, E, F;
		REMARK   3   R VALUE            (WORKING SET) : 0.161
		REMARK   3   FREE R VALUE                     : 0.178
		*/
		cout << "Annotating " << pdb << " " << i << "\n";
		if (pdbfile.exists)
		{
			bool foundR = false;
			bool foundFree = false;
			bool foundChain = false;
			bool foundName = false;
			bool foundTitle = false;
			unsigned int count = 0;
			while (!(foundName && foundR && foundFree && foundChain && foundTitle) && count < pdbfile.fileVector.size())
			{
				if (pdbfile.fileVector[count].size() > 0)
				{
					string line = pdbfile.fileVector[count][0];
					if (!foundR)
					{						
						int pos1 = line.find("R VALUE");
						int pos2 = line.find("WORKING + TEST SET");
						if (pos1 > -1 && pos2 > -1)
						{
							vector<string> rvalvec = StringManip::stringToVector(line, ":");
							rval = StringManip::trim(rvalvec[1]);
							foundR = true;
						}
					}
					
					if (!foundR)
					{					
						int pos1 = line.find("R VALUE");
						int pos2 = line.find("WORKING SET");
						if (pos1 > -1 && pos2 > -1)
						{
							vector<string> rvalvec = StringManip::stringToVector(line, ":");
							rval = StringManip::trim(rvalvec[1]);
							foundR = true;
						}
					}

					if (!foundFree && foundR)//must come after the R Value
					{
						int pos3 = line.find("FREE R VALUE");
						if (pos3 > -1)
						{
							vector<string> rvalvec = StringManip::stringToVector(line, ":");
							rfree = StringManip::trim(rvalvec[1]);
							foundFree = true;
						}
					}

					if (!foundChain)
					{
						int pos4 = line.find("3 CHAIN");
						if (pos4 > -1)
						{
							vector<string> chainvec = StringManip::stringToVector(line, ":");
							chains = StringManip::removeChar(chainvec[1], ",");
							chains = StringManip::removeChar(chains, ";");
							chains = StringManip::removeChar(chains, " ");
							foundChain = true;
						}
					}

					if (!foundName)
					{
						int pos5 = line.find("HEADER");
						if (pos5 > -1)
						{
							vector<string> namevec = StringManip::stringToVector(line, "  ");
							name = StringManip::removeChar(namevec[1], " ");							
							name = StringManip::removeChar(name, ",");
							date = StringManip::removeChar(namevec[2], " ");
							foundName = true;
						}
					}

					if (!foundTitle)
					{
						int pos5 = line.find("TITLE");
						if (pos5 > -1)
						{
							vector<string> namevec = StringManip::stringToVector(line, "  ");
							string title = namevec[1];
							int pos6 = title.find("COMPLEX");
							complex = pos6 > -1 ? "Y" : "N";							
							foundTitle = true;
						}
					}
				}
				++count;
			}
		}
		observation.push_back(pdb);
		observation.push_back(res);
		observation.push_back(name);
		observation.push_back(complex);
		observation.push_back(rval);
		observation.push_back(rfree);
		observation.push_back(sf);
		observation.push_back(chains);
		observation.push_back(date);
		observation.push_back("NA");
		annHighRes.fileVector.push_back(observation);
	}
	annHighRes.print();
}


