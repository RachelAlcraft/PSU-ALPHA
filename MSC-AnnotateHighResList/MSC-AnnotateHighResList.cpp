// MSC-AnnotateHighResList.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <CSVFile.h>
#include <DataFrame.h>
#include <StringManip.h>
#include <sstream>
#include <PDBFile.h>
//#include <LogFile.h>
#include <ProteinManager.h>
#include <LogFile.h>
using namespace std;

/*
This annotates the list of pdbs that are unique and high res with info from their downloaded PDB files
From Mark: their resolution, R and Rfree values and whether or not they have experimental data =structure factors
Also: How many chains, whether there are small molecules
*/

int main()
{
	/*
	This runs in 2 sections because of memroty problems which I will fix another time!
	*/

	string set_label = "2019";
	string runID = "B";

	string filePath = "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-RemoveSimilarity\\";
	
	string inputname = "2019nonsim90.csv";
	string outputname = "2019nonsim90_annotated";

	string dir = "F:\\PSUA\\Code\\PSU-ALPHA\\Project\\PDBData\\";
	string pdbdir = "F:\\PSUA\\ProjectData\\HighResFiles\\";
	
	string annPDBFile = outputname + runID + ".csv";

	CSVFile inPDBs(dir + inputname, ",", true);
	DataFrame annPDBs(dir + annPDBFile);


	
	
	
	
	
	bool success = LogFile::getInstance()->setLogFile("F:\\PSUA\\Code\\PSU-ALPHA\\MSC-AnnotateHighResList\\logger.txt", "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-AnnotateHighResList\\");	
	string aadata = "F:\\PSUA\\Code\\PSU-ALPHA\\Config\\data_aminoinfo.csv";
	ProteinManager::getInstance()->createAminoAcidData(aadata);
	

	annPDBs.headerVector.push_back("PDB"); //pdb code
	annPDBs.headerVector.push_back("RES"); //the resolution
	annPDBs.headerVector.push_back("CLASS"); //header class
	annPDBs.headerVector.push_back("COMPLEX"); //is the structure a complex?
	annPDBs.headerVector.push_back("RVALUE"); // the r value
	annPDBs.headerVector.push_back("RFREE");  // the r free value
	annPDBs.headerVector.push_back("OCCUPANCY"); //any atoms with occupancy less than 1?
	annPDBs.headerVector.push_back("BFACTOR"); //is there ever a b factor > 30? Y or N
	annPDBs.headerVector.push_back("HYDROGENS"); //level of detail of resolution such that hydrogen atoms are in the pdb structure
	annPDBs.headerVector.push_back("STRUCFAC"); //are there structure factors in the pdb
	annPDBs.headerVector.push_back("CHAINS"); //How many chains
	annPDBs.headerVector.push_back("RESIDUES"); //how many residies
	annPDBs.headerVector.push_back("NUCLEOTIDES"); //how many nucleotides
	annPDBs.headerVector.push_back("DATE"); //what date was it deposited
	annPDBs.headerVector.push_back("STATUS"); //NA here comments can be annoted later

	for (unsigned int i = 1; i < inPDBs.fileVector.size(); ++i)
	{
		vector<string> observation;
		string pdb = inPDBs.fileVector[i][0];

		/* shall we run this pdb?*/
		bool shallRun = false;
		string pdbNum = pdb.substr(0, 1);
		if (runID == "A")
		{			
			if (pdbNum == "1" || pdbNum == "2" || pdbNum == "3")
				shallRun = true;
		}
		if (runID == "B")
		{
			if (pdbNum == "4" || pdbNum == "5" || pdbNum == "6" || pdbNum == "7")
				shallRun = true;
		}

		if (shallRun)
		{
			string res = inPDBs.fileVector[i][1];
			// header problems mneans I cannot actually use my pdb class or the whole project stops building, needs fixing TODO
			CSVFile pdbfile(pdbdir + pdb + ".pdb", "@", true); //dummy seperator as I want the whole line
			CSVFile structurefile(pdbdir + pdb + "-sf.cif", "@", false); //dummy seperator only checking file exists
			string sf = structurefile.exists ? "Y" : "N";
			string rval = "NA";
			string rfree = "NA";
			string occ = "NA";
			string bfactor = "NA";
			string hyd = "NA";
			string chains = "NA";
			string name = "NA";
			string date = "NA";
			string complex = "NA";
			string residues = "NA";
			string nucleotides = "NA";
			
			/*
			HEADER    OXIDOREDUCTASE                          17-SEP-98   1BVR
			COMPND   3 CHAIN: A, B, C, D, E, F;
			REMARK   3   R VALUE            (WORKING SET) : 0.161
			REMARK   3   FREE R VALUE                     : 0.178
			*/
			stringstream ss;
			ss << "Annotating " << pdb << " " << i;
			LogFile::getInstance()->writeMessage(ss.str());
			//if (true)
			//{			


			//}
			if (pdbfile.exists)
			{
				//slowly put the functionality into the pdbfile class
				PDBFile* pdbf = ProteinManager::getInstance()->getOrAddPDBFile(pdb, pdbdir + pdb + ".pdb");
				pdbf->loadData();
				pdbf->loadAtoms();
				int iresidues = pdbf->residues;				
				stringstream ssres;
				ssres << iresidues;
				residues = ssres.str();

				int inucleotides = pdbf->nucleotides;
				stringstream ssnuc;
				ssnuc << inucleotides;
				nucleotides = ssnuc.str();


				//occupancy
				bool occupancy = ProteinManager::getInstance()->hasOccupancy(pdb,"A");
				occupancy ? occ = "Y" : occ = "N";

				//BFactor
				double bf = ProteinManager::getInstance()->maxBFactor(pdb,"A");
				stringstream ssbf;
				ssbf.precision(5);
				ssbf << bf;
				bfactor = ssbf.str();

				//hydrogens - level of detail of resolution such that hydrogen atoms are in the pdb structure
				bool hydrogens = ProteinManager::getInstance()->hasHydrogens(pdb,"A");
				hydrogens ? hyd = "Y" : hyd = "N";

				//tbhis needs to be moved into the pdbfile class, but memory and state need to be sorted out
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
			observation.push_back(occ);
			observation.push_back(bfactor);
			observation.push_back(hyd);
			observation.push_back(sf);
			observation.push_back(chains);
			observation.push_back(residues);
			observation.push_back(nucleotides);
			observation.push_back(date);
			observation.push_back(set_label);
			annPDBs.fileVector.push_back(observation);
			ProteinManager::getInstance()->deletePdbs();//keep memory clear
		}
	}
	LogFile::getInstance()->writeMessage("Success, now printing");
	annPDBs.print();
}


