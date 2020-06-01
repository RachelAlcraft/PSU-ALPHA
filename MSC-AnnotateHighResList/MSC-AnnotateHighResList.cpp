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

	//string set_label = "HIGHRES";
	string runID = "X";

	//string filePath = "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-RemoveSimilarity\\";
	
	// USER UBNPOUT FILE PATHS AND NAMES
	//string inputname = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_HighRes_90.csv";
	//string outputname = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_HighRes_90_annotated";

	string inputname = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_2019_95.csv";
	string outputname = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_2019_95_annotated";

	//string inputname = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_2018_95.csv";
	//string outputname = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_2018_95_annotated";
		
	string pdbdir = "F:\\Code\\BbkTransfer\\pdbfiles\\pdbdata\\";

	string aadata = "F:\\PSUA\\Code\\PSU-ALPHA\\Config\\data_aminoinfo.csv";

	// CODE BEGINS
	
	string annPDBFile = outputname + runID + ".csv";

	CSVFile inPDBs(inputname, ",", true);
	DataFrame annPDBs(annPDBFile);
	
	bool success = LogFile::getInstance()->setLogFile("logger.txt", "");
	
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
	annPDBs.headerVector.push_back("INSTITUTION"); //who deposited it
	annPDBs.headerVector.push_back("SOFTWARE"); //how was it refined?
	annPDBs.headerVector.push_back("SEQUENCE"); //how was it refined?
	annPDBs.headerVector.push_back("EXPMETHOD"); //always xray
	//annPDBs.headerVector.push_back("STATUS"); //NA here comments can be annoted later

	for (unsigned int i = 1; i < inPDBs.fileVector.size(); ++i)
	{
		vector<string> observation;
		string pdb = inPDBs.fileVector[i][0];

		//if (pdb == "6JX3")
		{

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
			if (runID == "X")
				shallRun = true;

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
				string complex = "NA";
				string residues = "NA";
				string nucleotides = "NA";
				string date = "NA";
				string institution = "NA";
				string software = "NA";
				string sequence = "NA";

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
			

					int inucleotides = pdbf->nucleotides;
					stringstream ssnuc;
					ssnuc << inucleotides;
					nucleotides = ssnuc.str();


					//occupancy
					bool occupancy = ProteinManager::getInstance()->hasOccupancy(pdb, "A");
					occupancy ? occ = "Y" : occ = "N";

					//BFactor
					double bf = ProteinManager::getInstance()->maxBFactor(pdb, "A");
					stringstream ssbf;
					ssbf.precision(5);
					ssbf << bf;
					bfactor = ssbf.str();

					//hydrogens - level of detail of resolution such that hydrogen atoms are in the pdb structure
					bool hydrogens = ProteinManager::getInstance()->hasHydrogens(pdb, "A");
					hydrogens ? hyd = "Y" : hyd = "N";

					institution = pdbf->institution;
					software = pdbf->software;

					bool nullmodel = pdbf->nullModel;
										
					vector<string> seqs = pdbf->getSequence();
					sequence = "";
					bool differ = false;
					unsigned int iresidues = 0;
					for (unsigned int r = 0; r < seqs.size(); ++r)
					{
						if (seqs[r].length() > iresidues)
							iresidues = seqs[r].length();
					}
					
					for (unsigned int r = 0; r < seqs.size(); ++r)
					{
						if (r == 0)
						{
							sequence = seqs[r];
						}
						else
						{
							if (seqs[r] != sequence)
								differ = true;
						}												
					}

					sequence = "";
					if (nullmodel || differ) //then they are independently verified chains or different
					{
						for (unsigned int r = 0; r < seqs.size(); ++r)						
							sequence += seqs[r];						
					}
					else // all the chains are the same and we only want 1
					{
						if (seqs.size() > 0)
							sequence = seqs[0];
					}

					//log any podb files that need chains removed
					if (seqs.size() > 1)
					{
						if (!differ && !nullmodel)
						{
							LogFile::getInstance()->writeMessage("This structure has redundant chains");
						}
					}



					stringstream ssres;
					ssres << iresidues; // this is the max chain not the total residues
					residues = ssres.str();

 					rval = pdbf->rvalue;
					rfree = pdbf->rfree;
					chains = pdbf->maxChain();
					name = pdbf->proteinclass;
					complex = pdbf->inComplex;
					date = pdbf->date;
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
				observation.push_back(institution);
				observation.push_back(software);
				observation.push_back(sequence);
				observation.push_back("XR");
				//observation.push_back(set_label);
				annPDBs.fileVector.push_back(observation);

				ProteinManager::getInstance()->deletePdbs();//keep memory clear
			}
		}
	}
	LogFile::getInstance()->writeMessage("Success, now printing");
	annPDBs.print();
}


