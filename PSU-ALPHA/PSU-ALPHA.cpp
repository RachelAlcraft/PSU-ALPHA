// PSU-ALPHA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "LogFile.h"
#include "InputParams.h"
#include "ProteinManager.h"
#include <RamaReport.h>

using namespace std;

int main()
{
	// Run is controlled by an input file that describes what the run should be
	//The params file is forced to be next to the exe as it is the only place we can guarantee
	InputParams::getInstance()->setInputFile("InputParams.txt");
	string OUTPATH = InputParams::getInstance()->getParam("OUTPUTPATH");
	string INPATH = InputParams::getInstance()->getParam("CONFIGPATH");

	// Set up log file for output
	string runId = LogFile::getInstance()->runId();
	OUTPATH = OUTPATH + runId + "\\";	
	bool success = LogFile::getInstance()->setLogFile(OUTPATH + "logger.pdb",OUTPATH);
	if (success)
	{
		LogFile::getInstance()->writeMessage("");
		LogFile::getInstance()->writeMessage("**********************************************************");
		LogFile::getInstance()->writeMessage("Starting Protein Struture Utility (PSU Rachel Alcraft 2020)");




		//The run parameters control the flow of the prgram
		string PDB1 = InputParams::getInstance()->getParam("PDB1");
		string PDB2 = InputParams::getInstance()->getParam("PDB2");
		string RAMA = InputParams::getInstance()->getParam("RAMA");
		string CALPHA = InputParams::getInstance()->getParam("CALPHA");
		string RMSDFIX = InputParams::getInstance()->getParam("RMSDFIX");
		string RMSDOPT = InputParams::getInstance()->getParam("RMSDOPT");
		string ALIGNMENT = InputParams::getInstance()->getParam("ALIGNMENT");


		//// INPUT CONSTRUCTION ////
		string inConfig = INPATH;
		string inPDB = INPATH + "PDB\\" + PDB1 + ".pdb";
		string ramareport = OUTPATH + "Reports\\" + PDB1 + "_torsion.txt";
		string calphareport = OUTPATH + "Reports\\" + PDB1 + "_calpha.txt";

		// Set up data manager config ////////////////////
		ProteinManager::getInstance()->createConfigData(inConfig);

		//Load pdb data if we have it first
		PDBFile* pdb1 = nullptr;
		if (PDB1 != "")
		{
			LogFile::getInstance()->writeMessage("Loading data fro PDB1, file=" + inPDB);
			pdb1 = ProteinManager::getInstance()->getOrAddPDBFile(PDB1, inPDB);
			pdb1->loadData();
		}

		//Shall we run a ramachandran report?
		if (PDB1 != "" && RAMA == "TRUE")
		{
			LogFile::getInstance()->writeMessage("RAMACHANDRAN PLOT REPORT, outfile=" + ramareport);
			
			RamaReport rr;
			rr.printReport(pdb1, ramareport);

			

		}
		//Shall we run a CAlpha report
		if (PDB1 != "" && CALPHA == "TRUE")
		{
			LogFile::getInstance()->writeMessage("CALPHA Contact Map REPORT, outfile=" + calphareport);

		}
		//Shall we run an RMSD report between 2 structures with fixed positions?
		if (PDB1 != "" && PDB2 != "" && RMSDFIX == "TRUE")
		{
			if (ALIGNMENT != "")//whether a fixed calpha or a pairing needed
			{

			}
			else
			{

			}
		}
		//Shall we run an RMSD report between 2 structures optimising the alignments?
		if (PDB1 != "" && PDB2 != "" && RMSDOPT == "TRUE")
		{
			if (ALIGNMENT != "")//whether a fixed calpha or a pairing needed
			{

			}
			else
			{

			}
		}
	}
	else
	{
		//It is not a good idea to proceed with a run that could take hours with no log files so don't even try
		cout << "ERROR!!! THE LOG FILE COULD NOT BE CREATED!!! Check paths in config file and write access to directories\n";
	}


}



