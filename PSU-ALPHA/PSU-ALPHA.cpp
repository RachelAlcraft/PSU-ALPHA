// PSU-ALPHA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "LogFile.h"
#include "InputParams.h"
#include "ProteinManager.h"
#include <RamaReport.h>
#include <CAlphaReport.h>
#include <FASTAFile.h>
#include <RMSD.h>
#include <RMSDReport.h>
#include <GeometricalDataReport.h>
#include <GeometricalAggregationReport.h>

using namespace std;

int main()
{
	// Run is controlled by an input file that describes what the run should be
	//The params file is forced to be next to the exe as it is the only place we can guarantee
	InputParams::getInstance()->setInputFile("InputParams.txt");
	string OUTPATH = InputParams::getInstance()->getParam("OUTPUTPATH");
	string INPATH = InputParams::getInstance()->getParam("CONFIGPATH");
	string RUNID = InputParams::getInstance()->getParam("RUNID");

	// Set up log file for output
	if (RUNID == "")
		RUNID = LogFile::getInstance()->runId();
	OUTPATH = OUTPATH + RUNID + "\\";
	bool success = LogFile::getInstance()->setLogFile(OUTPATH + "logger.txt",OUTPATH);
	if (success)
	{
		LogFile::getInstance()->writeMessage("");
		LogFile::getInstance()->writeMessage("**********************************************************");
		LogFile::getInstance()->writeMessage("Starting Protein Struture Utility (PSU Rachel Alcraft 2020)");




		//The run parameters control the flow of the prgram
		string PDB1 = InputParams::getInstance()->getParam("PDB1");
		string PDB2 = InputParams::getInstance()->getParam("PDB2");
		string RAMA = InputParams::getInstance()->getParam("RAMA");
		string CONTACT = InputParams::getInstance()->getParam("CONTACT");
		string CONTACTCHAIN1 = InputParams::getInstance()->getParam("CONTACTCHAIN1");
		string CONTACTCHAIN2 = InputParams::getInstance()->getParam("CONTACTCHAIN2");
		string RMSDFIX = InputParams::getInstance()->getParam("RMSDFIX");
		string RMSDOPT = InputParams::getInstance()->getParam("RMSDOPT");
		string ALIGNMENT = InputParams::getInstance()->getParam("ALIGNMENT");
		string RMSDCONTACT = InputParams::getInstance()->getParam("RMSDCONTACT");
		string GEODATABASE = InputParams::getInstance()->getParam("GEODATABASE");
		string PDBDATABASE = InputParams::getInstance()->getParam("PDBDATABASE");
		string GEOREPORT = InputParams::getInstance()->getParam("GEOREPORT");
		string GEOAGGREGATION = InputParams::getInstance()->getParam("GEOAGGREGATION");


		//// INPUT CONSTRUCTION ////
		string inConfig = INPATH;
		string inPDB = INPATH + "PDB\\" + PDB1 + ".pdb";
		string inPDB2 = INPATH + "PDB\\" + PDB2 + ".pdb";
		string ramareport = OUTPATH + "Reports\\" + PDB1 + "_torsion.txt";
		string calphareport = OUTPATH + "Reports\\" + PDB1 + "_calpha.txt";

		// Set up data manager config ////////////////////
		ProteinManager::getInstance()->createConfigData(inConfig);

		//Load pdb data if we have it first
		PDBFile* pdb1 = nullptr;
		if (PDB1 != "")
		{
			LogFile::getInstance()->writeMessage("Loading data for PDB1, file=" + inPDB);
			pdb1 = ProteinManager::getInstance()->getOrAddPDBFile(PDB1, inPDB);
			pdb1->loadData(false);
		}

		PDBFile* pdb2 = nullptr;
		if (PDB2 != "")
		{
			LogFile::getInstance()->writeMessage("Loading data for PDB2, file=" + inPDB2);
			pdb2 = ProteinManager::getInstance()->getOrAddPDBFile(PDB2, inPDB2);
			pdb2->loadData(false);
		}

		//Shall we run a ramachandran report?
		if (PDB1 != "" && RAMA == "TRUE")
		{
			LogFile::getInstance()->writeMessage("RAMACHANDRAN PLOT REPORT, outfile=" + ramareport);			
			RamaReport rr;
			rr.printReport(pdb1, ramareport);
		}
		//Shall we run a CAlpha report
		if (PDB1 != "" && CONTACT == "TRUE")
		{
			LogFile::getInstance()->writeMessage("CALPHA Contact Map REPORT, outfile=" + calphareport);
			CAlphaReport ca;
			ca.printReport(pdb1, CONTACTCHAIN1, CONTACTCHAIN2, calphareport);

		}
		//Shall we run an RMSD report between 2 structures with fixed positions?
		if (PDB1 != "" && PDB2 != "" && (RMSDFIX == "TRUE" || RMSDOPT == "TRUE"))
		{
			string rmsdreport = OUTPATH + "Reports\\" + PDB1 + "_" + PDB2 + "_rmsd.txt";
			string fileroot = OUTPATH + "Reports\\" + PDB1 + PDB2 + "_";
			LogFile::getInstance()->writeMessage("Running RMSD report to file=" + rmsdreport);

			bool alignment = false;
			bool optimised = false;
			FastaFile* ff = nullptr;
			if (ALIGNMENT != "")//whether a fixed calpha or a pairing needed
			{
				alignment = true;
				//and create the ff
			}
			if (RMSDOPT != "")//whether a fixed calpha or a pairing needed
			{
				LeastSquares* ls = new LeastSquares(pdb1,pdb2,alignment);
				RMSDReport rrmsd;
				rrmsd.printLeastSquaresReport(ls, rmsdreport, fileroot);
				
			}
			else
			{
				RMSD* rmsd = new RMSD(pdb1, pdb2, ff, alignment, optimised);
				RMSDReport rrmsd;
				rrmsd.printReport(rmsd, rmsdreport, optimised, fileroot);
			}
			
			

			if (RMSDCONTACT == "TRUE")
			{
				string rmsdcontact = OUTPATH + "Reports\\" + PDB1 + "_rmsdcontact.txt";
				LogFile::getInstance()->writeMessage("RMSD CALPHA Contact Map REPORT, outfile=" + rmsdcontact);
				CAlphaReport ca;
				ca.printMultiReport(pdb1,pdb2, rmsdcontact,true);
			}

			
		}
		//Shall we write out a database file of the geometric features of this pdb?
		if (GEOREPORT == "TRUE")
		{
			string geodata1 = OUTPATH + "Reports\\" + PDB1 + "_geometricfeatures.csv";
			string geodata2 = GEODATABASE + PDB1 + "_geometricfeatures.csv";
			LogFile::getInstance()->writeMessage("Outputting geometric features to, outfile=" + geodata1 + " and " + geodata2);
			GeometricalDataReport gdr;
			gdr.printReport(pdb1, geodata1, geodata2,PDBDATABASE,GEODATABASE);

		}
		//Shall we aggregate all of the existing geometrical data into a probability distribution file?
		if (GEOAGGREGATION == "TRUE")
		{	
			//string geodata1 = OUTPATH + "Reports\\geoprobdist.csv";
			string geodata = GEODATABASE;
			LogFile::getInstance()->writeMessage("Outputting geometric features to, outfile=" + geodata);
			GeometricalAggregationReport gar;
			gar.printReport(geodata);

		}

		LogFile::getInstance()->writeMessage("########### PSU-Alpha concluded ##############");
	}
	else
	{
		//It is not a good idea to proceed with a run that could take hours with no log files so don't even try
		cout << "ERROR!!! THE LOG FILE COULD NOT BE CREATED!!! Check paths in config file and write access to directories\n";
	}


}



