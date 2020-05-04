#include "GeometricalDataReport.h"
#include <LogFile.h>
#include <ProteinManager.h>
#include <FoldersFiles.h>

void GeometricalDataReport::printReport(PDBFile* pdb, string fileName1, /*string fileName2,*/ string directory, string geodir)
{//produce data frame report for R reporting
	
	if (directory != "")
	{
		vector<string> pdbs = FoldersFiles::getFilesWithinFolder(directory);
		for (int i = pdbs.size()-1; i >= 0; --i)
		{			
			stringstream status;
			status << i << " out of " << pdbs.size();
			string pdb = pdbs[i];
			try
			{
				string file = directory + pdb;				
				LogFile::getInstance()->writeMessage(status.str() + " - Loading data for PDB1, file=" + pdb);
				PDBFile* pf = ProteinManager::getInstance()->getOrAddPDBFile(pdb, file);
				pf->loadData();
				pf->loadAtoms();
				pf->loadBonds();
				if (pf != nullptr)
				{					
					map<string, ProteinStructure*> versions = pf->getStructureVersions();
					for (map<string, ProteinStructure*>::iterator iter = versions.begin(); iter != versions.end(); ++iter)
					{
						string geofile = geodir + pdb + "_" + iter->first + ".geo.txt";						
						printOneReportWithGeoDef(pf, iter->first, geofile);
					}
				}
			}
			catch (...)
			{
				LogFile::getInstance()->writeMessage("Error loading pdb file " + pdb);
			}
		}		
	}
	else
	{
		map<string, ProteinStructure*> versions = pdb->getStructureVersions();
		for (map<string, ProteinStructure*>::iterator iter = versions.begin(); iter != versions.end(); ++iter)
		{			
			//printOneReport(pdb, iter->first, fileName1 + "_" + iter->first + ".geo.txt", fileName2 + "_" + iter->first + ".geo.txt");
			//printOneReportWithGeoDef(pdb, iter->first, fileName1 + "_" + iter->first + "_geo.txt");
			printOneReportWithGeoDef(pdb, iter->first, fileName1 + "_" + iter->first + "_geo.txt");
		}		
	}
}

void GeometricalDataReport::printOneReportWithGeoDef(PDBFile* pdb, string occupant, string fileName1)
{
	LogFile::getInstance()->writeMessage("Starting Geometric Data report for " + pdb->pdbCode);

	stringstream report;	
	report << "PdbCode,Chain,AminoAcid,AminoNo,PdbAtoms,SecStruct,GeoType,ExperimentalMethod,GeoAtoms,AllAAs,Value\n";
	
	map<string, Chain*> chains = ProteinManager::getInstance()->getChains(pdb->pdbCode, occupant);
	for (map<string, Chain*>::iterator iter = chains.begin(); iter != chains.end(); ++iter)
	{
		Chain* ch = iter->second;
		map<int, AminoAcid*> aminos = ch->getAminoAcids();
		vector<AtomGeo*> v1;
		vector<AtomGeo*> v2;
		vector<AtomGeo*> v3;
		vector<AtomGeo*> v4;
		vector<AtomGeo*> v5;
		vector<AtomGeo*> v6;

		for (map<int, AminoAcid*>::iterator biter = aminos.begin(); biter != aminos.end(); ++biter)
		{
			AminoAcid* aa = biter->second;
			v1 = aa->getAtomBonds(getGeoDefinitions(aa->aminoCode, "BOND"));
			v2 = aa->getAtomCAlphas(getGeoDefinitions(aa->aminoCode, "CALPHA"));
			v3 = aa->getAtomOneFours(getGeoDefinitions(aa->aminoCode, "ONEFOUR"));
			v4 = aa->getAtomAngles(getGeoDefinitions(aa->aminoCode, "ANGLE"));
			v5 = aa->getAtomDihedrals(getGeoDefinitions(aa->aminoCode, "DIHEDRAL"));
			v6 = aa->getAtomImpropers(getGeoDefinitions(aa->aminoCode, "IMPROPER"));

			// we can delete the pointers to AtomGeo as we go
			for (unsigned int i = 0; i < v1.size(); ++i)
			{
				report << getReportString(v1[i], pdb, occupant, "BOND");
				delete v1[i];
				v1[i] = nullptr;
			}

			for (unsigned int i = 0; i < v2.size(); ++i)
			{
				report << getReportString(v2[i], pdb, occupant, "CALPHA");
				delete v2[i];
				v2[i] = nullptr;
			}

			for (unsigned int i = 0; i < v3.size(); ++i)
			{
				report << getReportString(v3[i], pdb, occupant, "ONEFOUR");
				delete v3[i];
				v3[i] = nullptr;
			}

			for (unsigned int i = 0; i < v4.size(); ++i)
			{
				report << getReportString(v4[i], pdb, occupant, "ANGLE");
				delete v4[i];
				v4[i] = nullptr;
			}

			for (unsigned int i = 0; i < v5.size(); ++i)
			{
				report << getReportString(v5[i], pdb, occupant, "DIHEDRAL");
				delete v5[i];
				v5[i] = nullptr;
			}

			for (unsigned int i = 0; i < v6.size(); ++i)
			{
				report << getReportString(v6[i], pdb, occupant, "IMPROPER");
				delete v6[i];
				v6[i] = nullptr;
			}
		}

	}

	ofstream outfile(fileName1);
	if (outfile.is_open())
	{
		outfile << report.str();
	}

	
}

string GeometricalDataReport::getReportString(AtomGeo* ab, PDBFile* pdb, string occupant, string geoType)
{
	stringstream report;
	report << pdb->pdbCode << "_" << occupant << ",";
	report << ab->getChain() << ",";
	report << ab->getAA() << ",";
	report << ab->getId() << ",";
	report << ab->getAtomNos() << ",";
	report << ab->getSS() << ",";
	report << geoType << ",";
	report << pdb->experimentalMethod << ",";
	//report << ab.getAtoms() << ",";
	report << ab->geoDef << ",";
	report << ab->allAAs << ",";
	report << ab->getValue() << "\n";
	return report.str();
}

void GeometricalDataReport::printOneReport(PDBFile* pdb, string occupant, string fileName1)//, string fileName2)
{//produce data frame report for R reporting

	
	LogFile::getInstance()->writeMessage("Starting Scoring Data report for " + pdb->pdbCode);

	stringstream report;	
	//report << "Chain,AminoAcid,Id,SecStruct,DataType,ExperimentalMethod,Atoms,Value\n";
	report << "PdbCode,Chain,AminoAcid,AminoNo,PdbAtoms,SecStruct,GeoType,ExperimentalMethod,GeoAtoms,Value\n";	
	vector<AtomBond> bonds = ProteinManager::getInstance()->getAtomBonds(pdb->pdbCode, occupant);
	vector<AtomAngle> angles = ProteinManager::getInstance()->getAtomAngles(pdb->pdbCode, occupant);
	vector<AtomTorsion> torsions = ProteinManager::getInstance()->getAtomTorsions(pdb->pdbCode, occupant);

	for (unsigned int i = 0; i < bonds.size(); ++i)
	{
		report << pdb->pdbCode << "_" << occupant << ",";
		report << bonds[i].getChain() << ",";
		report << bonds[i].getAA() << ",";
		report << bonds[i].getId() << ",";
		report << bonds[i].getAtomNos() << ",";
		report << bonds[i].getSS() << ",";
		report << "BOND,";
		report << pdb->experimentalMethod << ",";
		report << bonds[i].getAtoms() << ",";
		report << bonds[i].getValue() << "\n";

	}
	for (unsigned int i = 0; i < angles.size(); ++i)
	{
		report << pdb->pdbCode << "_" << occupant << ",";
		report << angles[i].getChain() << ",";
		report << angles[i].getAA() << ",";
		report << angles[i].getId() << ",";
		report << angles[i].getAtomNos() << ",";
		report << angles[i].getSS() << ",";
		report << "ANGLE,";
		report << pdb->experimentalMethod << ",";
		report << angles[i].getAtoms() << ",";
		report << angles[i].getValue() << "\n";

	}
	for (unsigned int i = 0; i < torsions.size(); ++i)
	{
		report << pdb->pdbCode << "_" << occupant << ",";
		report << torsions[i].getChain() << ",";
		report << torsions[i].getAA() << ",";		
		report << torsions[i].getId() << ",";
		report << torsions[i].getAtomNos() << ",";
		report << torsions[i].getSS() << ",";
		report << "TORSION,";
		report << pdb->experimentalMethod << ",";
		report << torsions[i].getAtoms() << ",";
		report << torsions[i].getValue() << "\n";
	}

	//Dual output to the temporary results directory, and also to a database directory for accumulation of data
	ofstream outfile(fileName1);
	if (outfile.is_open())
	{
		outfile << report.str();
	}

	/*if (fileName2 != "")
	{
		ofstream outfile2(fileName2);
		if (outfile2.is_open())
		{
			outfile2 << report.str();
		}
	}*/
}
vector<string> GeometricalDataReport::getGeoDefinitions(string aminoCode, string geoType)
{
	/*
	the file is of the format
	-------------------------------
	AminoAcid, GeoType, Atoms
	#Previousand next atoms, ,
	*,BOND,CP-N
	ILE,ONEFOUR,CD1-CG2
	-------------------------------
	Where a # means ignore, and a * means all amino acids
	*/
	vector<string> vatoms;
	for (unsigned int i = 1; i < geoDefinitions.size(); ++i)//skip the header
	{
		vector<string> row = geoDefinitions[i];		
		if (row.size() == 3)
		{
			string amino = row[0];
			string geo = row[1];
			string atoms = row[2];
			if ((amino == aminoCode || amino == "*") && geoType == geo)
				vatoms.push_back(atoms);
		}
	}		
	return vatoms;
}