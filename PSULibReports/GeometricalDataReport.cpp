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
	report << "PdbCode,Occupant,Chain,AminoNo,GeoAtoms,AminoCode,AminoNos,AminoCodes,AtomNos,SecStruct,GeoType,Value\n";
	
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
		vector<AtomGeo*> v7;

		for (map<int, AminoAcid*>::iterator biter = aminos.begin(); biter != aminos.end(); ++biter)
		{
			AminoAcid* aa = biter->second;
			v1 = aa->getAtomDistance(getGeoDefinitions(aa->aminoCode, "BOND"),"BOND");
			v2 = aa->getAtomDistance(getGeoDefinitions(aa->aminoCode, "CALPHA"),"CALPHA");
			v3 = aa->getAtomDistance(getGeoDefinitions(aa->aminoCode, "ONEFOUR"),"ONEFOUR");
			v4 = aa->getAtomAngles(getGeoDefinitions(aa->aminoCode, "ANGLE"),"ANGLE");
			v5 = aa->getAtomDihedrals(getGeoDefinitions(aa->aminoCode, "DIHEDRAL"),"DIHEDRAL");
			v6 = aa->getAtomDihedrals(getGeoDefinitions(aa->aminoCode, "IMPROPER"),"IMPROPER");
			v7 = aa->getAtomDistance(getGeoDefinitions(aa->aminoCode, "INTER"),"INTER");

			// we can delete the pointers to AtomGeo as we go
			for (unsigned int i = 0; i < v1.size(); ++i)
			{
				report << getReportString(v1[i], pdb, occupant);
				delete v1[i];
				v1[i] = nullptr;
			}

			for (unsigned int i = 0; i < v2.size(); ++i)
			{
				report << getReportString(v2[i], pdb, occupant);
				delete v2[i];
				v2[i] = nullptr;
			}

			for (unsigned int i = 0; i < v3.size(); ++i)
			{
				report << getReportString(v3[i], pdb, occupant);
				delete v3[i];
				v3[i] = nullptr;
			}

			for (unsigned int i = 0; i < v4.size(); ++i)
			{
				report << getReportString(v4[i], pdb, occupant);
				delete v4[i];
				v4[i] = nullptr;
			}

			for (unsigned int i = 0; i < v5.size(); ++i)
			{
				report << getReportString(v5[i], pdb, occupant);
				delete v5[i];
				v5[i] = nullptr;
			}

			for (unsigned int i = 0; i < v6.size(); ++i)
			{
				report << getReportString(v6[i], pdb, occupant);
				delete v6[i];
				v6[i] = nullptr;
			}
			
			for (unsigned int i = 0; i < v7.size(); ++i)
			{
				report << getReportString(v7[i], pdb, occupant);
				delete v7[i];
				v7[i] = nullptr;
			}
		}

	}

	ofstream outfile(fileName1);
	if (outfile.is_open())
	{
		outfile << report.str();
	}

	
}

string GeometricalDataReport::getReportString(AtomGeo* ab, PDBFile* pdb, string occupant)
{
	//PdbCode, Occupant, Chain, AminoNo, GeoAtoms, AminoCode, AminoNos, AminoCodes, AtomNos, SecStruct, GeoType, Value
	stringstream report;
	report << pdb->pdbCode << ",";
	report << occupant << ",";
	report << ab->getChain() << ",";
	report << ab->getAminoId() << ",";
	report << ab->getGeoDef() << ",";
	report << ab->getAA() << ",";
	report << ab->getAminoNos() << ",";
	report << ab->getAACodes() << ",";
	report << ab->getAtomNos() << ",";
	report << ab->getSS() << ",";
	report << ab->getAlias() << ",";
	report << quickRound(ab->getValue()) << "\n"; // 3 decimal places
	return report.str();
}

double GeometricalDataReport::quickRound(double val)
{	 
	double dVal = val * 1000;
	int iVal = (int)round(dVal);
	dVal = (double)iVal / 1000;
	return dVal;
}

void GeometricalDataReport::printOneReport(PDBFile* pdb, string occupant, string fileName1)//, string fileName2)
{//produce data frame report for R reporting

	
	LogFile::getInstance()->writeMessage("Starting Scoring Data report for " + pdb->pdbCode);

	stringstream report;	
	//report << "Chain,AminoAcid,Id,SecStruct,DataType,ExperimentalMethod,Atoms,Value\n";
	report << "PdbCode,Chain,AminoAcid,AminoNo,AtomNo,AtomNos,SecStruct,GeoType,ExperimentalMethod,GeoAtoms,Value\n";	
	vector<AtomBond> bonds = ProteinManager::getInstance()->getAtomBonds(pdb->pdbCode, occupant);
	vector<AtomAngle> angles = ProteinManager::getInstance()->getAtomAngles(pdb->pdbCode, occupant);
	vector<AtomTorsion> torsions = ProteinManager::getInstance()->getAtomTorsions(pdb->pdbCode, occupant);

	for (unsigned int i = 0; i < bonds.size(); ++i)
	{
		report << pdb->pdbCode << "_" << occupant << ",";
		report << bonds[i].getChain() << ",";
		report << bonds[i].getAA() << ",";
		report << bonds[i].getAminoId() << ",";		
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
		report << angles[i].getAminoId() << ",";	
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
		report << torsions[i].getAminoId() << ",";		
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
vector<pair<string,string>> GeometricalDataReport::getGeoDefinitions(string aminoCode, string geoType)
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
	vector<pair<string,string>> vatoms;
	for (unsigned int i = 1; i < geoDefinitions.size(); ++i)//skip the header
	{
		vector<string> row = geoDefinitions[i];		
		if (row.size() == 4)
		{
			string amino = row[0];
			string geo = row[1];
			string atoms = row[2];
			string alias = row[3];
			if ((amino == aminoCode || amino == "*") && geoType == geo)
				vatoms.push_back(pair<string,string>(atoms,alias));			
		}
	}		
	return vatoms;
}