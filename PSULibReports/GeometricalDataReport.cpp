#include "GeometricalDataReport.h"
#include <LogFile.h>
#include <ProteinManager.h>
#include <FoldersFiles.h>
#include <DataFrame.h>


void GeometricalDataReport::printReport(PDBFile* pdb,string directory, string geodir, string atomdir)
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
						string atomfile = atomdir + pdb + "_" + iter->first + ".atom.txt";
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
			printOneReportWithGeoDef(pdb, iter->first, geodir + "_" + iter->first + "_geo.txt");
			printAtomsForGeoDef(pdb, iter->first,atomdir + "_" + iter->first + "_atom.txt");
		}		
	}
}

void GeometricalDataReport::printOneReportWithGeoDef(PDBFile* pdb, string occupant, string fileName_geo)
{
	LogFile::getInstance()->writeMessage("Starting Geometric Data report for " + pdb->pdbCode);

	//stringstream report;	
	DataFrame geo_calcs(fileName_geo);
		
	geo_calcs.headerVector.push_back("PdbCode");
	geo_calcs.headerVector.push_back("Occupant");
	geo_calcs.headerVector.push_back("Chain");
	geo_calcs.headerVector.push_back("AminoNo");
	geo_calcs.headerVector.push_back("GeoAtoms");
	geo_calcs.headerVector.push_back("Alias");
	geo_calcs.headerVector.push_back("AminoCode");
	geo_calcs.headerVector.push_back("AminoNos");
	geo_calcs.headerVector.push_back("AminoCodes");
	geo_calcs.headerVector.push_back("AtomNos");
	geo_calcs.headerVector.push_back("SecStruct");
	geo_calcs.headerVector.push_back("GeoType");
	geo_calcs.headerVector.push_back("Value");

	//report << "PdbCode,Occupant,Chain,AminoNo,GeoAtoms,AminoCode,AminoNos,AminoCodes,AtomNos,SecStruct,GeoType,Value\n";
	
	
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
				geo_calcs.fileVector.push_back(getReportString(v1[i], pdb, occupant));				
				delete v1[i];
				v1[i] = nullptr;
			}

			for (unsigned int i = 0; i < v2.size(); ++i)
			{
				geo_calcs.fileVector.push_back(getReportString(v2[i], pdb, occupant));
				delete v2[i];
				v2[i] = nullptr;
			}

			for (unsigned int i = 0; i < v3.size(); ++i)
			{
				geo_calcs.fileVector.push_back(getReportString(v3[i], pdb, occupant));
				delete v3[i];
				v3[i] = nullptr;
			}

			for (unsigned int i = 0; i < v4.size(); ++i)
			{
				geo_calcs.fileVector.push_back(getReportString(v4[i], pdb, occupant));
				delete v4[i];
				v4[i] = nullptr;
			}

			for (unsigned int i = 0; i < v5.size(); ++i)
			{
				geo_calcs.fileVector.push_back(getReportString(v5[i], pdb, occupant));
				delete v5[i];
				v5[i] = nullptr;
			}

			for (unsigned int i = 0; i < v6.size(); ++i)
			{
				geo_calcs.fileVector.push_back(getReportString(v6[i], pdb, occupant));
				delete v6[i];
				v6[i] = nullptr;
			}
			
			for (unsigned int i = 0; i < v7.size(); ++i)
			{
				geo_calcs.fileVector.push_back(getReportString(v7[i], pdb, occupant));
				delete v7[i];
				v7[i] = nullptr;
			}
		}

	}

	geo_calcs.print();


	/*ofstream outfile(fileName1);
	if (outfile.is_open())
	{
		outfile << report.str();
	}*/

	
}


void GeometricalDataReport::printAtomsForGeoDef(PDBFile* pdb, string occupant, string fileName_atoms)
{
	LogFile::getInstance()->writeMessage("Starting Atoms Data report for " + pdb->pdbCode);

	DataFrame atom_file(fileName_atoms);
	
	atom_file.headerVector.push_back("pdb_code");
	atom_file.headerVector.push_back("atom_no");
	atom_file.headerVector.push_back("occupant");
	atom_file.headerVector.push_back("element");
	atom_file.headerVector.push_back("xcoord");
	atom_file.headerVector.push_back("ycoord");
	atom_file.headerVector.push_back("zcoord");
	atom_file.headerVector.push_back("occupancy");
	atom_file.headerVector.push_back("bfactor");

	vector<Atom*> atoms = ProteinManager::getInstance()->getAtoms(pdb->pdbCode, occupant);
	for (unsigned int a = 0; a < atoms.size(); ++a)
	{
		Atom* atm = atoms[a];
		if (atm)
		{
			if (atm->occupant == occupant)
				atom_file.fileVector.push_back(atm->getObservation());
		}
	}
	atom_file.print();
}

vector<string> GeometricalDataReport::getReportString(AtomGeo* ab, PDBFile* pdb, string occupant)
{
	//PdbCode, Occupant, Chain, AminoNo, GeoAtoms, AminoCode, AminoNos, AminoCodes, AtomNos, SecStruct, GeoType, Value
	vector<string> observation;	
	observation.push_back(pdb->pdbCode);
	observation.push_back(occupant);
	observation.push_back(ab->getChain());
	observation.push_back(quickInt(ab->getAminoId()));
	observation.push_back(ab->getGeoDef());
	observation.push_back(ab->getAlias());
	observation.push_back(ab->getAA());
	observation.push_back(ab->getAminoNos());
	observation.push_back(ab->getAACodes());
	observation.push_back(ab->getAtomNos());
	observation.push_back(ab->getSS());
	observation.push_back(ab->getGeoType());
	observation.push_back(quickRound(ab->getValue()));
	return observation;
}

string GeometricalDataReport::quickRound(double val)
{	 
	double dVal = val * 1000;
	int iVal = (int)round(dVal);
	dVal = (double)iVal / 1000;
	stringstream ss;
	ss << dVal;
	return ss.str();
}
string GeometricalDataReport::quickInt(int val)
{
	stringstream ss;
	ss << val;
	return ss.str();
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