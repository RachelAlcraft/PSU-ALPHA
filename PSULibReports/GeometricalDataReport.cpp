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
						printOneReport(pf,iter->first, geofile);
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
			printOneReport(pdb, iter->first, fileName1 + "_" + iter->first + "_geo.txt");
		}		
	}
}

void GeometricalDataReport::printOneReport(PDBFile* pdb, string occupant, string fileName1)//, string fileName2)
{//produce data frame report for R reporting

	
	LogFile::getInstance()->writeMessage("Starting Scoring Data report for " + pdb->pdbCode);

	stringstream report;
	report << "Chain,AminoAcid,Id,SecStruct,DataType,ExperimentalMethod,Atoms,Value\n";
	vector<AtomBond> bonds = ProteinManager::getInstance()->getAtomBonds(pdb->pdbCode, occupant);
	vector<AtomAngle> angles = ProteinManager::getInstance()->getAtomAngles(pdb->pdbCode, occupant);
	vector<AtomTorsion> torsions = ProteinManager::getInstance()->getAtomTorsions(pdb->pdbCode, occupant);

	for (unsigned int i = 0; i < bonds.size(); ++i)
	{
		report << bonds[i].getChain() << ",";
		report << bonds[i].getAA() << ",";
		report << bonds[i].getId() << ",";
		report << bonds[i].getSS() << ",";
		report << "BOND,";
		report << pdb->experimentalMethod << ",";
		report << bonds[i].getAtoms() << ",";
		report << bonds[i].getValue() << "\n";

	}
	for (unsigned int i = 0; i < angles.size(); ++i)
	{
		report << angles[i].getChain() << ",";
		report << angles[i].getAA() << ",";
		report << angles[i].getId() << ",";
		report << angles[i].getSS() << ",";
		report << "ANGLE,";
		report << pdb->experimentalMethod << ",";
		report << angles[i].getAtoms() << ",";
		report << angles[i].getValue() << "\n";

	}
	for (unsigned int i = 0; i < torsions.size(); ++i)
	{
		report << torsions[i].getChain() << ",";
		report << torsions[i].getAA() << ",";
		report << torsions[i].getId() << ",";
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