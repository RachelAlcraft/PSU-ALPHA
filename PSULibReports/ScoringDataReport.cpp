#include "ScoringDataReport.h"
#include <LogFile.h>
#include <ProteinManager.h>

void ScoringDataReport::printReport(PDBFile* pdb, string fileName)
{//produce data frame report for R reporting
	LogFile::getInstance()->writeMessage("Starting Scoring Data report for " + pdb->pdbCode);

	stringstream report;
	report << "Chain,AminoAcid,Id,SS,DataType,Atoms,Value\n";
	vector<AtomBond> bonds = ProteinManager::getInstance()->getAtomBonds(pdb->pdbCode);
	vector<AtomAngle> angles = ProteinManager::getInstance()->getAtomAngles(pdb->pdbCode);
	vector<AtomTorsion> torsions = ProteinManager::getInstance()->getAtomTorsions(pdb->pdbCode);

	for (unsigned int i = 0; i < bonds.size(); ++i)
	{
		report << bonds[i].getChain() << ",";
		report << bonds[i].getAA() << ",";
		report << bonds[i].getId() << ",";
		report << bonds[i].getSS() << ",";
		report << "BOND,";
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
		report << torsions[i].getAtoms() << ",";
		report << torsions[i].getValue() << "\n";
	}
	
	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		outfile << report.str();
	}
}