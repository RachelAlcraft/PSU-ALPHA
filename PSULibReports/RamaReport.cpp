
#include "RamaReport.h"
#include <LogFile.h>
#include <ProteinManager.h>



void RamaReport::printReport(PDBFile* pdb, string fileName)
{//produce data frame report for R reporting
	LogFile::getInstance()->writeMessage("Starting Torsion report");

	stringstream report;
	report << "Chain,AminoAcid,Id,Hydrophobocity,Hydropathy,Volume,Donicity,Chemical,Physio,Charge,Polar,Phi,Psi,Omega,Chi1,Chi2,Chi3,Chi4,Chi5\n";
	map<string, Chain*> chains = ProteinManager::getInstance()->getChains(pdb->pdbCode);
	for (map<string, Chain*>::iterator iter = chains.begin(); iter != chains.end(); ++iter)
	{
		Chain* chain = iter->second;
		string chainid = iter->first;
		map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode,chainid);
		for (map<int, AminoAcid*>::iterator iter = aminos.begin(); iter != aminos.end(); ++iter)
		{
			double phi = 0;
			double psi = 0;
			double omega = 0;
			double chi1 = 0;
			double chi2 = 0;
			double chi3 = 0;
			double chi4 = 0;
			double chi5 = 0;

			AminoAcid* aa = iter->second;
			BackboneTorsion* torsion = aa->getBackboneTorsion();
			SidechainTorsion* sidetorsion = aa->getSidechainTorsion();
			if (sidetorsion)
			{
				chi1 = sidetorsion->getChi1();
				chi2 = sidetorsion->getChi2();
				chi3 = sidetorsion->getChi3();
				chi4 = sidetorsion->getChi4();
				chi5 = sidetorsion->getChi5();
			}
			if (torsion)
			{
				double phi = torsion->getPhi();
				double psi = torsion->getPsi();
				double omega = torsion->getOmega();

				report << chainid << "," << aa->AminoCode << "," << aa->aminoId << ",";
				report << aa->Hydro << ",";
				report << aa->Hydropathy << ",";
				report << aa->Volume << ",";
				report << aa->Donicity << ",";
				report << aa->Chemical << ",";
				report << aa->Physio << ",";
				report << aa->Charge << ",";
				report << aa->Polar << ",";
				report << phi << "," << psi << "," << omega << ",";
				report << chi1 << "," << chi2 << "," << chi3 << "," << chi4 << "," << chi5 << "\n";
			}
		}
	}
	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		outfile << report.str();
	}
}

