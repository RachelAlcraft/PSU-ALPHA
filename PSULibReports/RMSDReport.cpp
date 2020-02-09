#include "RMSDReport.h"
#include <sstream>
#include <fstream>


void RMSDReport::printReport(RMSD* rmsd, string fileName, bool optimised, string fileroot)
{
	stringstream report;
	report << "RMSD Calculation\n";
	report << rmsd->getAtomMatches() << "\n";
	report << "Value calculated = " << rmsd->calculateRMSD();

	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		outfile << report.str();
	}

	if (optimised)
	{
		//If this is an optimsied report we want to also print out the pdb files of the shifted structures
		rmsd->PDB1->printShiftedFile(fileroot);
		rmsd->PDB2->printShiftedFile(fileroot);
	}
}
