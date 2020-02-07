#include "RMSDReport.h"
#include <sstream>
#include <fstream>


void RMSDReport::printReport(RMSD* rmsd, string fileName)
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
}
