#include <iostream>
#include <fstream>
#include <string>

#include "PoscarInfo.h"

double TOL = 0.01; //how close (fractional) two coords must be to be considered the same site
int NTRIES = 1000; //how many tries the TOL loop goes through before quitting

std::string POSCAR_LOC = "POSCAR0";
std::string VAC_LOC = "POSCAR_VACS";
std::string OUTFILE_LOC = "POSCAR_EXPLICIT_VAC";

int main(int argc, char* argv[]) {
	//Command line arguments.  Input and output file locations
	std::string infileLoc = POSCAR_LOC;
	std::string vacLoc = VAC_LOC;
	std::string outfileLoc = OUTFILE_LOC;
	for (int i = 0; i < argc; i++) {
		if (argv[i][0] == '-' && (argv[i][1] == 'i' || argv[i][2] == 'I')) ///infile location
			infileLoc = std::string(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'v' || argv[i][2] == 'V')) ///vacancy infile location
			vacLoc = std::string(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'o' || argv[i][2] == 'O')) ///outfile location
			outfileLoc = std::string(argv[i + 1]);
	}

	//Get all potential vacancy locations.  A dummy POSCAR file with only vacancies as it's atoms (labeled "Vac")
	Poscar pVac("readAll", vacLoc);
	pVac.convertToDirect();

	Poscar write;
	double dTol = TOL / 10.0;
	int n = 0;
	for (; n < NTRIES; n++) {
		//Get all original POSCAR's atom locations
		Poscar p("readAll", infileLoc);
		p.convertToDirect();

		//For each potential site (vacancies), if  there is no actual atom present, add a vacancy.  
		for (int i = 0; i < pVac.atomCoords.size(); i++) {
			bool alreadyExists = false;
			for (int j = 0; j < p.atomCoords.size(); j++) {
				if (dist_(pVac.atomCoords[i], p.atomCoords[j]) < TOL) {
					alreadyExists = true;
					break;
				}
			}
			if (!alreadyExists)
				p.atomCoords.push_back(pVac.atomCoords[i]);
		}

		//Check to make sure that there as many atoms in the file to write as there are in the "all potential locations" file.
		int notIodineSites = 0;
		for (int i = 0; i < p.atomCoords.size(); i++) {
			if (p.atomCoords[i].atomType != "nothing - I want to count iodine sites now lol")
				notIodineSites++;
		}
		if (notIodineSites == pVac.atomCoords.size()) {
			write = p;
			break;
		}
		//If there are more atoms than expected, TOL is too lax.  Make it smaller.  
		if(p.atomCoords.size() > pVac.atomCoords.size())
			TOL -= dTol;
		//Otherwise, there are less atoms than expected - TOL is too strict.  Make it bigger.
		else
			TOL += dTol;
		
		//Edit dTol to make finer and finer adjustments to converge to solution
		dTol *= dTol * 10.0 / (10.0 + double(n));
	}

	if (n >= NTRIES - 1) {
		std::cout << "WriteExplicit2x2x1Vacancies.cpp: Correct number of atoms could not be found!\n";
		return 1;
	}

	//Write new POSCAR file with vacancies in preperation for next program.  End
	write.updateAll();
	write.write(outfileLoc);
	return 0;
}
