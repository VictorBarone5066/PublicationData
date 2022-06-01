#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "PoscarInfo.h"
#include "Calculations.h"

std::string INFILE_LOC = "C://Users//baron//Desktop//CONTCAR_cubii4WithAgbii4Positions.vasp";
std::string BOND_FILE_LOC = "C://Users//baron//Desktop//bondInfo";

void PrintBondInfo(std::vector <atomPair>);

//argv[0] = name of exe, argv[1] = infile loc, argv[2] = outfile loc, argv[3] = SiSi bond dist from infile, argv[4] = multiplication factor
int main(int argc, char *argv[])
{
	Poscar pos("readAll", INFILE_LOC.c_str());

	std::cout << "Input File Information:\nAtom Count Info...\n";

	//Output info about number of atoms
	for (int i = 0; i < pos.atomTypes.size(); i++)
		std::cout << pos.atomTypes[i] << "\t" << pos.atomTypeNums[i] << "\n";

	//Output info about bond lengths
	std::cout << "\nBond Length Info...\n";
	pos.fetchAtomBonds(BOND_FILE_LOC.c_str());
	pos.atomBonds.pop_back();
	PrintBondInfo(pos.atomBonds);

	std::cout << "\nDone.  Press a key to continue.\n";
	std::cin.get();
	return 0;
}

void PrintBondInfo(std::vector <atomPair> pairs)
{
	//Output current bond Information
	for (int i = 0; i < pairs.size(); i++)
	{
		if (pairs[i].extraInfo != "counted")
		{
			int count = 0;
			double sumLen = 0;
			static int k = 0;
			std::vector <atomPair> tmp;
			for (int j = 0; j < pairs.size(); j++)
			{
				if ((pairs[i] ^= pairs[j]) && pairs[j].extraInfo != "counted")
				{
					pairs[j].extraInfo = "counted";
					count++;
					sumLen = sumLen + pairs[j].lenBetween;
					tmp.push_back(pairs[j]);
				}
			}
			k++;
			if (count != 0)
				std::cout << "(" << k << ") " << pairs[i].type << " :\t" << count << "\tAvg Len:\t" << sumLen / double(count) << "\n";
		}
	}
}