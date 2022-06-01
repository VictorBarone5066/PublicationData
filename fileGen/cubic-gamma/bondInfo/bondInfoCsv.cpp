#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "PoscarInfo.h"
#include "Calculations.h"

std::string INFILE_LOC = "POSCAR";
std::string BOND_FILE_LOC = "bondDataInfo";
std::string OUTFILE_LOC = "bondInfo.csv";

void PrintBondInfo(std::vector <atomPair>, std::string, std::string);


int main(int argc, char *argv[])
{
	for (int i = 0; i < argc; i++) {
		if (argv[i][0] == '-' && (argv[i][1] == 'i' || argv[i][1] == 'I'))
			INFILE_LOC = argv[i + 1];
		if (argv[i][0] == '-' && (argv[i][1] == 'b' || argv[i][1] == 'B'))
			BOND_FILE_LOC = argv[i + 1];
		if (argv[i][0] == '-' && (argv[i][1] == 'o' || argv[i][1] == 'O'))
			OUTFILE_LOC = argv[i + 1];
	}

	Poscar pos("readAll", INFILE_LOC.c_str());

	//Output info about bond lengths
	pos.fetchAtomBonds(BOND_FILE_LOC.c_str());
	pos.atomBonds.pop_back();
	PrintBondInfo(pos.atomBonds, OUTFILE_LOC.c_str(), pos.fileTitle);

	return 0;
}

void PrintBondInfo(std::vector <atomPair> pairs, std::string outfileLoc, std::string beg)
{
	//Output current bond Information
	std::ofstream outfile(outfileLoc.c_str(), std::ios::app);
	outfile << beg;

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
				outfile << ',' << pairs[i].type << ',' << count << ',' << sumLen / double(count);
		}
	}
	outfile << "\n";
	outfile.close();
}
