#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <string>
#include <vector>

#include "PoscarInfo.h"
#include "Calculations.h"

std::string POSCAR_LOC = "POSCAR_EXPLICIT_VAC";
std::string OUTFILE_LOC = "pairCorr";
int SCALE_TYPE = 0; //0 = no scale.  1 = normalize to number of that pair.  2 = sum of atoms.  3 = normalize to total sum of corr values
int PREC = 12; //number of decimals to write numbers as
double euler = 2.718281828;

//Values for correlation.  
double MAX_CORR = 1.0; //the maximum value that a corr pair can give.  Achieved by a pair dist of EXPECTED_NN_DIST
double EXPECTED_NN_DIST = 4.29000; //expected n.n. dist for any pair not including Iodine. Will give the max correlation val
double DELTA1 = 2.0 * EXPECTED_NN_DIST; double FRAC1 = 0.000025; //corr val of pair DELTA1 angst. from expected is FRAC1 of max val
double DELTA2 = 3.0 * EXPECTED_NN_DIST; double FRAC2 = 0.0000000005; //corr val of pair DELTA2 angst. from expected is FRAC2 of max val

//Pairs of interest.  Will be written to OUTFILE_LOC
const int NPAIRS = 8;
std::string POI[NPAIRS][2] = {{"Ag", "Ag"}, 
						      {"Ag", "Bi"},
							  {"Ag", "I"}, {"Ag", "Vac"},
						      {"Bi", "Bi"},
							  {"Bi", "I"}, {"Bi", "Vac"}, {"Vac", "Vac"}};


struct pairCorr
{
	std::string pair0;
	std::string pair1;

	double unscaledCorrVal;
	double scaledCorrVal;

	//Empty initialization
	pairCorr(std::string pair0_, std::string pair1_) {
		pair0 = pair0_;
		pair1 = pair1_;
		unscaledCorrVal = 0.0;
	}

	//Normal initialization
	pairCorr(atomPair pair, double ucv) {
		pair0 = pair.pairedAtoms[0].atomType;
		pair1 = pair.pairedAtoms[1].atomType;
		unscaledCorrVal = ucv;
	}

	void ScaleBy(double denom) {
		scaledCorrVal = unscaledCorrVal / denom;
	}

	//Operator declerations
	bool operator==(const pairCorr&) const;
};
//Compare by pair's atom identifier strings
bool pairCorr::operator ==(const pairCorr &paircorr) const {
	return ((pair0 == paircorr.pair0 && pair1 == paircorr.pair1) || (pair0 == paircorr.pair1 && pair1 == paircorr.pair0));
}

//Function Declerations
double GetPairCorrelationValue(double, double, double, double, double, double);
double GetNNPairDist(std::string, std::string);
std::string FloatToString(const double, int prec);
void WriteOutfile(std::vector <std::vector <pairCorr>>, std::string, int);

//cutoffVal is the fraction of maxVal that will be returned when scaledDist = cutoffDist (in angstroms)
double GetPairCorrelationValue(double scaledDist, double maxVal = MAX_CORR, double cutoffVal1 = FRAC1, double cutoffDist1 = DELTA1,
							   double cutoffVal2 = FRAC2, double cutoffDist2 = DELTA2) {
	double order = std::log((std::log(cutoffVal2)) / (std::log(cutoffVal1))) / std::log(cutoffDist2 / cutoffDist1); //log base b of a is log(a)/log(b)
	double bParam = cutoffDist1 * std::pow(std::log(1.0 / double(cutoffVal1)), -1.0 / double(order));
	return maxVal * std::pow(euler, -(std::pow(std::abs(scaledDist) / bParam, order)));
}

int main(int argc, char *argv[]) {
	//Command line arguments.  Input and output file locations
	std::string infileLoc = POSCAR_LOC;
	std::string outfileLoc = OUTFILE_LOC;
	for (int i = 0; i < argc; i++) {
		if (argv[i][0] == '-' && (argv[i][1] == 'i' || argv[i][2] == 'I')) ///infile location
			infileLoc = std::string(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'o' || argv[i][2] == 'O')) ///outfile location
			outfileLoc = std::string(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 's' || argv[i][2] == 'S')) ///scaling type
			SCALE_TYPE = std::stoi(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'p' || argv[i][2] == 'P')) ///output precision
			PREC = std::stoi(argv[i + 1]);
	}

	//Get all pairs
	Poscar p("readAll", infileLoc);
	p.convertToCartesian();
	p.fetchAtomPairs(18.0, 3); ///all atom pairs withing 18 angstroms of one another.  Second number = ceil(18 / smallest unit vector size)
	
	//Fill vector with unscaled instances of pair correlation values
	std::vector <std::vector <pairCorr>> pairCorrHolder;
	///Fill with initial values for each pair that I actually care about
	for (int i = 0; i < NPAIRS; i++) {
		pairCorrHolder.push_back(*new std::vector <pairCorr>);
		pairCorrHolder[i].push_back(*new pairCorr(POI[i][0], POI[i][1]));
	}
	///Now, fill the vector.  Also get the sum of corr values
	double corrValSum = 0;
	for (int i = 0; i < p.atomPairs.size(); i++) {
		double thisCorrVal = GetPairCorrelationValue(p.atomPairs[i].lenBetween - GetNNPairDist(p.atomPairs[i].pairedAtoms[0].atomType,
													 p.atomPairs[i].pairedAtoms[1].atomType));
		pairCorr thisPairCorr(p.atomPairs[i], thisCorrVal);

		for (int j = 0; j < pairCorrHolder.size(); j++) {
			if (thisPairCorr == pairCorrHolder[j][0]) {
				pairCorrHolder[j].push_back(thisPairCorr);
				break;
			}
		}
		corrValSum += thisCorrVal;
	}

	//Set the first index of each A-B set to its normalized value
	for (int i = 0; i < pairCorrHolder.size(); i++) {
		for (int j = 0; j < pairCorrHolder[i].size(); j++) {
			pairCorrHolder[i][0].unscaledCorrVal += pairCorrHolder[i][j].unscaledCorrVal;
		}
		if (SCALE_TYPE == 0) ///no scaling
			pairCorrHolder[i][0].ScaleBy(1.0);
		else if (SCALE_TYPE == 1) //scale by number of a-b pairs
			pairCorrHolder[i][0].ScaleBy(pairCorrHolder[i].size() - 1);
		else if (SCALE_TYPE == 2) { //scale by number a atoms + number of b atoms
			int nA = 0; int nB = 0;
			for (int k = 0; k < p.atomCoords.size(); k++) {
				if (p.atomCoords[k].atomType == pairCorrHolder[i][0].pair0)
					nA++;
				if (p.atomCoords[k].atomType == pairCorrHolder[i][0].pair1)
					nB++;
			}
			pairCorrHolder[i][0].ScaleBy(2.0*(double(nA) + double(nB)));
		}
		else if(SCALE_TYPE == 3) ///scale by total sum
			pairCorrHolder[i][0].ScaleBy(corrValSum);
	}

	//Give Results
	WriteOutfile(pairCorrHolder, outfileLoc, PREC);
	return 0;
}

double GetNNPairDist(std::string atomType1, std::string atomType2) {
	//Ag-Ag, Ag-Bi, Ag-Vac, Bi-Bi, Bi-Vac, Vac-Vac pairs 
	//if(atomType1 != "I" && atomType2 != "I")
		return EXPECTED_NN_DIST;
	//return -50.0; //for X-I bonds (will make correlation value zero)
}

std::string FloatToString(const double val, int prec = 10){
	std::ostringstream out;
	out.precision(prec);
	out << std::fixed << val;
	return out.str();
}

void WriteOutfile(std::vector <std::vector <pairCorr>> pairs, std::string outloc, int prec) {
	std::ofstream outfile(outloc);
	outfile << std::fixed << std::setprecision(prec);
	std::string lin = "CSV ";
	for (int i = 0; i < pairs.size(); i++) {
		outfile << "PAIRS " << pairs[i][0].pair0 << " " << pairs[i][0].pair1 << " " << pairs[i][0].scaledCorrVal << "\n";
		if (i != 0)
			lin += ',';
		lin += pairs[i][0].pair0 + ',' + pairs[i][0].pair1 + ',' + FloatToString(pairs[i][0].scaledCorrVal, prec);
	}
	outfile << lin << "\n" << "END\n";
	outfile.close();
}
