#Writes a file of normalized (between 0 and 1) intensities for each genNum
#Format is like this:
#genNum, theta0, theta1, theta2, ... thetaN
#num1  , int1_1, int1_2, ...
#num2  , int2_1, int2_2, ...
#and N is the same for each file.

from os import remove
from subprocess import run

import headerBoltzmann as hb

#Get the list of genNums to run
allDat = hb.GetOutputData("agbii4-cubic-cryst-444-gga-boltzmann//output.csv")

#Write the output file in-place
with open("xrd//allXrdGraphs.csv", 'w') as outfile:
    for n, dat in enumerate(allDat):
        print(n + 1, '/', len(allDat))

        genNum = dat.wDirectory

        ##Generate a file of theta, intensity where theta is evenly spaced and intensity ranges from 0 to 1
        ##for this genNum.  Read it in and prep for output to main output file
        run(".\\xrd_ser.exe -p agbii4-cubic-cryst-444-gga-contcars\CONTCAR-" + str(genNum) + \
            " -a xrd/AFF.csv -hkl 30 -wvln 1.5406 -spacing 0.01", shell=True)
        #Note a wavelength of 1.4506 angstroms = Cu K alpha radiation

        ##If this is the first line, we actually want to write the header and two-thetas at the top
        if(n == 0):
            twoThetaStr = ''
            with open("plt.csv", 'r') as infile:
                for line in infile:
                    twoThetaStr += ',' + line.split(',')[0].replace("\n", '')
                infile.close()
            outfile.write("genNum" + twoThetaStr + "\n")

        intStr = ''
        with open("plt.csv", 'r') as infile:
            for line in infile:
                intStr += ',' + line.split(',')[1].replace("\n", '')
            infile.close()

        #Write results to main output file
        outfile.write(str(genNum) + intStr + "\n")

    outfile.close()

#Clean up
remove("plt.csv")
remove("ints.csv")
