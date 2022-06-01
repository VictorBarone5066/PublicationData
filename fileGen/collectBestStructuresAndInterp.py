import os
import shutil
import headerPoscar as hp
from numpy import polyfit
from random import uniform
from random import randint

VOL_LOC = "cubic-struct//outputVol.csv"
ALL_FILE_LOC = "cubic-struct//all//"
NEW_FILE_LOC = "cubic-final//s1//"

inds = dict()
#Read file
with open(VOL_LOC, 'r') as infile:
    for n, line_ in enumerate(infile):
        if(n == 0):
            continue
        line = line_.split(',')

        wdir = int(line[0].split('/')[6])
        nrg = float(line[2])
        vol = float(line[11])

        if(wdir not in inds.keys()):
            inds[wdir] = [[vol], [nrg], None] ##None = tmp variable for min e lat const
        else:
            inds[wdir][0].append(vol)
            inds[wdir][1].append(nrg)
    infile.close()

#Get min energy lattice constant
for ke, va in inds.items():
    x, y = zip(*sorted(zip(va[0], va[1])))
    poly = polyfit(x, y, 2)
    minLC1 = (-poly[1]/(2.*poly[0]))**(1./3.) ##3rd root of -b/2a
    poly = polyfit(x[1:], y[1:], 2) ##first point is sometimes slightly outside of the elastic region
    minLC2 = (-poly[1]/(2.*poly[0]))**(1./3.) ##3rd root of -b/2a
    va[2] = (minLC1 + minLC2)/2.

#Make new files (dont include in above loop in case more preprocessing is needed later)
for ke, va in inds.items():
    ##Find best base poscar (probably 1.0, but whatever)
    bstErr, bstFol = 10000000000000000, "0.000"
    for strain in ["0.94", "0.97", "1.0", "1.03"]:
        tstPos = hp.Poscar(ALL_FILE_LOC + str(ke) + "//" + strain + "//CONTCAR")
        tstLat = (tstPos.superCellVecA[0] + tstPos.superCellVecB[1] + tstPos.superCellVecC[2])/3.
        if(abs(tstLat - va[2]) < bstErr):
            bstErr = abs(tstLat - va[2])
            bstFol = strain
    lastContcar = hp.Poscar(ALL_FILE_LOC + str(ke) + "//" + bstFol + "//CONTCAR")

    ##Edit supercell vectors, adding random noise to break symmetry
    ##Noise should be between (min, max) noise for off-diagonal portions
    minNoise = min([lastContcar.superCellVecA[i] for i in [1, 2]] + \
                    [lastContcar.superCellVecB[i] for i in [0, 2]] + \
                    [lastContcar.superCellVecC[i] for i in [0, 1]])
    maxNoise = max([lastContcar.superCellVecA[i] for i in [1, 2]] + \
                    [lastContcar.superCellVecB[i] for i in [0, 2]] + \
                    [lastContcar.superCellVecC[i] for i in [0, 1]])

    noise = uniform(-abs(minNoise), abs(maxNoise))
    noise = -1.0*noise if randint(0, 1) == 0 else 1.0*noise ##actually, these if stmnts are redundent.
    lastContcar.superCellVecA[0] = va[2] + noise
    noise = uniform(-abs(minNoise), abs(maxNoise))
    noise = -1.0*noise if randint(0, 1) == 0 else 1.0*noise
    lastContcar.superCellVecB[1] = va[2] + noise
    noise = uniform(-abs(minNoise), abs(maxNoise))
    noise = -1.0*noise if randint(0, 1) == 0 else 1.0*noise
    lastContcar.superCellVecC[2] = va[2] + noise

    #Write this contcar to a new folder as a poscar
    os.makedirs(NEW_FILE_LOC + str(ke))
    lastContcar.Write(NEW_FILE_LOC + str(ke) + "//POSCAR")
