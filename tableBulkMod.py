#Outputs bulk modulus
from numpy import polyfit

import headerBoltzmann as hb
from globalFig import FloatFormat as ff

BULK_DAT_IN = "agbii4-cubic-cryst-111-333-gga-bulkMod//output333.csv"
BULK_DAT_OUT = "bulkMods.csv"
EVperACUBED_TO_GPA = 160.2176487

allDatB = hb.GetOutputData("agbii4-cubic-cryst-444-gga-boltzmann//output.csv")
hb.SetAllProbs(datLis=allDatB, temp=600)

def VOfEMin(vols, nrgs, retE=False):
    eMin, eMinInd = nrgs[0], 0
    for i in range(1, len(nrgs)):
        if(nrgs[i] < eMin):
            eMin = nrgs[i]
            eMinInd = i
    if(retE):
        return vols[eMinInd], eMin
    return vols[eMinInd]

def BulkMod(vols, nrgs):
    V0 = VOfEMin(vols, nrgs)
    coeffA = polyfit(x=vols, y=nrgs, deg=2)[0]
    return 2.*V0*coeffA*EVperACUBED_TO_GPA


allDat = dict() ##allDat[idNum] = [[volumes], [energies]]
with open(BULK_DAT_IN, 'r') as infile:
    for n, line in enumerate(infile):
        if(n == 0):
            continue
        lin = line.split(',')

        key = int(lin[0].split('/')[-2])
        vol = float(lin[11]) ##A^3
        nrg = float(lin[2]) ##eV

        if(key not in allDat.keys()):
            allDat[key] = [[vol], [nrg]]
        else:
            allDat[key][0].append(vol)
            allDat[key][1].append(nrg)
    infile.close()
#Sort, clean data
for ke, va in allDat.items():
    volsS, nrgsS = (list(t) for t in zip(*sorted(zip(va[0], va[1]))))
    if("output111.csv" == BULK_DAT_IN.split("//")[-1]):
        pass
    elif("output333.csv" == BULK_DAT_IN.split("//")[-1]):
        nrgsS = nrgsS[1:-1] ##remove things outside of elastic limit
        volsS = volsS[1:-1]
    else:
        print(":(")
        exit(1)

    allDat[ke][0] = volsS
    allDat[ke][1] = nrgsS


#Individual
with open(BULK_DAT_OUT, 'w') as outfile:
    outfile.write("ID,nrgMin,B(GPa)\n")
    for ke, va in allDat.items():
        #print(va, "\n")
        outfile.write(str(ke) + ',' + str(VOfEMin(va[0], va[1], retE=True)[1]) + ',' + \
                      str(BulkMod(va[0], va[1])) + "\n")
    outfile.close()

#Boltzmann average
sum = 0.
for dat in allDatB:
    sum += dat.prob * BulkMod(allDat[int(dat.wDirectory)][0], allDat[int(dat.wDirectory)][1])
print(sum)
