
INFILE_LOC = "vxy"

allKs, allNrgs, allNrgsS = [], [], []

with open(INFILE_LOC, 'r') as infile:
    ks, nrgs, nrgsS = [], [], []
    for lin in infile:
        if("BEG" in lin):
            ks, nrgs, nrgsS = [], [], []
            continue
        if("END" in lin):
            allKs.append(ks)
            allNrgs.append(nrgs)
            allNrgsS.append(nrgsS)
            continue
        line = lin.split()

        ks.append(float(line[0]))
        nrgs.append(float(line[1]))
        nrgsS.append(float(line[2]))

    infile.close()


for i in range(0, len(allKs)):
    allKs[i], allNrgs[i], allNrgsS[i] = (list(t) for t in zip(*sorted(zip(allKs[i], allNrgs[i], allNrgsS[i]))))

from matplotlib import pyplot as plt
plt.title(INFILE_LOC)
for i in range(0, len(allKs)):
    plt.plot(allKs[i], allNrgs[i], 'r.')
    plt.plot(allKs[i], allNrgsS[i], 'k')


plt.show()
