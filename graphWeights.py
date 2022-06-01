#Outputs boltzmann-averaged things such as lattice constants, band gap

from matplotlib import pyplot as plt
import headerBoltzmann as hb
import globalFig as gf

allDat = hb.GetOutputData("agbii4-cubic-cryst-444-gga-boltzmann//output.csv")
hb.SetAllProbs(datLis=allDat, temp=600)
allDat.sort(key=lambda x: x.nrg)

probLo, probHi = [], []
for n, dat in enumerate(allDat):
    if(dat.prob*100. > 5):
        probLo.append(dat.prob)
    else:
        probHi.append(dat.prob)
    print(n+1, gf.FloatFormat(dat.prob*100., 2))

print("\nLow-energy sum:", gf.FloatFormat(float(sum(probLo)), 3) + " from " + str(len(probLo)) + " models",
      "\nHi-energy sum:", gf.FloatFormat(float(sum(probHi)), 3) + " from " + str(len(probHi)) + " models\n")

plt.plot([i for i in range(0, len(allDat))], [dat.prob*100. for dat in allDat], 'o')
plt.show()
