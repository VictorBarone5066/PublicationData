#Outputs boltzmann-averaged things such as lattice constants, band gap

import headerBoltzmann as hb
from globalFig import FloatFormat as ff

allDat = hb.GetOutputData("agbii4-cubic-cryst-444-gga-boltzmann//output.csv")
hb.SetAllProbs(datLis=allDat, temp=600)

params = {"aAvg": 0.0, "bAvg": 0.0, "cAvg": 0.0, "bgAvg": 0.0,
          "aMin": 0.0, "bMin": 0.0, "cMin": 0.0, "bgMin": 0.0}

mins = sorted(allDat, key=lambda x: x.nrg)[0]
params["aMin"] = mins.aV
params["bMin"] = mins.bV
params["cMin"] = mins.cV
params["bgMin"] = mins.bg

for dat in allDat:
    params["aAvg"] += dat.aV * dat.prob
    params["bAvg"] += dat.bV * dat.prob
    params["cAvg"] += dat.cV * dat.prob
    params["bgAvg"] += dat.bg * dat.prob

with open("params.csv", 'w') as outfile:
    outfile.write("-,a,b,c,bg\n")
    outfile.write("avg," + ff(params["aAvg"]) + ',' + ff(params["bAvg"]) + ',' + ff(params["cAvg"]) + ',' +\
                  ff(params["bgAvg"]) + "\n")
    outfile.write("min," + ff(params["aMin"]) + ',' + ff(params["bMin"]) + ',' + ff(params["cMin"]) + ',' +\
                  ff(params["bgMin"]) + "\n")
    outfile.close()
