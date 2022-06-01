import headerBoltzmann as hb


allDat = hb.GetOutputData("agbii4-cubic-cryst-444-gga-boltzmann//output.csv")
hb.SetAllProbs(datLis=allDat, temp=600)

#Get dict of {genNum: [bondLength Ag-I, bondlength Bi-I]}
bondInfo = dict()
with open("agbii4-cubic-cryst-444-gga-bonds//bondInfo.csv", 'r') as infile:
    for line in infile:
        lin = line.split(',')
        bondInfo[int(lin[0])] = [float(lin[3]), float(lin[6])]

    infile.close()

#Get average bond lengths
agi, bii = 0., 0.
for dat in allDat:
    agi += float(bondInfo[dat.wDirectory][0]) * dat.prob
    bii += float(bondInfo[dat.wDirectory][1]) * dat.prob

#Get min bond lengths
mins = sorted(allDat, key=lambda x: x.nrg)[0]
agiM = float(bondInfo[mins.wDirectory][0])
biiM = float(bondInfo[mins.wDirectory][1])

#Print results
with open("bonds.csv", 'w') as outfile:
    outfile.write("Avg Ag-I," + str(agi) + "\n")
    outfile.write("Avg Bi-I," + str(bii) + "\n")
    outfile.write("Min Ag-I," + str(agiM) + "\n")
    outfile.write("Min Bi-I," + str(biiM) + "\n")
    outfile.close()
