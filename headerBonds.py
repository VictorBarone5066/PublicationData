#-----FILE INPUT OUTPUT---------------------------------------------------------------------------------
#For Site Location Data
#Output file is written as:
#wDir, elem1elem2, count, avgDist, ...
class BondData:
    wDirectory = None
    allBondInfo = None

    def __init__(self, line, dataBeginOffset=1):
        self.wDirectory = None
        self.allBondInfo = dict()

        spl = line.split(',')
        self.wDirectory = int(spl[0])

        elemStr, count, dist = None, None, None
        for i in range(1, len(spl)):
            if((i-dataBeginOffset)%3 == 0):
                elemStr = spl[i]
            if((i-dataBeginOffset)%3 == 1):
                count = int(spl[i])
            if((i-dataBeginOffset)%3 == 2):
                dist = float(spl[i])
                self.allBondInfo[elemStr] = {'pairCount': count, 'avgDist': dist}


#For standard Output (from writeOutput.sh).  Assumes the bandgap is written after volume
class OutData:
    directory = None
    wDirectory = None
    nrg = None
    nConv = None
    irredKpts = None
    aV, bV, cV = None, None, None
    al, be, ga = None, None, None
    vol = None
    bg = None

    line = None

    def __init__(self, lin):
        # The fact that I have to do this is completly insane.
        self.directory = None
        self.wDirectory = None
        self.nrg = None
        self.nConv = None
        self.irredKpts = None
        self.aV, self.bV, self.cV = None, None, None
        self.al, self.be, self.ga = None, None, None
        self.vol = None
        self.bg = None
        self.line = None

        # Now to the actual initialization from csv file
        self.directory = lin.split(',')[0]
        self.wDirectory = int(lin.split(',')[1])
        self.nrg = float(lin.split(',')[2])
        self.nConv = int(lin.split(',')[3])
        self.irredKpts = int(lin.split(',')[4])
        self.aV, self.bV, self.cV = float(lin.split(',')[5]), float(lin.split(',')[6]), \
                                    float(lin.split(',')[7])
        self.al, self.be, self.ga = float(lin.split(',')[8]), float(lin.split(',')[9]), \
                                    float(lin.split(',')[10])
        self.vol = float(lin.split(',')[11])
        self.bg = float(lin.split(',')[12])

        self.line = lin

class PairedData(OutData):
    pairsIni = None
    pairsFin = None

    def __init__(self, line, pi, pf):
        OutData.__init__(self, line)
        self.pairsIni = pi
        self.pairsFin = pf

#Gets site data into a list
def GetBondData(infileLoc, dataBeginOffset=1):
    # Get pre-run info
    allSiteDat = []
    with open(infileLoc, 'r') as infile:
        for i, line in enumerate(infile):
            allSiteDat.append(BondData(line, dataBeginOffset=dataBeginOffset))
        infile.close()

    return allSiteDat

#Gets standard cell information into a list
def GetOutputData(infileLoc):
    # Get post-run info
    allOutDat = []
    with open(infileLoc, 'r') as infile:
        for lineNum, line_ in enumerate(infile):
            if (lineNum > 0):  ##don't read in the header
                allOutDat.append(OutData(line_))
    infile.close()

    return allOutDat

#Gets paired data as a list
def GetPairedInfo(outLoc, initPairLoc, finalPairLoc):
    ##Get outfile info and pair info for both initial and final structures
    out = GetOutputData(infileLoc=outLoc)
    init = GetBondData(infileLoc=initPairLoc)
    fina = GetBondData(infileLoc=finalPairLoc)

    #Put them in a dict to avoid O(n^3) garbage
    tmp = dict()
    for ou in out:
        tmp[ou.wDirectory] = []
    for elem in [out, init, fina]:
        for i in range(0, len(elem)):
            tmp[elem[i].wDirectory].append(elem[i]) ##order will be outfileInfo, initPairInfo, finPairInfo

    #Now return the list of paired info
    ret = []
    for va in tmp.values():
        ret.append(PairedData(va[0].line, va[1].allBondInfo, va[2].allBondInfo))
    return ret



#-----EVERYTHING ELSE---------------------------------------------------------------------------------
#Returns average pair length and pair counts as a dict of {pairType: [count, dist]}
def GetAvgDicts(allPairLis):
    ret = dict()

    for dat in allPairLis:
        for ke, va in dat.items():
            if(ke not in ret.keys()):
                ret[ke] = {"pairCount": va["pairCount"], "avgDist": va["avgDist"]}
            else:
                ret[ke]["pairCount"] = ret[ke]["pairCount"] + va["pairCount"]
                ret[ke]["avgDist"] = ret[ke]["avgDist"] + va["avgDist"]

    for va in ret.values():
        va["pairCount"] = va["pairCount"]/len(allPairLis)
        va["avgDist"] = va["avgDist"]/len(allPairLis)

    return ret
