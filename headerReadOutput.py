
#-----FILE INPUT OUTPUT---------------------------------------------------------------------------------
#For standard Output (from writeOutput.sh).  Assumes the bandgap is written after volume
class Dat:
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

    unnormalizedProb = None
    prob = None

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
        self.unnormalizedProb = None
        self.prob = None

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



#Gets standard cell information into a list
def GetOutputData(infileLoc):
    # Get post-run info
    allOutDat = []
    with open(infileLoc, 'r') as infile:
        for lineNum, line_ in enumerate(infile):
            if (lineNum > 0):  ##don't read in the header
                allOutDat.append(Dat(line_))
    infile.close()

    return allOutDat
