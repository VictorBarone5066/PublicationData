#Holds DOS information for a single atom.  If atomtype == full or atomNum == 0, corresponds to full dos
from copy import deepcopy
from scipy.signal import savgol_filter

class AtomDos:
    atomType = None
    atomNum = None

    eFermi = None

    energy = []
    dosInfo = {}
    summedDosInfo = {}

    hasF = None ##Is the l=3 quantum num written?
    hasISpin = None ##If ISPIN (spin up and down channels) are present
    hasNCL = None ##If noncollinear run (site projected channels) are present

    def __init__(self, data, atomType, atomNum, eFermi, hasF=False, hasISpin=False, hasNcl=False):
        self.atomType = None
        self.atomNum = None
        self.energy = []
        self.dosInfo = {}
        self.summedDosInfo = {}
        self.hasF = None
        self.hasISpin = None
        self.hasNCL = None
        self.eFermi = None

        self.atomType = atomType
        self.atomNum = atomNum
        self.eFermi = eFermi
        self.hasF = hasF
        self.hasISpin = hasISpin
        self.hasNCL = hasNcl

        #FULL DOS
        if(atomType == "full" and atomNum == 0):
            ##No Spin OR NCL
            if(not hasISpin):
                self.dosInfo = {"dos": []}
                for dat in data:
                    spl = dat.split()
                    self.energy.append(spl[0])
                    self.dosInfo["dos"].append(spl[1])
            ##Ispin
            if(hasISpin and not hasNcl):
                self.dosInfo = {"dosUp": [], "dosDn": []}
                for dat in data:
                    spl = dat.split()
                    self.energy.append(spl[0])
                    self.dosInfo["dosUp"].append(spl[1])
                    self.dosInfo["dosDn"].append(spl[2])

        #ATOM DOS
        if(atomType != "full" and atomNum != 0):
            ##ISPIN: No.  NCL: No. l=3: No
            if(not hasISpin and not hasNcl and not hasF):
                self.dosInfo = {"dosS": [], "dosP": [], "dosD": []}
                for dat in data:
                    spl = dat.split()
                    self.energy.append(spl[0])
                    self.dosInfo["dosS"].append(spl[1])
                    self.dosInfo["dosP"].append(spl[2])
                    self.dosInfo["dosD"].append(spl[3])
            ##ISPIN: No. NCL: No. l=3: Yes
            if(not hasISpin and not hasNcl and hasF):
                self.dosInfo = {"dosS": [], "dosP": [], "dosD": [], "dosF": []}
                for dat in data:
                    spl = dat.split()
                    self.energy.append(spl[0])
                    self.dosInfo["dosS"].append(spl[1])
                    self.dosInfo["dosP"].append(spl[2])
                    self.dosInfo["dosD"].append(spl[3])
                    self.dosInfo["dosF"].append(spl[4])
            ##ISPIN: Yes. NCL: No. l=3: No
            if(hasISpin and not hasNcl and not hasF):
                self.dosInfo = {"dosUpS": [], "dosDnS": [], "dosUpP": [], "dosDnP": [],
                                "dosUpD": [], "dosDnD": []}
                for dat in data:
                    spl = dat.split()
                    self.energy.append(spl[0])
                    self.dosInfo["dosUpS"].append(spl[1])
                    self.dosInfo["dosDnS"].append(spl[2])
                    self.dosInfo["dosUpP"].append(spl[3])
                    self.dosInfo["dosDnP"].append(spl[4])
                    self.dosInfo["dosUpD"].append(spl[5])
                    self.dosInfo["dosDnD"].append(spl[6])
            ##ISPIN: Yes. NCL: No. l=3: Yes
            if(hasISpin and not hasNcl and not hasF):
                self.dosInfo = {"dosUpS": [], "dosDnS": [], "dosUpP": [], "dosDnP": [],
                                "dosUpD": [], "dosDnD": [], "dosUpF": [], "dosDnF": []}
                for dat in data:
                    spl = dat.split()
                    self.energy.append(spl[0])
                    self.dosInfo["dosUpS"].append(spl[1])
                    self.dosInfo["dosDnS"].append(spl[2])
                    self.dosInfo["dosUpP"].append(spl[3])
                    self.dosInfo["dosDnP"].append(spl[4])
                    self.dosInfo["dosUpD"].append(spl[5])
                    self.dosInfo["dosDnD"].append(spl[6])
                    self.dosInfo["dosUpF"].append(spl[7])
                    self.dosInfo["dosDnF"].append(spl[8])
            ##ISPIN: No.  NCL: Yes. l=3: No
            if(not hasISpin and hasNcl and not hasF):
                self.dosInfo = {}


        #string to float
        self.energy = [float(e) for e in self.energy]
        for key, val in self.dosInfo.items():
            self.dosInfo[key] = [float(v) for v in val]

    #Scale energy to E' = E - Efermi
    def ScaleByEfermi(self):
        self.energy = [e - self.eFermi for e in self.energy]

    #Apply a smoothing filter to the DOS
    def ApplySavgol(self, winsize=11, poly=3):
        for key, val in self.dosInfo.items():
            self.dosInfo[key] = savgol_filter(val, window_length=winsize, polyorder=poly)
        if(self.atomType == "full"):
            return
        for key, val in self.summedDosInfo.items():
            self.summedDosInfo[key] = savgol_filter(val, window_length=winsize, polyorder=poly)

    #Creates useful (for plotting) sums of pdos
    def InitializeSums(self):
        ##FULL
        if(self.atomType == "full"):
            if(self.hasISpin == False):
                self.summedDosInfo = {"dosSum": deepcopy(self.dosInfo["dos"])}
            if(self.hasISpin == True):
                self.summedDosInfo = {"dosSum": []}
                for up, dn in zip(self.dosInfo["dosUp"], self.dosInfo["dosDn"]):
                    self.summedDosInfo["dosSum"].append(up + dn)
        ##ATOM
        if(self.atomType != "full"):
            ##ISPIN: No.  NCL: No. l=3: No |---> s + p + d
            if (not self.hasISpin and not self.hasNCL and not self.hasF):
                self.summedDosInfo = {"dosSum": []}
                for s, p, d in zip(self.dosInfo["dosS"], self.dosInfo["dosP"], self.dosInfo["dosD"]):
                    self.summedDosInfo["dosSum"].append(s + p + d)
            ##ISPIN: Yes. NCL: No. l=3: No |---> s(up/dn) + p(up/dn) + d(up/dn)
            if (self.hasISpin and not self.hasNCL and not self.hasF):
                self.summedDosInfo = {"dosSumUp": [], "dosSumDn": []}
                for su, pu, du in zip(self.dosInfo["dosUpS"], self.dosInfo["dosUpP"], self.dosInfo["dosUpD"]):
                    self.summedDosInfo["dosSumUp"].append(su + pu + du)
                for sd, pd, dd in zip(self.dosInfo["dosDnS"], self.dosInfo["dosDnP"], self.dosInfo["dosDnD"]):
                    self.summedDosInfo["dosSumDn"].append(sd + pd + dd)

    ##TODO: add spin gaps, atom-specific gaps, etc.  Returns a list of [vbm, cbm, bandgap] (if asked for)
    def GetGapInfo(self, tol=1e-8, retVbm=False, retCbm=False, retGap=True):
        if(self.atomType != "full"):
            return None
        ret = dict()

        if(len(list(self.summedDosInfo.keys())) == 0):
            self.InitializeSums()
        ##Get location (index) of eFermi
        eFermiLoc = 0
        bestDiff = abs(self.energy[0] - self.eFermi)
        for i in range(1, len(self.energy)):
            if(abs(self.energy[i] - self.eFermi) < bestDiff):
                bestDiff = abs(self.energy[i] - self.eFermi)
                eFermiLoc = i

        #Insulator check
        if(self.summedDosInfo["dosSum"][eFermiLoc] >= tol):
            print("GetGapInfo: No gap found.")
            ret["BG"] = 0.
            return ret

        vbm, cbm = self.summedDosInfo["dosSum"][-1], self.summedDosInfo["dosSum"][0]
        #Search backwards from eFermi to find VBM
        for i in range(eFermiLoc, 0, -1):
            if(self.summedDosInfo["dosSum"][i] >= tol):
                vbm = 1./2.*(self.energy[i] + self.energy[i+1])
                if(retVbm):
                    ret["vbm"] = vbm
                break
        #search forwards from eFermi to find CBM
        for i in range(eFermiLoc, len(self.energy), +1):
            if(self.summedDosInfo["dosSum"][i] >= tol):
                cbm = 1./2.*(self.energy[i] + self.energy[i-1])
                if(retCbm):
                    ret["cbm"] = cbm
                break
        if(retGap):
            ret["BG"] = cbm - vbm

        return ret



#returns a dict of DOSCAR header information
def DoscarInfo(inLoc, lineNum=5):
    ret = dict()
    with open(inLoc, 'r') as infile:
        for ln, lin in enumerate(infile):
            if(ln == lineNum):
                ret["eMax"] = float(lin.split()[0])
                ret["eMin"] = float(lin.split()[1])
                ret["nedos"] = int(lin.split()[2])
                ret["eFermi"] = float(lin.split()[3])
                break
        infile.close()
    return ret

#returns a 2D array, where each elem is a pair of [atom type, range], where range is [low, high] inclusive
def ScanPoscar(inLoc, atomTypeLineNum=5, atomNumLineNum=6):
    atomTypes, atomTypeNums, ret = [], [], []
    with open(inLoc, 'r') as infile:
        for i, line in enumerate(infile):
            if (i == atomTypeLineNum):
                for a in line.split():
                    atomTypes.append(a)
            if (i == atomNumLineNum):
                for a in line.split():
                    atomTypeNums.append(int(a))
                infile.close()
                break

    #Don't worry about how this works
    tot = sum(atomTypeNums)
    ret.append([atomTypes[0], [1, atomTypeNums[0]]])
    tot= tot - atomTypeNums[0]
    i = 1
    while(tot > 0):
        ret.append([atomTypes[i], [ret[i-1][1][1] + 1, ret[i-1][1][1] + atomTypeNums[i]]])
        tot = tot - atomTypeNums[i]
        i = i + 1

    return ret

#Returns a list of AtomDos objects.  The 0'th atom is the full DOS
def ApxEql(a, b, TOL=0.01):
    return (abs(a - b) <= TOL)
def ReadDoscar(doscarLoc, poscarLoc=None, iSpin=False, ncl=False, hasF=False):
    def IsAtomHead(dosInf, lin):
        return (ApxEql(dosInf["eMax"], float(lin.split()[0])) and \
                ApxEql(dosInf["eMin"], float(lin.split()[1])) and \
                ApxEql(dosInf["nedos"], int(lin.split()[2])) and \
                ApxEql(dosInf["eFermi"], float(lin.split()[3])) and \
                ApxEql(1.0, float(lin.split()[4])))

    #Set up atom ranges from POSCAR if available
    atomRanges = [["full", [0, 0]]]
    if(poscarLoc != None):
        ranges = ScanPoscar(poscarLoc)
        for r in ranges:
            atomRanges.append(r)

    #Get header info
    doscarInfo = DoscarInfo(doscarLoc)

    ret = []
    with open(doscarLoc, 'r') as infile:
        atomNum = -1
        dosInfo = []
        header = False
        for num, line in enumerate(infile):
            ##Determine if we should start reading new dos info
            if(num > 4 and len(line.split()) == 5):
                if(IsAtomHead(doscarInfo, line)):
                    header = True
            ##Skip any bullshit empty lines that sometimes appear for no reason
            if(len(line.split()) == 0):
                continue

            if(not header):
               dosInfo.append(line)

            #Set flag to begin reading dos info on the next loop
            if(header):
                ##Append info to AtomDos, add to list
                if(len(dosInfo) != 0):
                    for r in atomRanges:
                        ###r[0] = atomtype, r[1][0] = atomRangeMin, r[1][1] = atomRangeMax
                        if(r[1][0] <= atomNum <= r[1][1]):
                            ret.append(AtomDos(data=dosInfo, atomType=r[0], atomNum=atomNum,
                                               eFermi=doscarInfo["eFermi"]))
                            break
                atomNum += 1
                dosInfo = []
                header = False

        #Add last atom, finish up
        for r in atomRanges:
            if (r[1][0] <= atomNum <= r[1][1]):
                ret.append(AtomDos(data=dosInfo, atomType=r[0], atomNum=atomNum,
                                   eFermi=doscarInfo["eFermi"]))
        infile.close()

    return ret

#Sums the corresponding DOS info element-wise
#Takes a list of AtomDos objects, returns a dict of summed AtomDos objects with elements = keys
def SumElementDos(atomDosList):
    ret = dict()
    for atomDos in atomDosList:
        if(atomDos.atomType == "full"):
            continue
        if(atomDos.atomType not in ret.keys()):
            ret[atomDos.atomType] = deepcopy(atomDos)
            ret[atomDos.atomType].atomNum = -1
        else:
            for key, val in atomDos.dosInfo.items():
                for i, v in enumerate(val):
                    ret[atomDos.atomType].dosInfo[key][i] += v
            for key, val in atomDos.summedDosInfo.items():
                for i, v in enumerate(val):
                    ret[atomDos.atomType].summedDosInfo[key][i] += v
    return ret
