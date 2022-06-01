from copy import deepcopy
import os
import re
import time
from shutil import copy
from random import choice
import numpy as np
import subprocess
from subprocess import run as externRun

import headerPoscar as hp

#Map from element(string) to element(number).
possibleElems = ["Ag", "Bi", 'I', "Vac"]
ELEM_NUM_MAP = dict()
for i in range(1, len(possibleElems) + 1):
    ELEM_NUM_MAP[possibleElems[i - 1]] = i

#Stoichiometry rules
STOICH = {"AgBiI4":  {"Ag": 1, "Bi": 1, 'I': 4},
          "Ag2BiI5": {"Ag": 2, "Bi": 1, 'I': 5},
          "Ag3BiI6": {"Ag": 3, "Bi": 1, 'I': 6},
          "AgBi2I7": {"Ag": 1, "Bi": 2, 'I': 7}}

#---Functions---#
def Fmt(s, prec=5):
    return f'{float(s):.{prec}f}'

def AlmostEqual(a, b, tol=0.00001):
    return abs(a - b) < tol
def Dist(atomA, atomB):
    return ((atomB.a - atomA.a)**2. + (atomB.b - atomA.b)**2. + (atomB.c - atomA.c)**2.)**(1./2.)

def Log(logLoc, message):
    with open(logLoc, 'a') as outfile:
        outfile.write(str(message))
        outfile.close()

def GetUniqueGenNums(runningLis, pendingLis, finishedLis, nReturns, a, b):
    def GetAllowedList(runLis, penLis, finLis, nMax, nMin=0):
        ret = []
        for i in range(nMin, nMax + 1):
            if ((i not in runLis) and (i not in penLis) and (i not in finLis)):
                ret.append(i)

        return ret

    allowed = GetAllowedList(runningLis, pendingLis, finishedLis, nMax=b, nMin=a)
    ##return entire list if requested number of returns >= possible returns
    if(len(allowed) <= nReturns):
        return allowed

    ret = []
    successfulAdditions = 0
    ##else generate random choices from allowed numbers
    while(True):
        trial = choice(allowed)
        if(trial not in ret):
            ret.append(trial)
            successfulAdditions += 1

        if(successfulAdditions == nReturns):
            return ret

#Returns a mapping of random numbers to stoich-following usable binary representations for POSCAR generation
def InitBinMap(binMapLoc, name):
    ret = dict()

    #file format is just a list of binary numbers
    if(name == "AgBiI4" or name == "AgBi2I7"):
        with open(binMapLoc, 'r') as infile:
            for num, line in enumerate(infile):
                ret[num] = str(line).replace("\n", '') ##keep as string to avoid losing prefixed zeros
        infile.close()
    #file format is eg 4 8 1 7 9999999 per line
    if(name == "Ag2BiI5" or name == "Ag3BiI6"):
        with open(binMapLoc, 'r') as infile:
            for num, line in enumerate(infile):
                ls = line.split()
                ret[num] = str(ls[0]) + '_' + str(ls[1]) + '_' + str(ls[2]) + '_' + str(ls[3])
        infile.close()

    return ret

def FileWordCount(infileLoc):
    cnt = 0
    with open(infileLoc, 'r') as infile:
        for li in infile:
            cnt += 1
        infile.close()
    return cnt

#Quick end of file
def Tail(filePath, nLines=10, buff=1028):
    if (buff < 1):
        print("Setting buffer size < 1 is a bad idea\n")
        return ["Tail():  Bad buffer size"]

    infile = open(filePath, 'r')
    foundLines = []

    blockCnt = -1
    # keep going until there are n lines
    while (len(foundLines) < nLines):
        try:
            infile.seek(blockCnt * buff, os.SEEK_END)
        except IOError:  # just in case you have a small file
            infile.seek(0)
            foundLines = infile.readlines()
            break

        foundLines = infile.readlines()
        blockCnt = blockCnt - 1
    infile.close()

    return foundLines[-nLines:]

def CheckOutcarForTime(filePath = "OUTCAR"):
    end = Tail(filePath, nLines = 32, buff = 4096)
    for words in end:
        if(re.search("Voluntary context switches", words)):
            return True
    return False

def AvgDiff(lisA, lisB):
    ret = 0.
    for i in range(0, len(lisA)):
        ret += abs(lisA[i] - lisB[i])

    return ret / len(lisA)

#Creates a directory capable of running a vasp job
def PrepVaspRun(pathTo, incarLoc, kptsLoc, potcarLoc, runLoc,
                poscarGenCommand, genNum):

    ##move necessary runfiles to new directory after making it
    if(not os.path.exists(pathTo)):
        os.makedirs(pathTo)

    ##remove files if they exist because copy can't overwrite files for some reason
    for file in ["INCAR", "KPOINTS", "POTCAR", "RUN_" + str(genNum)]:
        if(os.path.exists(os.path.join(pathTo, file))):
            os.remove(os.path.join(pathTo, file))

    copy(incarLoc, os.path.join(pathTo, "INCAR"))
    copy(kptsLoc, os.path.join(pathTo, "KPOINTS"))
    copy(potcarLoc, os.path.join(pathTo, "POTCAR"))
    copy(runLoc, os.path.join(pathTo, "RUN" + '_' + str(genNum)))

    ##generate a new poscar
    externRun(str(poscarGenCommand), shell=True) ##shell=True is dangerous - dont mess up :)

    return

#Submit a vasp run - assumes all necessary VASP files are in the folder moveTo
def SubVaspRun(cmd, moveTo, moveFrom):
    os.chdir(moveTo)

    #If POSCAR0 does not exist, this is the first attempt at a run
    if(not os.path.exists("POSCAR0")):
        copy("POSCAR", "POSCAR0")
    if(os.path.exists("CONTCAR")):
        if(FileWordCount("CONTCAR") > 5):
            copy("CONTCAR", "POSCAR")
    externRun(cmd, shell=True)

    os.chdir(moveFrom)
    return

def IsSuspended(runNum):
    cmd = "qstat -u vkb5066 | grep RUN_" + str(runNum)
    psuLine = externRun(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

    if(psuLine == '' or psuLine == " "):
        return False

    if('S' in psuLine.split()): ##psu will report an 'S' if the job is suspended.  Otherwise, a C or R or Q.
        return True
    return False

#Returns "finished" if vasp run is thought to be finished,
#        "needs resubbed" if vasp run is though to need resubmission,
#and     "running" if vasp run is thought to be running
def CheckVaspRun(outcarLoc, oszicarLoc, runNumber, secsBeforeConsideredFrozen=900):
    #if outcar and oszicar don't exist, the job hasn't even started yet (sitting in queue ie "running")
    if((not os.path.exists(outcarLoc)) and (not os.path.exists(oszicarLoc))):
        return "running"

    #if the job is suspended, consider it running - it will restart eventually, ... right? :)
    if(IsSuspended(runNum=runNumber)):
        return "running"

    #if outcar has timing information, job is finished
    if(CheckOutcarForTime(outcarLoc)):
        return "finished"

    #if outcar has no timing info, check oszicar...
    ##if oszicar has been updated recently, file is still running
    if(float(time.time()) - float(os.path.getmtime(oszicarLoc)) < secsBeforeConsideredFrozen):
        return "running"
    ##otherwise, the job is frozen
    return "needs resubbed"

#Checks if a poscar file has consistent stoichiometry
def PoscarMatchesStoich(pos, name):
    counts = {"Ag": 0, "Bi": 0, 'I': 0}

    for i in range(0, len(pos.atoms)):
        if(pos.atoms[i].atomType in ["Ag", "Bi", 'I']):
            counts[pos.atoms[i].atomType] += 1

    ok = True
    for elem, elemCount in STOICH[name].items():
        if(round(counts[elem] / min(counts.values())) != elemCount):
            ok = False

    return ok

#Makes a unique, consistent label given an atom
#Combines the atom coords into a string
def SitesToLabel(atom):
    return Fmt(atom.a) + '_' + Fmt(atom.b) + '_' + Fmt(atom.c)

#Returns the integer that corresponds to the element
def ElemToNum(elem):
    return ELEM_NUM_MAP[elem]

#Transforms a poscar atom to a (label, value) pair that can be used for ml
def SiteElemPair(atom):
    return SitesToLabel(atom), ElemToNum(atom.atomType)

#Returns the poscar of all site potential site locations (all initialized to Vac)
def InitAllSiteData(allVacLoc):
    p = hp.Poscar(allVacLoc)
    p.ConvertToCartesian()
    p.ConvertToDirect()

    for i in range(0, len(p.atoms)):
        p.atoms[i].atomType = "Vac"

    return p

#Returns a dict of (label, value) pairs.  Initializes any empty site to vacancy
def GetAllSiteElemPairs(allSitePoscar, thisPoscar, nTries = 100, tol=0.01):
    thisPoscar.ConvertToCartesian()
    thisPoscar.ConvertToDirect()

    ##allSitePoscar will have at least as many atoms as thisPoscar
    ok = False
    ret = dict()
    for n in range(0, nTries):
        ret.clear()

        ##Attempt to account for all sites, vacant or otherwise
        for i in range(0, len(allSitePoscar.atoms)):
            includedInThisPoscar = False
            for j in range(0, len(thisPoscar.atoms)):
                if(Dist(allSitePoscar.atoms[i], thisPoscar.atoms[j]) < tol):
                    key, val = SiteElemPair(thisPoscar.atoms[j])
                    ret[key] = val

                    includedInThisPoscar = True
                    break

            if(not includedInThisPoscar):
                key, val = SiteElemPair(allSitePoscar.atoms[i])
                ret[key] = val

        ##Tol modifier:
        ###If we successfully found all sites, finish up
        if(len(ret.keys()) == len(allSitePoscar.atoms)):
            ok = True
            break
        ###If we have too few sites, our tolerence was too tight - loosen it
        if(len(ret.keys()) < len(allSitePoscar.atoms)):
            tol += tol / float(nTries + 1)
            continue
        ###If we have too many sites, our tolerence was too loose - tighten it
        if(len(ret.keys()) < len(allSitePoscar.atoms)):
            tol -= tol / float(nTries + 1)
            continue

    if(ok):
        return ret
    else:
        print("GetAllSiteElemPairs(): Failed to account for all sites!\n")
        return False

#Returns a list of consistently sorted values given a dict of {labels, values}
#Make sure to call AFTER verifying that all expected sites are found, otherwise the labels will be out of
#order!
def GetSiteList(siteElemPairs):
    class helper:
        lab, val = None, None

        def __init__(self, la, va):
            self.lab = la
            self.val = va

    pairs = []
    for key, val in siteElemPairs.items():
        pairs.append(helper(la=key, va=val))
    pairs = sorted(pairs, key=lambda x: x.lab)

    return [p.val for p in pairs]

#Returns the final converged energy from OUTCAR.  No error checking :)
def GetFinEnergy(outcarLoc):
    ret = 0.
    with open(outcarLoc, 'r') as infile:
        for line in infile:
            if("TOTEN" in line.split()): ##will only run when outcar has time info so TOTEN -> fin conv nrg
                ret = float(line.split()[4].replace("eV", ''))
        infile.close()

    return ret

def DumpTrainingData(outfileLoc, Xlis, ylis, dec=8):
    if(len(Xlis) == 0):
        return

    head = ""
    for i in range(0, len(Xlis[0])):
        head += 'X' + str(i) + ' '*(dec - len(str(i)) + 2)
    head += 'y\n'

    with open(outfileLoc, 'w') as outfile:
        outfile.write(head)
        for i in range(0, len(Xlis)):
            thisStr = ""
            for j in range(0, len(Xlis[i])):
                thisStr += Fmt(Xlis[i][j], prec=dec) + ' '
            outfile.write(thisStr + Fmt(ylis[i]) + "\n")
        outfile.close()

    return
