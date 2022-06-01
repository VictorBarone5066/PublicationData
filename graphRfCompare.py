#!/usr/bin/python3.9

from sklearn.ensemble import RandomForestRegressor
import os
import numpy as np
from shutil import copy
from shutil import rmtree
from time import sleep
import numpy as np
from random import randint
from random import shuffle

import headerPoscar as hp
import headerRndmFrst as hrf
import globalFig as gf
join = os.path.join
rfr = RandomForestRegressor(random_state=468) ##random state for reproducability (but it dosn't work!!!)
					      ##numpy prob has its own random seed or something
from matplotlib import pyplot as plt

#Random Forest Constants
N_TOT = 12870
N_LOOPS = 5 ##stopping condition - if more than N_LOOPS are done, main loop ends
N_JOBS_PER_LOOP = 100 ##number of files choosen from potential jobs to run before further tree training
SPLIT_FRAC = 0.328 ##reduction of potential files to choose per loop (smaller = more reduction)
                   ## 0.328 = the exact number needed to train 500 points
N_REQ_MINS = 40 ##number of mins that we hope to find

RANGE_FROM_MIN = 0.4 ##eV - I'll print out the number of files prediced up to this number from the min energy

#Control Constants
SORT = True
PAIRCOUNTS = True
CUTBYSPLIT = False #SORT and PAIRCOUNTS

#File i/o constants
DATA_FOLDER = "agbii4-cubic-cryst-111-gga-rfCompare//"

def SettingsToCol(sett):
    if(sett == [True, False, False]): #No pair counts, no search space reduction
        return 0
    if(sett == [True, True, False]): #pair counts, no search space reduction
        return 1
    if(sett == [True, True, True]): #pair counts, search space reduction
        return 2

fig, ax = plt.subplots(1, 3, sharey=True, sharex=False, figsize=(gf.WIDTH(2), gf.HEIGHT(1)))

                                #Sort, paircnt, cutbysplit
for col, settings in enumerate([[True, False, False], [True, True, False], [True, True, True]]):
    SORT = settings[0]
    PAIRCOUNTS = settings[1]
    CUTBYSPLIT = settings[2]

    class Dat:
        x = None
        py = None
        ay = None

        def __init__(self, x_, ay_):
            self.x = x_
            self.ay = ay_

    def Fmt(s, d=5):
        return f'{float(s):.{d}f}'

    def GetTrainingData(fullList, boundL=0, boundH=0):
        return np.array([f.x for f in fullList[boundL:boundH]]), \
               np.array([f.ay for f in fullList[boundL:boundH]])

    ##prev* = lists.  new* = numpy arrays
    def AppendToTotalTrainingData(prevX, prevY, newX, newY):
        for i in range(0, len(list(newX))):
            prevX.append(newX[i])
            prevY.append(newY[i])

    def ReinitPreds(fullList):
        if(fullList == []):
            return

        predictions = list(rfr.predict(np.array([f.x for f in fullList])))
        for i in range(0, len(predictions)):
            fullList[i].py = predictions[i]

    def CutList(fullList, frac=0.0, indList=[], type='f'):
        if(type == 'f'):
            indMax = int(len(fullList) * frac)
            return fullList[:indMax]
        if(type == 'i'):
            ret = fullList[:]
            for i in indList:
                ret[i].py = 987654321 ##god help you if your predicted energy is this high
            return [r for r in ret if r.py != 987654321]

    #Make sure bestList, actList are the same size or this gives garbage
    def MinAcc(bestList, actList):
        acc = 0.
        for b in bestList:
            if(b.x in [a.x for a in actList]):
                acc += 1.
        return acc/len(bestList)

    def UpdateBestMinList(currMinList, allTestedList):
        atlS = sorted(allTestedList, key = lambda x: x.ay)
        currMinList = atlS[0:N_REQ_MINS]
        return currMinList

    def DoesMinExist(fullList, trueMin, tol=0.000001):
        for elem in fullList:
            if(elem.ay - trueMin <= tol):
                return True
        return False


    #Read all data
    print("reading bondInfo...")
    allBondInfo = dict()
    with open(DATA_FOLDER + "pairInfoIni.csv", 'r') as infile:
        for line in infile:
            spl = line.split(',')

            num = int(spl[0])
            agag = int(spl[2])
            agbi = int(spl[5])
            agi = int(spl[8])
            bibi = int(spl[11])
            bii = int(spl[14])
            ii = int(spl[17])

            if(PAIRCOUNTS):
                allBondInfo[num] = [agag, agbi, agi, bibi, bii, ii]
            else:
                allBondInfo[num] = [] ##uncomment to remove pair counts from training data

    print("reading trainingDat...")
    allDat = []
    with open(DATA_FOLDER + "trainingDat", 'r') as infile:
        for nu, line in enumerate(infile):
            if(nu == 0):
                continue
            lin = line.split()
            allDat.append(Dat(x_=[int(li.replace('0', '').replace('.', '')) for li in lin[:-1]] + \
                                  allBondInfo[nu - 1],
                              ay_=float(lin[-1])
                              )
                          )
    allDatCOPY = allDat[:]

    ##Get true min and list of mins
    sortedDat = sorted(allDat, key=lambda x: x.ay)
    trueMin = sortedDat[0].ay
    trueMins = sortedDat[:N_REQ_MINS]

    #Simulate DaQ loops
    bestMin = 987654321
    bestMins = allDat[0:N_REQ_MINS]
    tX, tY = [], []
    trainedList = []
    for rflCnt in range(0, N_LOOPS):
        nJobsToRun = min(N_JOBS_PER_LOOP, len(allDat))

        ##randomly run 100-or-so, fit, predict, remove them from list of potential jobs to run
        shuffle(allDat)
        tX_, tY_ = GetTrainingData(allDat, boundL=0, boundH=nJobsToRun)
        AppendToTotalTrainingData(tX, tY, tX_, tY_)
        rfr.fit(tX, tY)
        ReinitPreds(allDat)
        for i in range(0, nJobsToRun):
            trainedList.append(allDat[i])
        allDat = CutList(allDat, indList=[i for i in range(0, nJobsToRun)], type='i')

        ##Update minimum data
        bestMin = min(tY) if min(tY) < bestMin else bestMin
        bestMins = UpdateBestMinList(bestMins, trainedList)

        ##print training data
        """
        print("Loop " + str(rflCnt) + "...")
        print("   Num training pts:\t\t", len(tX))
        print("   Files considered:\t\t", len(allDat) + nJobsToRun)
        print("   Error in minimum:\t\t", Fmt(bestMin - trueMin))
        print("   Error in min", N_REQ_MINS, "files:\t", Fmt(1. - MinAcc(bestMins, trueMins)))
        print("   Min in remaining files?\t", DoesMinExist(allDat, trueMin=trueMin))
        """

        ##sort predicted lowest to highest, cut list
        if(len(allDat) > N_JOBS_PER_LOOP):
            if(SORT):
                allDat.sort(key=lambda x: x.py)
            if(CUTBYSPLIT):
                allDat = CutList(allDat, frac=SPLIT_FRAC, type='f')
            #print("   List cut to length:\t\t", len(allDat))
        else:
        ##break condition: If we are near the end, run all remaining jobs
            tX_, tY_ = GetTrainingData(allDat, boundL=0, boundH=nJobsToRun)
            AppendToTotalTrainingData(tX, tY, tX_, tY_)
            rfr.fit(tX, tY)
            ReinitPreds(allDat)
            break

    #Then go back and estimate energies of ALL models with the training data
    print("Applying training data to all models...")
    x = [adc.ay for adc in allDatCOPY]
    y = list(rfr.predict(np.array([adc.x for adc in allDatCOPY])))

    #Plot things
    index = SettingsToCol(settings)

    letter = ":("
    if(index == 0):
        letter = "a)"
    if(index == 1):
        letter = "b)"
    if(index == 2):
        letter = "c)"
    ax[index].text(0.03, 0.95, letter, fontsize=gf.FONTSIZE_TEXT, transform=ax[index].transAxes)

    if(index != 0):
        ax[index] = gf.RmAxTicks(ax[index], x=False, y=True)
    ax[index] = gf.SetTickSize(ax[index])

    ax[index].set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax[index].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax[index].grid(alpha=gf.GRID_ALPHA)

    minEnergy = min(x)
    ax[index].plot([x_ - minEnergy for x_ in x], [y_ - minEnergy for y_ in y],
                 color="black", marker='o', linewidth=0.0, markersize=gf.MKSIZE_SML)

    limMax, limMin = max([x_ - minEnergy for x_ in x] + [y_ - minEnergy for y_ in y]), \
                     min([x_ - minEnergy for x_ in x] + [y_ - minEnergy for y_ in y])
    ax[index].plot([limMin, limMax], [limMin, limMax], color=gf.RUD_BI, linewidth=gf.LNSIZE_MED)
    print("num training points:", len(tX))

    #Get number of points between 0 and RANGE FROM MIN ev above the minimum energy
    #Also calculate the average error within this range
    counter = 0
    error = 0.
    for ptx, pty in zip([x_ - minEnergy for x_ in x], [y_ - minEnergy for y_ in y]):
        if(pty <= RANGE_FROM_MIN):
            error += abs(ptx - pty)
            counter += 1
    #print("num data points up to", RANGE_FROM_MIN, "ev from min:", counter)
    #print("avg error in this range:", error/counter, "\n")

ax[1].set_xlabel("DFT Energies (eV)", fontsize=12, fontweight="bold")
ax[0].set_ylabel("Predicted Energies (eV)", fontsize=12, fontweight="bold")
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.13, top=0.97, wspace=0.05, hspace=0.00)

plt.savefig("rfCompare.png", dpi=gf.DPI)
plt.savefig("rfCompare.pdf")
