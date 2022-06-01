from os import remove
from subprocess import run
from matplotlib import pyplot as plt

import headerBoltzmann as hb
import globalFig as gf

BOLTZMANN_WEIGHT_LOC = "agbii4-cubic-cryst-444-gga-boltzmann//output.csv"
TEMPERATURE = 600*5000 #kelvin
SAVE_LOC = "xrd"

#Get the dict of {GenNum: probability}
allDat = hb.GetOutputData(BOLTZMANN_WEIGHT_LOC)
hb.SetAllProbs(datLis=allDat, temp=TEMPERATURE)
probs = dict()
for dat in allDat:
    probs[str(dat.wDirectory)] = dat.prob

#Read twoThetas and intensities, apply weights to intensities
twoThetas, intensities = [], []
maxTwoThetaIndex = 0
maxTwoTheta = 0. ##the smallest value of twoTheta that we reached (different structures allow for different
                 ##evaluation limits.  This is the highest angle we can plot
with open("xrd//allXrdGraphs.csv", 'r') as infile:
    for n, line in enumerate(infile):
        if(n == 0): ###read two-thetas, initialize intensities to 0
            for lin in line.split(',')[1:]:
                twoThetas.append(float(lin))
                intensities.append(0.0)
            continue

        #Now, read in actual intensities and apply weights.  Keep track of maximum
        lin = line.split(',')
        weight = probs[lin[0]]
        for n_, li in enumerate(lin[1:]):
            if(n_ > len(twoThetas) - 1):
                break
            intensities[n_] += weight * float(li)

            if(twoThetas[n_] > maxTwoTheta):
                maxTwoThetaIndex = n_
                maxTwoTheta = twoThetas[n_]

    infile.close()

#Plot the thing
##(re)Normalize
maxInt = intensities[0]
for inten in intensities:
    if(inten > maxInt):
        maxInt = inten
x = twoThetas[:maxTwoThetaIndex]
y = [inten/maxInt for inten in intensities[:maxTwoThetaIndex]]

##Plot
fig, ax = plt.subplots(1, 1, figsize=(gf.WIDTH(1), gf.HEIGHT(1)))
ax.plot(x, y, color='k', linewidth=gf.LNSIZE_BIG)

ax.set_xlabel(r"$2\theta$ (degrees)", fontsize=gf.FONTSIZE_LABELS, fontweight="bold")
ax.set_ylabel("Intensity", fontsize=gf.FONTSIZE_LABELS, fontweight="bold")

ax.set_xlim(10, 50)
ax.set_ylim(0., 1.1)

ax = gf.SetTickSize(ax)

#Hardcoded miller indeces because automatically doing it isnt worth the trouble
deg12 = r"$(\!1\!1\!1\!)$"
deg23 = r"$(\!3\!1\!1\!)$"
deg28 = r"$(\!4\!0\!0\!)$"
deg40 = r"$(\!4\!4\!0\!)$"
ax.text(12, 0.53, deg12, fontsize=gf.FONTSIZE_TEXT, fontweight="normal")
ax.text(23 - 1, 0.31, deg23, fontsize=gf.FONTSIZE_TEXT, fontweight="normal")
ax.text(28, 1.025, deg28, fontsize=gf.FONTSIZE_TEXT, fontweight="normal")
ax.text(40, 0.78, deg40, fontsize=gf.FONTSIZE_TEXT, fontweight="normal")

fig.subplots_adjust(left=0.15, right=0.97, top=0.99, bottom=0.14, hspace=0., wspace=0.)
plt.savefig(SAVE_LOC + ".png", dpi=gf.DPI)
plt.savefig(SAVE_LOC + ".pdf")
