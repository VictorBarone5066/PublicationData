from matplotlib import pyplot as plt
from numpy import polyfit

import headerBonds as hb
import globalFig as fg

#IO globals
INFILE_FOLDER = "agbii4-cubic-cryst-111-gga-bondInfo//"
OUT_DAT_LOC = INFILE_FOLDER + "output.csv"
INIT_PAIR_LOC = INFILE_FOLDER + "pairInfoIni.csv"
FINAL_PAIR_LOC = INFILE_FOLDER + "pairInfoFin.csv"

SAVE_LOC = "nrgvsPairs"

#Plotting globals
IMPORTANT_PAIRS = ["AgAg"]
TEXT_BOX_PARAMS = dict(facecolor='white', edgecolor='white', alpha=0.95, pad=0.15)

def PrettyPair(pair):
    ret = pair[0]
    for i in range(1, len(pair)):
        if(pair[i].isupper()):
            ret += '-'
        ret += pair[i]
    return ret

def GetSlope(x, y):
    return polyfit(x, y, 1)[0]

#Get all important data
dat = hb.GetPairedInfo(outLoc=OUT_DAT_LOC, initPairLoc=INIT_PAIR_LOC, finalPairLoc=FINAL_PAIR_LOC)
avgIni = hb.GetAvgDicts([d.pairsIni for d in dat])
avgFin = hb.GetAvgDicts([d.pairsFin for d in dat])
allElemPairs = []
for d in dat:
    for ke in list(d.pairsIni.keys()) + list(d.pairsFin.keys()):
        allElemPairs = allElemPairs + [ke] if ke not in allElemPairs else allElemPairs
eMin = min([d.nrg for d in dat])

#Plot
fig, ax = plt.subplots(1, 1, sharex=False, sharey=True, figsize=(fg.WIDTH(1), fg.HEIGHT(1)))

for col, pair in enumerate(IMPORTANT_PAIRS):
    x, y = [], []
    counts = dict()
    for d in dat:
        if(pair in d.pairsIni.keys()):
            x.append(d.pairsIni[pair]["pairCount"])
            y.append(d.nrg - eMin)
            if(x[-1] in counts.keys()):
                counts[x[-1]] += 1
            else:
                counts[x[-1]] = 0
    ax.plot(x, y, 'o', color=fg.COLORS_DICT["black"], markersize=fg.MKSIZE_MED)
    xLine, yLine = fg.PolyInterp(x=x, y=y, deg=1)
    ax.plot(xLine, yLine, color=fg.RUD_BI, linewidth=fg.LNSIZE_MED)

    y_, x_ = (list(t) for t in zip(*sorted(zip(y, x))))
    for elem in range(0, len(x_)):
        print(elem, x_[elem], y_[elem])

    #Individual customization
    ax.grid(True, alpha=fg.GRID_ALPHA)
    ax.set_xlim(7, 17)

    ##ticks
    if(col != 0):
        ax = fg.RmAxTicks(ax=ax, x=False, y=True)
    ax = fg.SetTickSize(ax=ax, fontsize=fg.FONTSIZE_TICKS)

    ##add average nn dist and slope to plot
    slopeTxt = r"Slope: " + fg.FloatFormat(GetSlope(x, y), dec=2) + " eV"
    ax.text(0.05, 0.90, slopeTxt, fontsize=fg.FONTSIZE_TEXT, transform=ax.transAxes, ha='left',
                 bbox=TEXT_BOX_PARAMS)



#General customization
##labels
ax.set_xlabel("Pair Count", fontsize=fg.FONTSIZE_LABELS)
ax.set_ylabel(r"$E-E_{\mathrm{min}}$ (eV)", fontsize=fg.FONTSIZE_LABELS)

##padding
plt.subplots_adjust(left=0.15, right=0.99, bottom=0.13, top=0.99, hspace=0.0, wspace=0.0)

#Save, end
plt.savefig(SAVE_LOC + ".pdf")
plt.savefig(SAVE_LOC + ".png", dpi=fg.DPI)
plt.clf()
