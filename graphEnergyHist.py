from matplotlib import pyplot as plt
import matplotlib.patches as mpat

import headerReadOutput as hro
import globalFig as gf

BOLTZMANN_LOC = "agbii4-cubic-cryst-111-gga-bondInfo//output.csv"
SAVE_LOC = "histogram"

BIN_WIDTH = 0.25 #eV

fig, ax = plt.subplots(1, 1, figsize=(gf.WIDTH(1), gf.HEIGHT(1)))
allDat = hro.GetOutputData(BOLTZMANN_LOC)

scaledEnergies = sorted([d.nrg for d in allDat])
scaledEnergies = [e - scaledEnergies[0] for e in scaledEnergies]

width, bins = 0., []
while(width < scaledEnergies[-1]):
    bins.append(width)
    width += BIN_WIDTH
bins.append(width)
res = ax.hist(scaledEnergies, facecolor=gf.RUD_BI, edgecolor=gf.RUD_I, bins=bins)

for i in range(0, len(bins) - 1):
    yOffset = -125 + -100*(len(str(int(res[0][i]))) - 2)
    if(res[0][i] < 300):
        yOffset = 80
    ax.text(res[1][i] + BIN_WIDTH/4., res[0][i]+yOffset, str(int(res[0][i])), fontsize=gf.FONTSIZE_TEXT,
            rotation=90)

ax.set_xticks([bins[b] for b in range(0, len(bins)) if b%2 == 0])
ax = gf.SetTickSize(ax)

ax.set_xlabel(r"$E-E_\mathrm{min}$ (eV)", fontsize=gf.FONTSIZE_LABELS)
ax.set_ylabel("Count", fontsize=gf.FONTSIZE_LABELS)


fig.subplots_adjust(bottom=0.14, top=0.99, left=0.185, right=0.98, hspace=0., wspace=0.)
plt.savefig(SAVE_LOC + ".png", dpi=gf.DPI)
plt.savefig(SAVE_LOC + ".pdf")
