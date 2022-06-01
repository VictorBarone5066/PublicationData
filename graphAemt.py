from matplotlib import pyplot as plt
from headerAemt import *
from globalFig import *


nDat = InfileToEigDat("agbii4-cubic-cryst-202020-gga-meff//2612//out.aem", p=False)
pDat = InfileToEigDat("agbii4-cubic-cryst-202020-gga-meff//7497//out.aem", n=False)
pDatHi=InfileToEigDat("agbii4-cubic-cryst-202020-gga-meff//1778//out.aem", n=False)
SAVE_LOC = "meff"
WHICH = "cond" ##"cond" or "doss" for conductivity mass (arithmatic avg) or dos mass (geometric avg)

BAND_EDGES = dict()
BAND_EDGES = {'n': 3.596355, ##from 2612 - the lowest energy structure from the "high energy" group
              'p': 1.775} ##from 7497 - the lowest energy structure from the "low energy" group
"""
if(INFILE.split("//")[1] == "1778"):
    BAND_EDGES = {'n': 3.414767,
                  'p': 1.720072}
if(INFILE.split("//")[1] == "2612"):
    BAND_EDGES = {'n': 3.596355,
                  'p': 1.933096}
if(INFILE.split("//")[1] == "7497"):
    BAND_EDGES = {'n': 3.839319,
                  'p': 1.775}
"""

MIDGAP = (BAND_EDGES['n'] + BAND_EDGES['p'])/2
BANDGAP = 1.694695
BANDGAP = 1.8


COLORS = {1: COLORS_DICT["red"],
          2: COLORS_DICT["blue"],
          3: COLORS_DICT["black"],
          4: COLORS_DICT["purple"],
          5: COLORS_DICT["dGreen"]}
h = {0: "dotted",
     1: "dashdot",
     2: "dashed",
     3: "solid",
     4: (0, (8, 3)),
     5: (0, (1, 1))}

fig, ax = plt.subplots(2, 1, figsize=(WIDTH(1), HEIGHT(2)), sharex=True)



#(a, top): electrons
type_ = 'n'
n = 1
for n_, mu in enumerate(nDat.keys()):
    #if(mu < MIDGAP):
    #    continue
    x = nDat[mu][WHICH]["temp"]
    y = nDat[mu][WHICH]["mass"]

    lab = FloatFormat(abs(mu - BAND_EDGES['n']), 2)
    ax[0].plot(x, y, label=lab,
               color=COLORS[n],
               linewidth=LNSIZE_MED,
               linestyle=h[n])
    n += 1

#(b, bottom): holes
##Average the two low energy and high energy hole masses with specially weighted boltzmann avg
pdiffs = {"0.00": [], "0.01": [], "0.03": [], "0.10": [], "0.50": []}
x = []
for ke in pDat.keys():
    x = pDat[ke][WHICH]["temp"]
    kd = FloatFormat(abs(ke - 1.775), 2)
    pdiffs[kd] = [data*0.846 for data in pDat[ke][WHICH]["mass"]]
for ke in pDatHi.keys():
    kd = FloatFormat(abs(ke - 1.720072), 2)
    for i in range(0, len(pdiffs[kd])):
        pdiffs[kd][i] += 0.154*pDatHi[ke][WHICH]["mass"][i]

type_ = 'p'
n = 1
for n_, mu in enumerate(pdiffs.keys()):
    #if(mu > MIDGAP):
    #    continue

    y = pdiffs[mu]

    lab = FloatFormat((abs(float(mu))), 2) ##this is retarded
    ax[1].plot(x, y, label=lab,
               color=COLORS[n],
               linewidth=LNSIZE_MED,
               linestyle=h[n])
    n += 1


#General customization
for i in range(0, 2):
    ax[i].set_xlim(200, 600)
    ax[i].set_xticks([200, 300, 400, 500, 600])
    ax[i] = SetTickSize(ax[i])
    ax[i].grid(alpha=GRID_ALPHA)
##legend
ax[0].legend(fontsize=FONTSIZE_TEXT-2, fancybox=False, labelspacing=0.5, borderpad=0.5, loc=4,
                framealpha=1.0, title=r"CBM$-\mu\,$(eV)")
ax[1].legend(fontsize=FONTSIZE_TEXT-2, fancybox=False, labelspacing=0.5, borderpad=0.5, loc=4,
                framealpha=1.0, title=r"$\mu-$VBM$\,$(eV)")


ax[0] = RmAxTicks(ax[0], x=True, y=False)
ax[1].set_xlabel("Temperature (K)", fontsize=FONTSIZE_LABELS)
ax[0].set_ylabel(r"$\langle m^*_\mathrm{e}/m_0\rangle$", fontsize=FONTSIZE_LABELS)
ax[1].set_ylabel(r"$\langle m^*_\mathrm{h}/m_0\rangle$", fontsize=FONTSIZE_LABELS)
if(WHICH == "cond"):
    ax[0].set_ylim(0.95, 1.65)
    ax[0].set_yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
    ax[1].set_ylim(5, 17)
    ax[1].set_yticks([6, 8, 10, 12, 14, 16])
    pass
if(WHICH == "doss"):
    ax[0].set_ylim(1.2, 2.2)
    ax[0].set_yticks([1.3, 1.5, 1.7, 1.9, 2.1])
    ax[1].set_ylim(6, 18)
    ax[1].set_yticks([7, 9, 11, 13, 15, 17])

ax[0].text(x=0.04, y=0.935, s="a)", fontsize=FONTSIZE_LABELS, fontweight="bold", transform=ax[0].transAxes)
ax[1].text(x=0.04, y=0.935, s="b)", fontsize=FONTSIZE_LABELS, fontweight="bold", transform=ax[1].transAxes)

fig.subplots_adjust(bottom=0.07, top=0.99, left=0.16, right=0.96, hspace=0.045, wspace=0.0)
plt.savefig(SAVE_LOC + ".png", dpi=DPI)
plt.savefig(SAVE_LOC + ".pdf")
#plt.show()


##Print out a bunch of stuff
for mu in nDat.keys():
    plt.plot(nDat[mu][WHICH]["temp"], nDat[mu][WHICH]["mass"], label=str(mu))
    for ke in ['1', '2', '3']:
        for m, t, v in zip(nDat[mu][ke]["mass"], nDat[mu][ke]["temp"], nDat[mu][ke]["vec"]):
            print(mu, ke, "   ", f'{m:.3f}', '/', f'{t:.3f}',
                  '(', f'{v[0]:.3f}', f'{v[1]:.3f}', f'{v[2]:.3f}', ')')
    print("\n-------------------\n")
print("\n\n\n")

for mu in pDat.keys():
    #if(1.61 < mu < 1.63):
    #    continue
    plt.plot(pDat[mu][WHICH]["temp"], pDat[mu][WHICH]["mass"], label=str(mu))
    for ke in ['1', '2', '3']:
        for m, t, v in zip(pDat[mu][ke]["mass"], pDat[mu][ke]["temp"], pDat[mu][ke]["vec"]):
            print(mu, ke, "   ", f'{m:.3f}', '/', f'{t:.3f}',
                  '(', f'{v[0]:.3f}', f'{v[1]:.3f}', f'{v[2]:.3f}', ')')
    print("\n-------------------\n")
print("\n\n\n")

IND = 0
print("elec")
for mu in nDat.keys():
    for d in ["xx", "yy", "zz"]:
        print(mu, d, nDat[mu][d]["temp"][IND], nDat[mu][d]["mass"][IND])
print("\n\nhole")
for mu in pDat.keys():
    for d in ["xx", "yy", "zz"]:
        print(mu, d, pDat[mu][d]["temp"][IND], pDat[mu][d]["mass"][IND])
