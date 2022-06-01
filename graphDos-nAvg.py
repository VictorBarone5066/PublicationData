from matplotlib import pyplot as plt
import matplotlib.patches as mpat

import headerBoltzmann as hb
import headerDos as hd
import headerInterpolate as hi
import globalFig as gf

TEMPERATURE = 600*5000 #kelvin
BOLTZMANN_LOC = "agbii4-cubic-cryst-444-gga-boltzmann//output.csv"
DOSCAR_LOC = "agbii4-cubic-cryst-444-gga-doscars//"
CONTCAR_LOC = "agbii4-cubic-cryst-444-gga-contcars//"
SAVE_LOC = "dos"

INCLUDE_FULL = True

fig, ax = plt.subplots(1, 1, figsize=(gf.WIDTH(1), gf.HEIGHT(1)))

            #color, hatch, fill?, alpha, linestyle
col = {'Ag': [gf.RUD_AG, None, True, 0.5, "dotted"],
       'Bi': [gf.RUD_BI, None, True, 0.5, "dashed"],
       'I': [gf.RUD_I, None, False, 1.0, "dashdot"]
       }
if(INCLUDE_FULL):
    col={"Full": ["green",  None, False, 1.0, "solid"],
         'Ag':   [gf.RUD_AG, None, True, 0.5, "dotted"],
         'Bi':   [gf.RUD_BI, None, True, 0.5, "dashed"],
         'I':    [gf.RUD_I,  None, False, 1.0, "dashdot"]
         }

energy = None
elemScaledDos = dict()

outputDat = hb.GetOutputData(BOLTZMANN_LOC)
hb.SetAllProbs(datLis=outputDat, temp=TEMPERATURE)
for n, dat in enumerate(outputDat):
    print(n + 1, '/', len(outputDat))
    #if(n != 0):
    #    continue
    genNum = dat.wDirectory

    #Get this file's doscar data, initialize s+p+d, apply savgol
    doscarData = hd.ReadDoscar(doscarLoc=DOSCAR_LOC + "DOSCAR-" + str(genNum),
                                poscarLoc=CONTCAR_LOC + "CONTCAR-" + str(genNum),
                                iSpin=False, ncl=False, hasF=False)
    for d in doscarData:
        d.InitializeSums()
        d.ApplySavgol(winsize=35, poly=5)
        d.ScaleByEfermi()
    if(INCLUDE_FULL):
        full = doscarData[0] if (doscarData[0].atomType == "full") else None
    sums = hd.SumElementDos(doscarData)

    #Get a consistent grid to eval the DOS on, apply boltzmann scaling
    ##Full
    if(INCLUDE_FULL):
        energy, full.dosInfo["dos"] = hi.ReEval(x=full.energy, y=full.dosInfo["dos"], min=-5, max=5,
                                                step=0.01)
        if("Full" not in  elemScaledDos.keys()):
            elemScaledDos["Full"] = [d*dat.prob for d in full.dosInfo["dos"]]
        else:
            for i in range(0, len(elemScaledDos["Full"])):
                elemScaledDos["Full"][i] += full.dosInfo["dos"][i]*dat.prob
    ##Atom
    for ke, va in sums.items():
        energy, va.summedDosInfo["dosSum"] = hi.ReEval(x=va.energy, y=va.summedDosInfo["dosSum"],
                                                       min=-5, max=5, step=0.01)
        if(ke not in elemScaledDos.keys()):
            elemScaledDos[ke] = [v*dat.prob for v in va.summedDosInfo["dosSum"]]
        else:
            for i in range(0, len(elemScaledDos[ke])):
                elemScaledDos[ke][i] += va.summedDosInfo["dosSum"][i]*dat.prob


#Plot things
maxDos = 0.
for ke, va in elemScaledDos.items():
    for v in va:
        if (v > maxDos):
            maxDos = v
    #"Full" DOS
    ax.plot(energy, va, color=col[ke][0], linewidth=gf.LNSIZE_MED, linestyle=col[ke][4], label=ke)
    #ax[0].fill(energy, va, color=col[ke][0], fill=col[ke][2], hatch=col[ke][1], alpha=col[ke][3])



ax.grid(alpha=gf.GRID_ALPHA)

fig.supxlabel(r"$E-E_\mathrm{F}$ (eV)", fontsize=gf.FONTSIZE_LABELS, fontweight="bold")
ax.set_ylabel(r"DOS", fontsize=gf.FONTSIZE_LABELS, fontweight="bold")
ax = gf.RmAxTicks(ax, x=False, y=True)


ax.set_xlim(-4, 4)
ax.set_ylim(0, maxDos + 0.5)


ax = gf.SetTickSize(ax, x=True, y=False)


#ax.text(x=-3.8, y=(maxDos + 0.5)*93.5/100, s="a)", fontsize=gf.FONTSIZE_TEXT, fontweight="bold",
#           backgroundcolor="white")


#legend
#pats = []
#for k in col.keys():
#    if(False):#k == 'I'): ##save in case we want to do something special with iodine
#        pats.append(mpat.Patch(color=col[k][0], label=k, hatch='.', fill=col[k][2], alpha=col[k][3]))
#    else:
#        pats.append(mpat.Patch(color=col[k][0], label=k, hatch=col[k][1], fill=col[k][2], alpha=col[k][3]))
#fig.legend(handles=pats, loc="center", ncol=len(pats), bbox_to_anchor=(0.5, 0.97),
#           fontsize=gf.FONTSIZE_TEXT, frameon=False)
fig.legend(loc="center", ncol=len(col.keys()), bbox_to_anchor=(0.5, 0.97),
           fontsize=gf.FONTSIZE_TEXT - 1, frameon=False)

fig.subplots_adjust(bottom=0.13, top=0.93, left=0.1, right=0.99, hspace=0., wspace=0.075)
plt.savefig(SAVE_LOC + ".png", dpi=gf.DPI)
plt.savefig(SAVE_LOC + ".pdf")
