from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import linspace
import numpy as np
from matplotlib import cm
from globalFig import *

INFILE_LOC = "agbii4-cubic-cryst-202020-gga-meff//"
VAL, CON = "7497//", "1778//"
WHICHCON, WHICHVAL = 0, 4
SAVE_LOC = "bzInfo"
MOD = 150
CUT = True


AZIMUTH = 315
ELEVATION = 25

def DSqd(x1, x2, x3, y1, y2, y3):
    return (y1-x1)**2 + (y2-x2)**2 + (y3-x3)**2

#Read, find extrema for valance band
vMaxE, vMinE = -5000000, 5000000
vkx, vky, vkz, ve = [], [], [], []
with open(INFILE_LOC + VAL + "fullV", 'r') as infile:
    n = 0
    which = -1
    for lin in infile:
        if("BEG" in lin):
            which += 1
            vMaxE, vMinE = -5000000, 5000000
            continue
        if(which != WHICHVAL):
            continue
        if("END" in lin):
            if(which == WHICHVAL):
                break
            continue
        if(n%MOD != 0):
            n += 1
            #continue
        line = lin.split()

        nrg = float(line[3])
        if(nrg < vMinE):
            vMinE = nrg
        if(vMaxE < nrg):
            vMaxE = nrg

        vkx.append(float(line[0]))
        vky.append(float(line[1]))
        vkz.append(float(line[2]))
        ve.append(nrg)

        n += 1

    infile.close()

#Read, find extrema for conduction band
cMaxE, cMinE = -5000000, 5000000
ckx, cky, ckz, ce = [], [], [], []
with open(INFILE_LOC + CON + "fullC", 'r') as infile:
    n = 0
    which = -1
    for lin in infile:
        if("BEG" in lin):
            which += 1
            cMaxE, cMinE = -5000000, 5000000
            continue
        if(which != WHICHCON):
            continue
        if("END" in lin):
            if(which == WHICHCON):
                break
            continue
        if(n%MOD != 0):
            n += 1
            #continue
        line = lin.split()

        nrg = float(line[3])
        if(nrg < cMinE):
            cMinE = nrg
        if(cMaxE < nrg):
            cMaxE = nrg

        ckx.append(float(line[0]))
        cky.append(float(line[1]))
        ckz.append(float(line[2]))
        ce.append(nrg)

        n += 1

    infile.close()

print(vMinE, vMaxE, len(ve))
print(cMinE, cMaxE, len(ce))

a = 0.24600000
exX = [-a, -a, -a, +a, +a, +a, 0., +a, +a, 0., 0., 0., +a]
exY = [-a, -a, +a, +a, +a, -a, -a, -a, 0., 0., 0., -a, 0.]
exZ = [-a, +a, +a, +a, -a, -a, 0., 0., 0., +a, 0., +a, +a]

#Remove some unwanted stuff for better visualization
#Also re-scale energies to fractions between maxima and minima (for multi-colorbar reasons)
##Valance
tx, ty, tz, te = [], [], [], []
for i in range(0, len(ve)):
    if(vkx[i] > 0 and vky[i] < 0 and vkz[i] > 0):
        continue
    for k in range(0, len(exX)):
        if (DSqd(vkx[i], vky[i], vkz[i], exX[k], exY[k], exZ[k]) < 0.0001):
            vkx[i] = -0.01
            vky[i] = -0.01
            vkz[i] = -0.01
    tx.append(vkx[i])
    ty.append(vky[i])
    tz.append(vkz[i])
    te.append((ve[i] - vMinE) / (vMaxE - vMinE))
vkx, vky, vkz, ve = tx, ty, tz, te
##Conduction
tx, ty, tz, te = [], [], [], []
for i in range(0, len(ce)):
    if(ckx[i] > 0 and cky[i] < 0 and ckz[i] > 0):
        continue
    for k in range(0, len(exX)):
        if (DSqd(ckx[i], cky[i], ckz[i], exX[k], exY[k], exZ[k]) < 0.0001):
            ckx[i] = -0.01
            cky[i] = -0.01
            ckz[i] = -0.01
    tx.append(ckx[i])
    ty.append(cky[i])
    tz.append(ckz[i])
    te.append((ce[i] - cMinE) / (cMaxE - cMinE))
ckx, cky, ckz, ce = tx, ty, tz, te


#Plot the two band edges
fig = plt.figure(figsize=(3.5, 3.5))
Z = 0.24
##Valance
ax = fig.add_subplot(2, 2, (1), projection="3d")
inicmap = cm.get_cmap('jet')
revcmap = inicmap.reversed()
ax.scatter(exX, exY, exZ, c='k', s=5.5, alpha=1, depthshade=False)
img = ax.scatter(vkx, vky, vkz, c=ve, cmap=inicmap, s=1.5, zorder=-100)
ax.text(0.0075, 0.0075, 0.0075, r"$\Gamma$", fontsize=14)
ax.text(a + 0.0075, a + 0.0075, a + 0.0075, r"R", fontsize=14)
ax.text(0.07, -a - 0.2, -a, r"$\vec{k_\mathrm{x}\!\!}$", fontsize=14)
ax.text(a + 0.15, 0.00 - a / 3, -a, r"$\vec{k_\mathrm{y}\!\!\!}$", fontsize=14)
ax.text(0.07, 2 * a, -a, r"$\vec{k_\mathrm{z}\!\!}$", fontsize=14)
ax.axis('off')
ax.set_xlim(-Z, Z)
ax.set_ylim(-Z, Z)
ax.set_zlim(-Z, Z)
ax.view_init(azim=AZIMUTH, elev=ELEVATION)
##Conduction
ax = fig.add_subplot(2, 2, (2), projection="3d")
inicmap = cm.get_cmap('jet')
revcmap = inicmap.reversed()
ax.scatter(exX, exY, exZ, c='k', s=5.5, alpha=1, depthshade=False)
img = ax.scatter(ckx, cky, ckz, c=ce, cmap=revcmap, s=1.5, zorder=-100)
ax.text(0.0075, 0.0075, 0.0075, r"$\Gamma$", fontsize=14)
ax.text(a + 0.0075, a + 0.0075, a + 0.0075, r"R", fontsize=14)
ax.text(0.07, -a - 0.2, -a, r"$\vec{k_\mathrm{x}\!\!}$", fontsize=14)
ax.text(a + 0.15, 0.00 - a / 3, -a, r"$\vec{k_\mathrm{y}\!\!\!}$", fontsize=14)
ax.text(0.07, 2 * a, -a, r"$\vec{k_\mathrm{z}\!\!}$", fontsize=14)
ax.axis('off')
ax.set_xlim(-Z, Z)
ax.set_ylim(-Z, Z)
ax.set_zlim(-Z, Z)
ax.view_init(azim=AZIMUTH, elev=ELEVATION)



#Plot band structures
from scipy.signal import savgol_filter as sf
##Valance
###Read
FILES = [INFILE_LOC + VAL + "vx", INFILE_LOC + VAL + "vz"]
EXTENTS = [[0, 0.5], [0.5, 1]]
HOWMANY = 5
disps = dict()
maxE = -500000
kExts = dict()
for i in range(0, HOWMANY):
    disps[i] = {"x": [[], []], "z": [[], []]}
    kExts[i] = {"x": [500, -500], "z": [500, -500]}
for infileN in FILES:
    key = infileN[len(INFILE_LOC + VAL) + 1:]
    n = -1
    with open(infileN, 'r') as infile:
        for lin in infile:
            if ("BEG" in lin):
                n += 1
                continue
            if ("END" in lin):
                continue

            line = lin.split()
            nrg = float(line[1])
            k = float(line[0])

            if(k < kExts[n][key][0]): ##minimum k
                kExts[n][key][0] = k
            if(kExts[n][key][1] < k): ##maximum k
                kExts[n][key][1] = k
            if(maxE < nrg):
                maxE = nrg

            disps[n][key][0].append(k)
            disps[n][key][1].append(nrg)
###Plot
ax = fig.add_subplot(2, 2, (3))
for n in range(0, HOWMANY):
    for keI, ke in enumerate(disps[n].keys()):
        k = [((k_ - kExts[n][ke][0]) / (kExts[n][ke][1] - kExts[n][ke][0]) - 1) * \
             (EXTENTS[keI][1] - EXTENTS[keI][0]) + EXTENTS[keI][1] for k_ in disps[n][ke][0]]
        e = [1000*(e_ - maxE) for e_ in disps[n][ke][1]]

        k, e = (list(t) for t in zip(*sorted(zip(k, e))))
        e = sf(e, 7, 1)
        ax.plot(k, e, 'b')
for ext in EXTENTS:
    ax.axvline(x=ext[0], color='k')
    ax.axvline(x=ext[1], color='k')
##Ticks
ax.set_xlim(0, 1)
ax.set_yticks([0, -20, -40, -60, -80, -100])
ax.set_ylim(-0.1*1000, 0.005*1000)
ax.set_aspect(0.0115)
ax = RmAxTicks(ax=ax, y=False)
##Labels
ax.tick_params(axis='y', labelrotation = 0)
ax.set_ylabel(r"$E-$Band Edge (meV)")
ax.text(EXTENTS[0][0] + (EXTENTS[0][1] - EXTENTS[0][0])/2, -0.1125*1000, "[100]", ha="center", fontsize=12)
ax.text(EXTENTS[1][0] + (EXTENTS[1][1] - EXTENTS[1][0])/2, -0.1125*1000, "[001]", ha="center", fontsize=12)

##Conduction
###Read
FILES = [INFILE_LOC + CON + "cx", INFILE_LOC + CON + "cz"]
EXTENTS = [[0, 0.5], [0.5, 1]]
HOWMANY = 5
minE = 5000000
disps = dict()
kExts = dict()
for i in range(0, HOWMANY):
    disps[i] = {"x": [[], []], "z": [[], []]}
    kExts[i] = {"x": [500, -500], "z": [500, -500]}
for infileN in FILES:
    key = infileN[len(INFILE_LOC + VAL) + 1:]
    n = -1
    with open(infileN, 'r') as infile:
        for lin in infile:
            if ("BEG" in lin):
                n += 1
                continue
            if ("END" in lin):
                continue

            line = lin.split()
            nrg = float(line[1])
            k = float(line[0])

            if(k < kExts[n][key][0]): ##minimum k
                kExts[n][key][0] = k
            if(kExts[n][key][1] < k): ##maximum k
                kExts[n][key][1] = k
            if(nrg < minE):
                minE = nrg

            disps[n][key][0].append(k)
            disps[n][key][1].append(nrg)
###Plot
ax = fig.add_subplot(2, 2, (4))
for n in range(0, HOWMANY):
    for keI, ke in enumerate(disps[n].keys()):
        k = [((k_ - kExts[n][ke][0]) / (kExts[n][ke][1] - kExts[n][ke][0]) - 1) * \
             (EXTENTS[keI][1] - EXTENTS[keI][0]) + EXTENTS[keI][1] for k_ in disps[n][ke][0]]
        e = [1000*(e_ - minE) for e_ in disps[n][ke][1]]

        k, e = (list(t) for t in zip(*sorted(zip(k, e))))
        e = sf(e, 7, 1)
        ax.plot(k, e, 'b')
for ext in EXTENTS:
    ax.axvline(x=ext[0], color='k')
    ax.axvline(x=ext[1], color='k')
##Ticks
ax.set_xlim(0, 1)
ax.set_ylim(-0.005*1000, 0.500*1000)
ax.set_yticks([0, 125, 250, 375, 500])
ax.set_aspect(0.00240)
ax = RmAxTicks(ax=ax, y=False)
##Labels
#ax.set_ylabel(r"$E-$CBM (eV)")
ax.text(EXTENTS[0][0] + (EXTENTS[0][1] - EXTENTS[0][0])/2, 1000*-0.065, "[100]", ha="center", fontsize=12)
ax.text(EXTENTS[1][0] + (EXTENTS[1][1] - EXTENTS[1][0])/2, 1000*-0.065, "[001]", ha="center", fontsize=12)


#Plot color bar
inicmap = cm.get_cmap('jet')
revcmap = inicmap.reversed()
ax = fig.add_subplot(1, 2, (1, 2))
im = ax.imshow(np.array([[0, 0], [1, 1]]), vmin=0, vmax=1, aspect=0.05, cmap=revcmap)

ax.text(0.01, 1.25, "VBM", fontsize=12)
ax.text(0.78, 1.25, FloatFormat(vMinE - vMaxE, 2) +  " eV", fontsize=12)
ax.text(0.5, 1.25, "V.B. Energy", fontsize=12, ha="center")
ax.text(0.01, -1.05, "CBM", fontsize=12)
ax.text(0.80, -1.05, FloatFormat(cMaxE - cMinE, 2) + " eV", fontsize=12)
ax.text(0.5, -1.05, "C.B. Energy", fontsize=12, ha="center")

colorbar = fig.colorbar(im, cax=ax, orientation="horizontal")
#colorbar.set_label("C.B. Energy", fontsize=12)
colorbar.set_ticks([])

#Final Updates
fig.text(0.15, 0.95, "a)", fontsize=12)
fig.text(0.57, 0.95, "b)", fontsize=12)
fig.text(0.175, 0.415, "c)", fontsize=12)
fig.text(0.611, 0.415, "d)", fontsize=12)

#plt.subplot_tool()
plt.subplots_adjust(hspace=-0.05, wspace=0.348, left=0.17, right=0.93, top=1.15, bottom=-0.05)
#plt.show()
plt.savefig(SAVE_LOC + ".pdf")
plt.savefig(SAVE_LOC + ".png", dpi=350)
