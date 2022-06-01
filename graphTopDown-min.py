from matplotlib import pyplot as plt

import headerPoscar as hp
import headerTopDown as htd
import headerBoltzmann as hb
import globalFig as gf

#figure data
            #elem:  color, radius
POINT_MAP = {"Ag": [gf.RUD_AG, 0.05],
             "Bi": [gf.RUD_BI, 0.05],
             "I":  [gf.RUD_I, 0.03]}

#file i/o data
IDEAL_POSCAR_LOC = "agbii4-cubic-cryst-idealPositions"
POSCARS_LOC = "agbii4-cubic-cryst-111-gga-poscars//"
CONTCARS_LOC = "agbii4-cubic-cryst-444-gga-contcars//"
MIN_MODEL_NUM = 7497

SAVE_LOC = "topDownMin"

def IndexToRowCol(index):
    if(index == 0):
        return 0, 0
    if(index == 1):
        return 0, 1
    if(index == 2):
        return 1, 0
    if(index == 3):
        return 1, 1

def SiteEqualAtom(site, atom, tol=0.1):
    return (abs(site.a - atom.a) < tol and abs(site.b - atom.b) < tol and abs(site.c - atom.c) < tol)

class site:
    a, b, c = None, None, None
    elemInfo = None ##a dict of {genNum: [element, boltzmann probability]} for this site

    def __init__(self, a_, b_, c_):
        self.a = None
        self.b = None
        self.c = None
        self.elemInfo = dict()

        self.a = a_
        self.b = b_
        self.c = c_

        return

    def AddToElemInfo(self, genNum, elem, boltzmann):
        self.elemInfo[genNum] = [elem, boltzmann]
        return

    #Returns a dict of {elem: probability}
    def GiveFractionalOccupancies(self):
        toRet = dict()
        for va in self.elemInfo.values():
            if(va[0] not in toRet.keys()):
                toRet[va[0]] = va[1]
            else:
                toRet[va[0]] += va[1]
        return dict(sorted(toRet.items()))

#Base poscar: need to connect all sites to this POSCAR's sites
base = hp.Poscar(IDEAL_POSCAR_LOC)
base.MoveAtomsToUnitCell()
base.GetAtomEquivPositions()

##Fill base sites
sites = []
for atom in base.atoms:
    for equiv in atom.equivPositions:
        if(-0.1 < equiv.a < 1.1 and -0.1 < equiv.b < 1.1 and -0.1 < equiv.c < 1.1):
            sites.append(site(equiv.a, equiv.b, equiv.c))


##Get the CONTCAR corresponding to this file
thisContcar = hp.Poscar(POSCARS_LOC + "POSCAR0-" + str(MIN_MODEL_NUM))
thisContcar.MoveAtomsToUnitCell()
thisContcar.GetAtomEquivPositions()

nPlanes, fracs = htd.GetNumLayers(poscar=base, along='c', returnFractions=True)
fig, ax = plt.subplots(int(nPlanes/2), int(nPlanes/2), figsize=(gf.WIDTH(1), gf.HEIGHT(1)))
for index in range(0, nPlanes):
    row, col = IndexToRowCol(index)
    ax[row][col].text(0.5, 1.05, r"$c$ = " + gf.FloatFormat(fracs[index], dec=2), fontsize=gf.FONTSIZE_LABELS,
                      horizontalalignment='center')
    ax[row][col] = htd.GetPlane(ax=ax[row][col], referencePoscar=base, nGridLines=3, along='c')

    atoms = htd.GetAtomsInPlane(poscar=thisContcar, along='c', layer=index, tol=0.05)

    ##And then fill all equivalent positions with element and boltzmann data
    for atom in atoms:
        for equiv in atom.equivPositions:
            if (-0.1 < equiv.a < 1.1 and -0.1 < equiv.b < 1.1 and -0.1 < equiv.c < 1.1):
                ##Here, we've id'd a site that belongs to this plane.
                thisX, thisY = htd.DirectToPlot(xDirect=equiv.a, yDirect=equiv.b, referencePoscar=base,
                                                along='c')

                ##Prep points for plotting
                colors, fOccs, radaii = [], [], []
                fOccs.append(1.0)
                colors.append(POINT_MAP[equiv.atomType][0])
                radaii.append(POINT_MAP[equiv.atomType][1])

                ##Plot the points
                ax[row][col] = htd.MarkFractionalPoint(ax=ax[row][col], where=(thisX, thisY),
                                                       rad=radaii, fractions=fOccs, colors=colors)

#Final plot adjustments
fig.text(0.001, 0.95, "b)", fontsize=gf.FONTSIZE_LABELS)
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.97, wspace=0.01, hspace=0.03)


plt.savefig(SAVE_LOC + ".png", dpi=gf.DPI)
plt.savefig(SAVE_LOC + ".pdf")
