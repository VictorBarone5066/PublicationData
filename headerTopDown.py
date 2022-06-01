#For creating top-down views of crystals along an axis
import headerPoscar as hp
import globalFig as gf
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge as wedge

###-----------------------------------###
# Functions for determining planes, etc #
###-----------------------------------###

#Gives the number of layers needed to capture all atoms in seperate planes for a poscar
#tol = in fractional
#will fail if your material is not easily split into layers along the chosen plane, probably
#If returnFractions, returns the fractional portions used to differentiate the planes
def GetNumLayers(poscar, along='c', tol=0.05, returnFractions=False):
    poscar.ConvertToDirect()

    fracs = []
    if(along == 'a' or along == 'x'):
        fracs = [atom.a for atom in poscar.atoms]
    elif(along == 'b' or along == 'y'):
        fracs = [atom.b for atom in poscar.atoms]
    elif(along == 'c' or along == 'z'):
        fracs = [atom.c for atom in poscar.atoms]

    ##move stuff bc of periodic bcs
    for n, f in enumerate(fracs):
        if (1.0 - f < tol):  ##1-sf is always >=0 unless you mess up
            fracs[n] = 1.0 - f ##f = 1 - f does NOT change fracs, for some insane reason.  Need this instead
    sortedFracs = sorted(fracs)

    #Count planes by differentiating between 'large gaps' in the sorted fractions.  Each time we cross a
    #gap, we're on a new plane
    planeCount = 0
    currentFrac = -2.*tol
    fracsToReturn = []
    for sf in sortedFracs:
        if(sf - currentFrac > tol):
            planeCount += 1
            currentFrac = sf
            fracsToReturn.append(sf)

    if(returnFractions):
        return planeCount, fracsToReturn
    return planeCount

#Get a list of hp.Atom instances for a given plane on a given direction.  layer number is zero-indexed
def GetAtomsInPlane(poscar, along='c', layer=0, tol=0.05):
    poscar.ConvertToDirect()
    poscar.MoveAtomsToUnitCell()

    atoms = [hp.deepcopy(atom) for atom in poscar.atoms]
    #For the a direction
    if (along == 'a' or along == 'x'):
        for atom in atoms:
            if (1.0 - atom.a < tol): ##The layer at r=0 is the same as r=1.  Move all to r=0 for my sanity
                atom.a = 1.0 - atom.a
        atoms.sort(key=lambda x: x.a)
        sortedFracs = [atom.a for atom in atoms]
    #For the b direction
    elif (along == 'b' or along == 'y'):
        for atom in atoms:
            if (1.0 - atom.b < tol): ##The layer at r=0 is the same as r=1.  Move all to r=0 for my sanity
                atom.b = 1.0 - atom.b
        atoms.sort(key=lambda x: x.b)
        sortedFracs = [atom.b for atom in atoms]
    #For the c direction
    elif (along == 'c' or along == 'z'):
        for atom in atoms:
            if (1.0 - atom.c < tol): ##The layer at r=0 is the same as r=1.  Move all to r=0 for my sanity
                atom.c = 1.0 - atom.c
        atoms.sort(key=lambda x: x.c)
        sortedFracs = [atom.c for atom in atoms]

    #Count planes by differentiating between 'large gaps' in the sorted fractions.  Each time we cross a
    #gap, we're on a new plane.
    planeCount = 0
    currentFrac = -2.*tol
    toRet = [atoms[0]]
    for sf, atom in zip(sortedFracs, atoms):
        toRet.append(atom)
        if(sf - currentFrac > tol):
            if(planeCount == layer + 1):
                return toRet[:-1] #we've overcounted by 1 -- don't return the last atom
            planeCount += 1
            currentFrac = sf
            toRet = [toRet[-1]]

    return toRet ##Here, we've run out of atoms, so the remining ones must be what were asked for



###-------------------------------###
# Functions for plotting stuff, etc #
###-------------------------------###

#Returns an axis with a 2D representation of a cut on a plane for a given poscar
#Optionally sets nGridLines evenly spaced gridlines (if nGridLines > 0)
def GetPlane(ax, referencePoscar, along='c', axColor=gf.COLORS_DICT["black"], axLineWidth=gf.LNSIZE_MED,
             nGridLines=0, gridLineColor=gf.COLORS_DICT["black"], gridLineWidth=gf.LNSIZE_SML,
             gridlineAlpha=gf.GRID_ALPHA, offset=0.1):
    # Draw my own vectors
    ##Normalize
    normA, normB, normC = 0., 0., 0.
    for i in range(0, 3):
        normA += referencePoscar.superCellVecA[i] ** (2)
        normB += referencePoscar.superCellVecB[i] ** (2)
        normC += referencePoscar.superCellVecC[i] ** (2)
    for i in range(0, 3):
        referencePoscar.superCellVecA[i] /= normA ** (1. / 2.)
        referencePoscar.superCellVecB[i] /= normB ** (1. / 2.)
        referencePoscar.superCellVecC[i] /= normC ** (1. / 2.)

    ##Decide which ones are important
    x0, y0 = 0.0, 0.0
    if(along == 'a' or along == 'x'):
        xf1 = referencePoscar.superCellVecB[1]
        xf2 = referencePoscar.superCellVecC[1]
        yf1 = referencePoscar.superCellVecB[2]
        yf2 = referencePoscar.superCellVecC[2]
    elif(along == 'b' or along == 'y'):
        xf1 = referencePoscar.superCellVecA[0]
        xf2 = referencePoscar.superCellVecC[0]
        yf1 = referencePoscar.superCellVecA[2]
        yf2 = referencePoscar.superCellVecC[2]
    elif(along == 'c' or along == 'z'):
        xf1 = referencePoscar.superCellVecA[0]
        xf2 = referencePoscar.superCellVecB[0]
        yf1 = referencePoscar.superCellVecA[1]
        yf2 = referencePoscar.superCellVecB[1]

    ##draw the vectors on matplotlib
    # ax.plot([x0, x1], [y0, y1])
    ax.plot([x0, xf1], [y0, yf1], color=axColor, linewidth=axLineWidth)  ##1
    ax.plot([x0, xf2], [y0, yf2], color=axColor, linewidth=axLineWidth)  ##2
    ax.plot([xf2, xf2 + xf1], [yf2, yf2 + yf1], color=axColor, linewidth=axLineWidth)  ##1's mirror
    ax.plot([xf1, xf1 + xf2], [yf1, yf1 + yf2], color=axColor, linewidth=axLineWidth)  ##2's mirror

    #Add gridlines if necessary
    if(nGridLines > 0):
        spacing = 1./float(nGridLines + 1)
        current = spacing
        while (current < 1.0):
            ax.plot([xf2*current, xf2*current + xf1], [yf2*current, yf2*current + yf1], ##1 cmpnt
                color=gridLineColor, alpha=gridlineAlpha, linewidth=gridLineWidth)
            ax.plot([xf1*current, xf1*current + xf2], [yf1*current, yf1*current + yf2], ##2 cmpnt
                color=gridLineColor, alpha=gridlineAlpha, linewidth=gridLineWidth)
            current += spacing

    #Make the box square, why is this so diffucult?
    ax.axis("off")
    maxLen = max(xf1 + xf2, yf1 + yf2)
    ax.set_aspect(1.0)
    ax.set_xlim(-offset, maxLen + offset)
    ax.set_ylim(-offset, maxLen + offset)

    return ax

#Marks a point on the plane with a pie-chart style circle.
#where = a tuple of (x, y) where to mark the point
#fractions = a list of how much of a circle gets filled with a certian color
# (make sure all fractions add to 1)
def MarkFractionalPoint(ax, where, rad, fractions, colors):
    currTheta = 0.
    for i in range(0, len(fractions)):
        nextTheta = currTheta + 360.*fractions[i]
        ax.add_artist(wedge(center=where, r=rad[i], theta1=currTheta, theta2=nextTheta, color=colors[i]))
        currTheta = nextTheta

    return ax

#Transforms direct coords to plotting coords
#Returns xPlot, yPlot
def DirectToPlot(xDirect, yDirect, referencePoscar, along='c'):
    if(along == 'a' or along == 'x'):
        x1 = referencePoscar.superCellVecB[1]
        x2 = referencePoscar.superCellVecC[1]
        y1 = referencePoscar.superCellVecB[1]
        y2 = referencePoscar.superCellVecC[2]
    if(along == 'b' or along == 'y'):
        x1 = referencePoscar.superCellVecA[0]
        x2 = referencePoscar.superCellVecC[0]
        y1 = referencePoscar.superCellVecA[2]
        y2 = referencePoscar.superCellVecC[2]
    if(along == 'c' or along == 'z'):
        x1 = referencePoscar.superCellVecA[0]
        x2 = referencePoscar.superCellVecB[0]
        y1 = referencePoscar.superCellVecA[1]
        y2 = referencePoscar.superCellVecB[1]

    return xDirect*x1 + yDirect*x2, xDirect*y1 + yDirect*y2



    #Plot the atoms
    atoms = GetAtomsInPlane(poscar=ideal, along='c', layer=col, tol=0.05)
    for atom in atoms:
        if(atom.atomType == 'Anion'):
            color='blue'
        else:
            color='k'
        for equiv in atom.equivPositions:
            if(-0.01 < equiv.a < 1.01 and -0.01 < equiv.b < 1.01):
                plotPoint = [equiv.a, equiv.b]
                plotPoint = DtoP(plotPoint)
                ax[col].plot(plotPoint[0], plotPoint[1], marker='o', color=color)
