#Calculates single band effective masses for AgBiI4.  Too specific to generalize - don't even try
from numpy import polyfit as pf
from numpy import linspace as ls
from matplotlib import pyplot as plt

DIRS = ['x', 'y', 'z']
TYPES = ['c', 'v']
TYPES = ['v']
PLOT = True

fig, ax = plt.subplots(nrows=2, ncols=3)

for row, TYPE in enumerate(TYPES):
    for col, DIR in enumerate(DIRS):
        INFILE_LOC = "agbii4-cubic-cryst-202020-gga-meff//7497//" + TYPE + DIR + '_'

        PTS_AROUND = 4 ##number of points to the left or right of maxima to use for fit
        EXL, EXR = 0, 0 ##number of additional points around extrema to include
        #Very specific choices of extra points to include that provide reasonable poly fits
        if(DIR == 'x' and TYPE == 'c'):
            EXL, EXR = 2, 3
        if(DIR == 'x' and TYPE == 'v'):
            EXL, EXR = 1, 1
        if(DIR == 'y' and TYPE == 'c'):
            EXL, EXR = 2, 3
        if(DIR == 'y' and TYPE == 'v'):
            EXL, EXR = 1, 0
        if(DIR == 'z' and TYPE == 'c'):
            EXL, EXR = 3, 3
        if(DIR == 'z' and TYPE == 'v'):
            EXL, EXR = 1, 1


        def FindMaxIndex(x, y):
            ind, yHi = 0, y[0]
            for i in range(1, len(x)):
                tstY = y[i]
                if(tstY > yHi):
                    yHi = tstY
                    ind = i
            return ind

        def FindMinIndex(x, y):
            ind, yLo = 0, y[0]
            for i in range(1, len(x)):
                tstY = y[i]
                if(tstY < yLo):
                    yLo = tstY
                    ind = i
            return ind

        ##Read file, get all data
        allXPts, allYPts = [], []
        with open(INFILE_LOC, 'r') as infile:
            for lin in infile:
                line = lin.split()

                allXPts.append(float(line[0]))
                allYPts.append(float(line[2]))

            infile.close()
        allXPts, allYPts = (list(t) for t in zip(*sorted(zip(allXPts, allYPts))))

        #Get points to the left and right of the extrema
        pivot = FindMaxIndex(allXPts, allYPts) if (TYPE == 'v') else FindMinIndex(allXPts, allYPts)
        extXPts, extYPts = [], []
        ##left
        i = 0
        while(i != None):
            if(pivot - i < 0 or pivot - i > len(allXPts) - 1 or i > PTS_AROUND + EXL):
                break

            extXPts.append(allXPts[pivot - i])
            extYPts.append(allYPts[pivot - i])

            i += 1
        ##right
        i = 1
        while(i != None):
            if(pivot + i < 0 or pivot + i > len(allXPts) - 1 or i > PTS_AROUND + EXR):
                break

            extXPts.append(allXPts[pivot + i])
            extYPts.append(allYPts[pivot + i])

            i += 1
        i = None
        extXPts, extYPts = (list(t) for t in zip(*sorted(zip(extXPts, extYPts))))

        #Fit polynomial, print effective masses
        fit = pf(x=extXPts, y=extYPts, deg=2, full=False)
        a, b, c = fit[0], fit[1], fit[2]
        if(TYPE == 'v'):
            a = a*-1
        print(TYPE + DIR, 3.81/a)
        mass = 3.81/a
        if(TYPE == 'v'):
            a = a*-1

        ##Plot
        #black = all points, red = used points, blue = interpolation from polyfit over used points
        ax[row][col].set_title(TYPE + DIR)
        ax[row][col].text(0.4, 0.5, f'{mass:.2f}', transform=ax[row][col].transAxes)

        ax[row][col].plot(allXPts, allYPts, 'k.')
        ax[row][col].plot(extXPts, extYPts, 'r.')
        intXPts = ls(extXPts[0], extXPts[-1], 100)
        intYPts = [a*x*x + b*x + c for x in intXPts]
        ax[row][col].plot(intXPts, intYPts, 'b')


if(PLOT):
    plt.show()
