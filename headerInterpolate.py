#Header file for interpolating (densely packed) grid-sampled functions

#Given points 0 and f, returns the linear interpolated point at xi; x0 < xi < xf
def LinInterp(x0, xf, y0, yf, xi):
    if(x0 < xi < xf):
        slope = (yf-y0)/(xf-x0)
        intercept = yf - slope*xf
        return slope*xi + intercept
    else:
        print("LinInterp: Bad input...")
        return None

#Given a min, max value + spacing, returns an evenly spaced grid of points (to eval a function on)
#min and max are guarenteed to be in the list
#Spacing should roughly equal the spacing between the known sample function's x points
def LinSpace(min, max, step=0.1):
    ret = [float(min)]
    for i in range(0, int((max-min)/step)):
        ret.append(ret[-1] + step)
    return ret

#Returns the index of an eval grid whose value is closest to x
#i.e. if min, max, step = 1, 3, 0.25 then the eval grid is [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0].
#x = 1.374 would return 1, but x = 1.376 would return 2
def NearestIndex(min, max, step, x):
    return round((x - min)/step)

#Given a sampled x, y and an eval point xi, returns the two points closest to xi
#Return is xBelow, yBelow, xAbove, yAbove
def GetNeighbors(x, y, xi):
    nearestBelowX, nearestAboveX = x[0], x[-1]
    bestBelowDiffX, bestAboveDiffX = xi - nearestBelowX, nearestAboveX - xi
    nearestBelowY, nearestAboveY = y[0], y[-1]

    for i in range(1, len(x) - 1):
        if(x[i] < xi and xi - x[i] < bestBelowDiffX):
            bestBelowDiffX = xi - x[i]
            nearestBelowX = x[i]
            nearestBelowY = y[i]

        if(xi < x[i] and x[i] - xi < bestAboveDiffX):
            bestAboveDiffX = x[i] - xi
            nearestAboveX = x[i]
            nearestAboveY = y[i]

    return nearestBelowX, nearestBelowY, nearestAboveX, nearestAboveY

#Takes a function (defined in x, y) as well as LinSpace() arguments to evaluate the function on an evenly
#spaced eval grid
#Returns the grid itself and then the y values associated with it
def ReEval(x, y, min, max, step):
    assert x[0] < min
    assert max < x[-1]

    gridX = LinSpace(min=min, max=max, step=step)
    gridY = [0. for i in range(0, len(gridX))]
    nInterps = [0 for i in range(0, len(gridX))] ##number of times grid's index has been linear interped

    for i in range(0, len(gridX)):
        x0, y0, xf, yf = GetNeighbors(x=x, y=y, xi=gridX[i]) ##get neighbors for this grid point
        gridY[i] += LinInterp(x0=x0, xf=xf, y0=y0, yf=yf, xi=gridX[i])
        nInterps[i] += 1


    return gridX, [gridY[i]/float(nInterps[i]) for i in range(0, len(gridY))]
