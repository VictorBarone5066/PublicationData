from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
from numpy import linspace as ls
from numpy import polyfit as pf

# Figure Options
rcParams["font.family"] = "Times New Roman"
rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
rcParams["figure.autolayout"] = False

rcParams['axes.linewidth'] = 1.5
rcParams["font.weight"] = "bold"
rcParams["axes.labelweight"] = "bold"

#Common constants
def WIDTH(nCol=1): ##inches
    if(nCol == 1):
        return 3.5
    if(nCol == 2):
        return 7.2
    return 3.5
def HEIGHT(nCol=1): ##inches
    if(nCol == 1):
        return 3.5
    if(nCol == 2):
        return 7.2
    return 3.5


GRID_ALPHA = 0.5

FONTSIZE_LABELS = 12
FONTSIZE_TICKS = 12
FONTSIZE_TEXT = 12

MKSIZE_SML = 1.50
MKSIZE_MED = 3.25
MKSIZE_BIG = 5.00

LNSIZE_SML = 1.5
LNSIZE_MED = 1.5
LNSIZE_BIG = 1.5

DPI = 1000

#Useful functions
def FloatFormat(d, dec=3):
    return f'{d:.{dec}f}'

def RmAxTicks(ax, x=True, y=True):
    if(x):
        ax.tick_params(axis='x', colors="white")
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(0.001)
    if(y):
        ax.tick_params(axis='y', colors="white")
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(0.001)
    return ax

def SetTickSize(ax, fontsize=FONTSIZE_TICKS, x=True, y=True):
    if(x):
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
    if(y):
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
    return ax

def PolyInterp(x, y, deg=1, nPts=1000):
    xI, yI = (list(a) for a in zip(*sorted(zip(x, y))))
    fit = pf(xI, yI, deg=deg)

    xI = ls(start=xI[0], stop=xI[-1], num=nPts, endpoint=True)
    yI = []
    for x_ in xI:
        sum = 0.
        for i in range(0, deg + 1):
            sum += fit[i]*x_**(deg - i)
        yI.append(sum)

    return xI, yI


#Common colors
COLORS_DICT =  {"red":    "#ff0000",  # red
                "orange": "#eb7c00",  # orange
                "lime":   "#93ff27",  # lime
                "dGreen": "#00b423",  # d green
                "black":  "black",    # black
                "cyan":   "#00fdea",  # cyan
                "blue":   "#0058fd",  # blue
                "purple": "#cb34ff",  # purple
                "pink":   "#ff04b7"}  # pink

#               red        blue       d green    black    orange     purple    lime        pink
COLORS_LIST = ["#ff0000", "#0058fd", "#00b423", "black", "#eb7c00", "#cb34ff", "#93ff27", "#ff04b7",
#               cyan
               "#00fdea"]

#Useful symbols
ANGSTROM = 'Å'

#Rudorffite Specific
RUD_AG = "#9D9D9D"
RUD_BI = "#2CADF3"
RUD_I = "#000000"
