import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import matplotlib.pyplot as plt
from myPlotHandler import myPlotHandler
import numpy as np

if __name__=="__main__":
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    P = d["power"]
    R = d["radius"]
    x0 = (d["Right"] + d["Left"])/2
    r = lambda x : 2*P/np.pi/R*np.exp( -((x - x0)/R**2)**2 )

    xRange = np.linspace( d["Left"], d["Right"], 1000 )

    plt.figure(dpi=200)
    ax = plt.gca()
    plt.plot( xRange, r(xRange),
            label="Gaussian source term",
            )
    plt.xlabel(r"x")
    plt.ylabel(r"Linear heat source")
    ax.text(0, 100, 'P = {}\nR= {}\n$x_0 = {}$'.format( P, R, x0 ),
            bbox=dict(boxstyle="square,pad=0.3", fc="w", ec="k",lw=2))
    ax.text(40, 280, r'$r(x, t) = \frac{2 P }{\pi R^2} \exp( -(\frac{x-x_0(t)}{R^2})^2 )$',
            fontsize=18,
            bbox=dict(boxstyle="square,pad=0.3", fc="w", ec="k",lw=2))
    plt.show()
