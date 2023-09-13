import MovingHeatSource as mhs
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

inputFile = "input.yaml"

# Read input
problemInput = mhs.readInput( inputFile )
L = problemInput["L"]
Tfinal = problemInput["Tfinal"]
nels = problemInput["nels"]
maxIter = problemInput["maxIter"]


def getMesh(elementSize=0.1):
    cell_type="line2"
    nels = int(L / elementSize)
    points = np.linspace( -L/2, +L/2, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    cells = np.transpose( np.vstack( (np.arange(0, nels), np.arange(1, nels+1) ), dtype=int ) )
    dimension = 1
    return mhs.Mesh( {
                    "points" : points,
                    "cells" : cells,
                    "cell_type" : cell_type,
                    "dimension" : 1,
                    })

def plotProblem( p ):
    # Default settings
    mpl.rcParams.update(mpl.rcParamsDefault)

    plt.style.use("bmh")
    mpl.rcParams.update({'font.size':22})

    #sns.set_palette( sns.color_palette("bright") )
    plt.plot( p.domain.posLab[:, 0],
             p.unknown.values,
             #label="T",
             )

    plt.gca().set_xlim( [p.domain.posLab[0, 0],
                         p.domain.posLab[-1,0],] )
    plt.xlabel("x")
    plt.ylabel("u")
    #plt.legend()
    plt.pause(0.2)



def main(plot=False):
    mesh = getMesh(elementSize=(L/nels))

    # Initialize problems
    p  = mhs.Problem(mesh, problemInput)

    for _ in range(maxIter):
        p.iterate()
        p.writepos()
        if plot:
            plt.clf()
            plotProblem( p )
    if plot:
        plt.show()

if __name__=="__main__":
    main()
