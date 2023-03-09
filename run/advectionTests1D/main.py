import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
from wrapper import Problem
import numpy as np
import meshzoo
import pdb
import meshio

def mesh(leftEnd, rightEnd, nels):
    cell_type="line2"
    points = np.linspace( leftEnd, rightEnd, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    auxCells = []
    for iel in range(nels):
        auxCells.append( [iel, iel+1] )
    cells = np.array( auxCells, dtype=int )
    numberOfGaussPoints = 3
    return points, cells, cell_type, numberOfGaussPoints

if __name__=="__main__":
    inputFile = "input.txt"

    pUnstable = Problem("onlyAdvectionUnstable")
    pStable = Problem("onlyAdvectionStable")

    nels = 1000
    leftEnd = -50.0
    rightEnd = +50.0
    points, cells, cell_type, ngpoins = mesh(leftEnd, rightEnd, nels)

    for p in [pUnstable, pStable]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]=cell_type
        p.input["numberOfGaussPoints"]=3

    #read input
    for p in [pUnstable, pStable]:
        p.parseInput( inputFile )

    pUnstable.input["isStabilized"] = 0
    pStable.input["isStabilized"] = 1

    #set advection
    for p in [pUnstable, pStable]:
        p.input["advectionSpeedX"] = 5.0

    for p in [pUnstable, pStable]:
        p.initialize()

    # Manufactured Initial Condition
    L = rightEnd - leftEnd
    center = (rightEnd + leftEnd) / 2
    f = lambda pos : 2*(pos[0] >= center)
    for p in [pUnstable, pStable]:
        p.forceState( f )

    maxIter = pStable.input["maxIter"]

    # FORWARD
    while ( pStable.time < pStable.input["Tfinal"] ):
        ##fine problem
        #for istep in range(fineStepsPerStep):
            #pFineFRF.iterate()
        #pFineFRF.writepos()

        for p in [pStable, pUnstable]:
            p.iterate()#assembly + solve
            print("max T = {}".format( max(p.unknown.values) ))
            p.writepos()
