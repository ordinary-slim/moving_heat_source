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

    p = Problem("onlyAdvection")

    nels = 1000
    leftEnd = -50.0
    rightEnd = +50.0
    points, cells, cell_type, ngpoins = mesh(leftEnd, rightEnd, nels)
    p.input["points"] = points
    p.input["cells"] = cells
    p.input["cell_type"]=cell_type
    p.input["numberOfGaussPoints"]=3

    #read input
    p.parseInput( inputFile )

    #set advection
    p.input["advectionSpeedX"] = -5.0

    for p in [p]:
        p.initialize()

    # Manufactured Initial Condition
    L = rightEnd - leftEnd
    f = lambda pos : np.sin( 4 * pos[0] / L * np.pi )
    for p in [p,]:
        p.forceState( f )

    maxIter = p.input["maxIter"]

    # FORWARD
    for iteration in range(maxIter):
        ##fine problem
        #for istep in range(fineStepsPerStep):
            #pFineFRF.iterate()
        #pFineFRF.writepos()

        #iter FRF
        p.iterate()#assembly + solve
        print("max T = {}".format( max(p.unknown.values) ))
        p.writepos()
