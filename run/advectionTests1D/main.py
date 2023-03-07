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
    return points, cells, cell_type

def setAdimR( adimR, p ):
    r = p.input["radius"]
    speedX = max( abs(p.input["speedX"]), abs(p.input["advectionSpeedX"]))
    speedY = max( abs(p.input["speedY"]), abs(p.input["advectionSpeedY"]))
    speedZ = max( abs(p.input["speedZ"]), abs(p.input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [speedX, speedY, speedZ] ) )
    return (adimR * r / speed)

if __name__=="__main__":
    inputFile = "input.txt"

    adimR = 1

    p = Problem("FRF")

    nels = 500
    leftEnd = 0.0
    rightEnd = 100.0
    points, cells, cell_type = mesh(leftEnd, rightEnd, nels)
    p.input["points"] = points
    p.input["cells"] = cells
    p.input["cell_type"]=cell_type

    #read input
    p.parseInput( inputFile )

    # set dt
    dt = setAdimR( adimR, p )
    for p in [p]:
        p.input["dt"] = dt

    for p in [p]:
        p.initialize()

    # Different IC
    # Manufactured Initial Condition
    #f = lambda pos : abs(pos[0]+pos[1])
    #f = lambda pos : 25
    #for p in [p,]:
        #p.forceState( f )

    maxIter = p.input["maxIter"]

    # FORWARD
    for iteration in range(maxIter):
        ##fine problem
        #for istep in range(fineStepsPerStep):
            #pFineFRF.iterate()
        #pFineFRF.writepos()

        #iter FRF
        p.iterate()#assembly + solve
        p.writepos()
