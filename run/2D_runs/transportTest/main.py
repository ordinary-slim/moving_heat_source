import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../../Release/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem
import math
import pdb

def mesh(box, meshDen=4):
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], math.ceil(meshDen*(box[1]-box[0])+1)),
        np.linspace(box[2], box[3], math.ceil(meshDen*(box[3]-box[2])+1)),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def isInsideBox( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        pos = mesh.posFRF[el.con]
        xmin = min(pos[:, 0])
        xmax = max(pos[:, 0])
        ymin = min(pos[:, 1])
        ymax = max(pos[:, 1])

        isInside = 1*(xmin>=box[0] and xmax <= box[1] and ymin >= box[2] and ymax <= box[3])
        activeElements.append(isInside)
    return activeElements

def printSolutionActiveNodes( p ):
    ''' Debugging function '''
    for inode in range(p.mesh.nnodes):
        if (p.mesh.activeNodes[inode]):
            print("T({},{}) = {}".format( *p.mesh.pos[inode][0:2], p.unknown.values[inode]) )

def setAdimR( adimR, p ):
    r = p.input["radius"]
    HeatSourceSpeedX = max( abs(p.input["HeatSourceSpeedX"]), abs(p.input["advectionSpeedX"]))
    HeatSourceSpeedY = max( abs(p.input["HeatSourceSpeedY"]), abs(p.input["advectionSpeedY"]))
    HeatSourceSpeedZ = max( abs(p.input["HeatSourceSpeedZ"]), abs(p.input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [HeatSourceSpeedX, HeatSourceSpeedY, HeatSourceSpeedZ] ) )
    return (adimR * r / speed)

if __name__=="__main__":
    inputFile = "input.txt"
    mrfBox = [-36, 36, -5, 5]
    frfBox = [-16, 16, -5, 5]
    adimR = 2

    mrfProblem = Problem("mrfProblem")
    mrfTransporter = Problem("mrfTransport")
    # INITIALIZATIONS
    # Mesh BG
    points, cells, cell_type = mesh(mrfBox, meshDen=4)
    mrfProblem.setMesh( points, cells, cell_type )

    for p in [mrfTransporter, mrfProblem]:
        p.parseInput( inputFile )
        setAdimR(adimR, p)

    # Mesh FG
    points, cells, cell_type = mesh(frfBox, meshDen=4)
    mrfTransporter.setMesh( points, cells, cell_type )

    # Set advection FG
    mrfProblem.input["isAdvection"] = 1
    mrfProblem.input["advectionSpeedX"] = -mrfProblem.input["HeatSourceSpeedX"]
    mrfProblem.input["HeatSourceSpeedX"] = 0

    for p in [mrfTransporter, mrfProblem]:
        p.initialize()

    # Manufactured Initial Condition
    f = lambda pos : abs(pos[0]+pos[1])


    # CORRECT
    #mrfTransporter.unknown.getFromExternal(  mrfProblem.unknown )
    # DEBUGGING
    mrfTransporter.forceState( f )

    # MARCH FG, INTERPOLATE BACK TO BG
    for it in range(mrfTransporter.input["maxIter"]):
        mrfProblem.updateFRFpos()
        mrfProblem.unknown.getFromExternal( mrfTransporter.unknown )
        activeElements = isInsideBox( mrfProblem.mesh, frfBox )
        mrfProblem.activate( activeElements )

        mrfProblem.fakeIter()

        mrfTransporter.fakeIter()
        mrfTransporter.unknown.getFromExternal( mrfProblem.unknown )
        mrfProblem.writepos()
        mrfTransporter.writepos()
