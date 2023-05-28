import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
from wrapper import Problem, readInput, meshio_comparison
import numpy as np

Tfinal = 5.0

def mesh(leftEnd, rightEnd, elDen):
    cell_type="line2"
    nels = int((rightEnd-leftEnd)*elDen)
    points = np.linspace( leftEnd, rightEnd, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    auxCells = []
    for iel in range(nels):
        auxCells.append( [iel, iel+1] )
    cells = np.array( auxCells, dtype=int )
    return points, cells, cell_type

def run():
    inputFile = "input.txt"
    # read input
    problemInput = readInput( inputFile )

    # MESH
    elDen = 1
    leftEnd = -25.0
    rightEnd = +25.0
    boxDomain = [leftEnd, rightEnd]
    meshDict = {}
    meshDict["points"], meshDict["cells"], meshDict["cell_type"] = mesh(boxDomain[0], boxDomain[1], elDen)
    myMesh = mhs.Mesh(meshDict)

    p = Problem(myMesh, problemInput, caseName="neumann")

    # Neumann condition
    neumannVal = 10
    planeNormal  = np.array( [-1.0, 0.0, 0.0] )
    pointInPlane = np.array( [-25.0, 0.0, 0.0] )
    p.setNeumann( pointInPlane, planeNormal, neumannVal )
    planeNormal  = np.array( [+1.0, 0.0, 0.0] )
    pointInPlane = np.array( [+25.0, 0.0, 0.0] )
    p.setNeumann( pointInPlane, planeNormal, neumannVal )


    # FORWARD
    while (p.time < Tfinal-1e-7):
        p.iterate()#assembly + solve
    p.writepos()


def test():
    run()
    refds = "post_neumann_reference/neumann_1.vtu"
    newds = "post_neumann/neumann_1.vtu"

    # COMPARISON
    assert meshio_comparison(refds, newds)

if __name__=="__main__":
    test()
