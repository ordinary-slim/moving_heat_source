import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
from wrapper import Problem
import numpy as np
import meshzoo
import pdb
import meshio

speed = 10
bumper = 20
Tfinal = bumper/speed/2

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

def isInside( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        pos = mesh.posFRF[el.con]
        xmin = min(pos[:, 0])
        xmax = max(pos[:, 0])

        isInside = 1*(xmin>=box[0] and xmax <= box[1])
        activeElements.append(isInside)
    return activeElements

if __name__=="__main__":
    inputFile = "input.txt"

    p = Problem("myDirichletProblem")

    elDen = 4
    leftEnd = -50.0
    rightEnd = +50.0
    boxDomain = [leftEnd, rightEnd]
    # Meshing
    points, cells, cell_type = mesh(boxDomain[0], boxDomain[1], elDen)
    for p in [p,]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]=cell_type
        #p.input["numberOfGaussPointsCells"]=3
    # Dirichlet condition
    dirichletNodes = [0, len(points)-1]
    dirichletValues = [0.5, 0.5]

    for p in [p,]:
        p.input["dirichletNodes"] = dirichletNodes
        p.input["dirichletValues"] = dirichletValues

    #read input
    for p in [p,]:
        p.parseInput( inputFile )

    for p in [p,]:
        p.initialize()

    # FORWARD
    while (p.time < Tfinal):
        p.iterate()#assembly + solve
        p.writepos()
