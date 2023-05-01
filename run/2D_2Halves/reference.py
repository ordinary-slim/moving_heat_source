import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
import numpy as np
from wrapper import Problem
import pdb

def meshTri(box, meshDen=1, variant="up"):
    '''
    Variant = "zigzag",  or "up", "down", "center"
    '''
    cell_type="triangle3"
    points, cells = meshzoo.rectangle_tri(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        #cell_type=cell_type,
        variant=variant
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def meshQuad(box, meshDen=1, variant="up"):
    '''
    Variant useless here
    '''
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type="quad4",
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def exactSol22( x ):
    return (1 - np.power(x[0], 2) - np.power(x[1],2))

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

if __name__=="__main__":
    inputFile = "input.txt"
    box = [-1, 1, -1, 1]

    p  = Problem("reference")

    points, p.input["cells"], p.input["cell_type"] = meshQuad(box, meshDen=2, variant="zigzag")
    p.input["points"] = points

    #set Dirichlet BC. boundary nodes to 0
    dirichletNodes = []
    dirichletValues = []
    tol = 1e-7
    for inode in range(points.shape[0]):
        pos = points[inode, :]
        if (abs(abs(pos[0]) - 1) < tol) or (abs(abs(pos[1]) - 1) < tol):
            dirichletNodes.append( inode )
            dirichletValues.append( exactSol22(pos) )
    for p in [p,]:
        p.input["dirichletNodes"] = dirichletNodes
        p.input["dirichletValues"] = dirichletValues

    #read input
    for p in [p,]:
        p.parseInput( inputFile )

    for p in [p,]:
        p.initialize()


    #solve
    p.iterate()

    # post
    p.writepos()
