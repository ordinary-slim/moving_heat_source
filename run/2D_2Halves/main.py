import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem
import pdb

def mesh(box, meshDen=1, variant="up"):
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

def exactSol22( x ):
    return (1 - np.power(x[0], 2) - np.power(x[1],2))

def setDirichlet22( p ):
    points = p.input["points"]
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

    activeElements = mhs.MeshTag( mesh, mesh.dim, activeElements )
    return activeElements

if __name__=="__main__":
    inputFile = "input.txt"
    box = [-1, 1, -1, 1]
    boxLeft = [-1, 0, -1, 1]

    pLeft  = Problem("left")
    pRight  = Problem("right")

    pLeft.input["points"], pLeft.input["cells"], pLeft.input["cell_type"] = mesh(box, meshDen=2, variant="zigzag")
    pRight.input["points"], pRight.input["cells"], pRight.input["cell_type"] = mesh(box, meshDen=3, variant="up")

    meshLeft = mhs.Mesh()
    meshLeft.initializeMesh( pLeft.input )
    meshRight = mhs.Mesh()
    meshRight.initializeMesh( pRight.input )

    #read input
    for p in [pLeft, pRight,]:
        p.parseInput( inputFile )

    # set dirichlet BC
    for p in [pLeft, pRight,]:
        setDirichlet22( p )

    pLeft.initialize(meshLeft)
    pRight.initialize(meshRight)

    # Deactive left of pLeft
    leftActiveEls = isInsideBox( pLeft.domain.mesh, boxLeft )
    pLeft.domain.setActivation( leftActiveEls )

    pRight.deactivateFromExternal( pLeft )

    # Set up left problem with dirichlet condition at right boundary

    # post
    pLeft.writepos()
    pRight.writepos(functions={
        }
                    )
