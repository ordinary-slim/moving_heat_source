import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
import numpy as np
from wrapper import Problem, readInput
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

def exactSol( x ):
    return (1 - np.power(x[0], 2) - np.power(x[1],2))

def setDirichlet( p ):
    #set Dirichlet BC. boundary nodes to 0
    dirichletNodes = []
    dirichletValues = []
    tol = 1e-7
    for inode in range(p.domain.mesh.nnodes):
        pos = p.domain.mesh.pos[inode, :]
        if (abs(abs(pos[0]) - 1) < tol) or (abs(abs(pos[1]) - 1) < tol):
            dirichletNodes.append( inode )
            dirichletValues.append( exactSol(pos) )
    p.setDirichlet( dirichletNodes, dirichletValues )

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

    # Read input
    problemInput = readInput( inputFile )

    # Mesh
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = meshQuad(box, meshDen=2, variant="zigzag")
    m = mhs.Mesh( meshInput )

    # Initialize problems
    p  = Problem(m, problemInput, caseName="reference")

    # Set Dirichlet
    setDirichlet( p )

    #solve
    p.iterate()

    #exact sol
    vals = np.zeros( p.domain.mesh.nnodes )
    for inode in range( p.domain.mesh.nnodes):
        pos = p.domain.mesh.pos[inode, :]
        vals[inode] = exactSol( pos )
    es = mhs.Function( p.domain.mesh, vals )


    # post
    p.writepos(functions={
        "exactSol":es,
        })
