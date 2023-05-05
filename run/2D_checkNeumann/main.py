import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem, readInput
import pdb

maxIter = 50
dt = 0.5
Tfinal = 20
def mesh(box, meshDen=1):
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def specActivate( p ):
    activeEls = np.ones( p.domain.mesh.nels, dtype=int )
    for iel in range( p.domain.mesh.nels ):
        e = p.domain.mesh.getElement( iel )
        centroid = e.getCentroid()

        if (centroid[0] < 0.0):
            activeEls[iel] = 0

    activeEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    return activeEls

if __name__=="__main__":
    inputFile = "input.txt"
    box = [-16, 16, -5, 5]

    # Read input
    problemInput = readInput( inputFile )
    # Mesh
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = mesh(box, meshDen=2)
    m = mhs.Mesh( meshInput )

    # Initialize problems
    p  = Problem(m, problemInput, caseName="case")

    activeElements = specActivate( p )
    p.domain.setActivation( activeElements )

    pointInPlane = np.array([0.0, 0.0, 0.0])
    planeNormal  = np.array([-1.0, 0.0, 0.0])
    def flux(point):
        return np.array([-10.0, 0.0, 0.0])
    p.setNeumann( p.domain.justActivatedBoundary.getTrueIndices(), flux )

    it = 0
    time = 0.0
    while (( time < Tfinal) and ( it < maxIter )):
        time += dt
        it += 1
        p.iterate()
        p.writepos()
