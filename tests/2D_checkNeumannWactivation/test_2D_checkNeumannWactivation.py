import MovingHeatSource as mhs
import numpy as np
import meshzoo

maxIter = 1
dt = 1
Tfinal = 5
def mesh(box, meshDen=1):
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( np.uint32 )
    return points, cells, cell_type

def specActivate( p ):
    activeEls = []
    for iel in range( p.domain.mesh.nels ):
        e = p.domain.mesh.getElement( iel )
        centroid = e.getCentroid()

        if (centroid[0] >= 0.0):
            activeEls.append( iel )

    activeEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    return activeEls

def run():
    inputFile = "input.yaml"
    box = [-5, 5, -5, 5]

    # Read input
    problemInput = mhs.readInput( inputFile )
    # Mesh
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = mesh(box, meshDen=1)
    m = mhs.Mesh( meshInput )

    # Initialize problems
    p  = mhs.Problem(m, problemInput, caseName="2d_neumannWactivation")

    activeElements = specActivate( p )
    p.domain.setActivation( activeElements )

    def flux(point):
        return np.array([-10.0, 0.0, 0.0])
    p.setNeumann( p.domain.justActivatedBoundary.getIndices(), flux )

    it = 0
    time = 0.0
    while (( time < Tfinal) and ( it < maxIter )):
        time += dt
        it += 1
        p.iterate()
    p.writepos()

def test():
    run()
    refds = "post_2d_neumannWactivation_reference/2d_neumannWactivation_1.vtu"
    newds = "post_2d_neumannWactivation/2d_neumannWactivation_1.vtu"

    assert mhs.meshio_comparison( refds, newds )

if __name__=="__main__":
    test()

