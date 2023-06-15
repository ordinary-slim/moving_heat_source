import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
import meshio
from wrapper import Problem, readInput, meshio_comparison

def mesh():
    #cell_type="quad4"
    meshDen = 1
    points, cells = meshzoo.ngon(9, 8)
    print( "nels = ", cells.shape[0] )
    cells = cells.astype( np.uint32 )
    return points, cells, "triangle3"

def specActive( p, center, d ):
    activeEls = np.zeros( p.domain.mesh.nels, dtype=int )
    for iel in range( p.domain.mesh.nels ) :
        e = p.domain.mesh.getElement( iel )
        centroid = e.getCentroid()

        if (np.linalg.norm( centroid  - center ) < d):
            activeEls[iel] = 1

    activeElements = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( activeElements )


def specWritePos( p ):
    # write post
    postMesh = meshio.Mesh(
        p.domain.mesh.pos,
        #[ ("triangle", p.mesh.con_CellPoint.con), ],
        [ ("line", p.domain.mesh.con_FacetPoint.con), ],
        #point_data={"BoundaryNodes": bNodes},
        cell_data={"isBoun":[p.domain.boundaryFacets.x]},
    )

    postMesh.write(
        "tmp.vtk",  # str, os.PathLike, or buffer/open file
    )

def run():
    inputFile = "input.txt"
    problemInput = readInput( inputFile )

    points, cells, cell_type = mesh()
    meshInput = {}
    meshInput["points"] = points
    meshInput["cells"] = cells
    meshInput["cell_type"]= cell_type
    m = mhs.Mesh( meshInput )

    p = mhs.Problem( m, problemInput )

    R0 = 0.0
    Rfinal = 1.0
    R = R0 + (5+1.0)/10*(Rfinal - R0)
    specActive(p, np.zeros(3), R)
    specWritePos( p )

def test():
    run()
    ref = "tmp_reference.vtk"
    new = "tmp.vtk"
    assert (meshio_comparison(ref, new))

if __name__=="__main__":
    test()
