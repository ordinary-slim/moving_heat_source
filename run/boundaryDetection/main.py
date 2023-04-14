import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
import meshio
from wrapper import Problem
import pdb

def mesh():
    #cell_type="quad4"
    meshDen = 1
    points, cells = meshzoo.ngon(9, 8)
    print( "nels = ", cells.shape[0] )
    cells = cells.astype( int )
    return points, cells

def specActive( p, center, d ):
    activeEls = np.zeros( p.mesh.nels, dtype=int )
    for iel in range( p.mesh.nels ) :
        e = p.mesh.getElement( iel )
        centroid = e.getCentroid()

        if (np.linalg.norm( centroid  - center ) < d):
            activeEls[iel] = 1

    return activeEls


def specWritePos( p, iter ):
    auxIsBoun = np.zeros( p.mesh.con_FacetPoint.nels_oDim )
    auxIsBoun[ p.mesh.boundaryFacets ] = 1.0
    # write post
    postMesh = meshio.Mesh(
        p.mesh.pos,
        #[ ("triangle", p.mesh.con_CellPoint.con), ],
        [ ("line", p.mesh.con_FacetPoint.con), ],
        #point_data={"BoundaryNodes": bNodes},
        cell_data={"isBoun":[auxIsBoun]},
    )

    postMesh.write(
        "tmp_{}.vtu".format( iter ),  # str, os.PathLike, or buffer/open file
    )

if __name__=="__main__":
    inputFile = "input.txt"

    p     = Problem("FRF")

    points, cells = mesh()
    p.input["points"] = points
    p.input["cells"] = cells
    p.input["cell_type"]="triangle3"

    #read input
    p.parseInput( inputFile )

    p.initialize()

    '''
    isPointInside = np.zeros( p.mesh.nels )
    point = np.array(
            [0.75,
             0.42,
             0.0] )
    owner = p.mesh.findOwnerElement( point )
    if owner >= 0:
        isPointInside[owner] = 1
    '''

    numIter = 10
    R0 = 0.0
    Rfinal = 1.0
    for iter in range( 10 ):
        R = R0 + (iter+1.0)/numIter*(Rfinal - R0)
        activeEls = specActive(p, np.zeros(3), R)
        p.activate( activeEls )
        specWritePos( p, iter )
