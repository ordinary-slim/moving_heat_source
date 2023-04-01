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

    isPointInside = np.zeros( p.mesh.nels )
    point = np.array(
            [0.75,
             0.42,
             0.0] )
    owner = p.mesh.findOwnerElement( point )
    if owner >= 0:
        isPointInside[owner] = 1


    '''
    #build boun nodes
    bNodes = np.zeros( p.mesh.nnodes )
    #pdb.set_trace()
    bFacets = np.zeros(p.mesh.con_FacetCell.con.shape[0])
    for ifacet in range( p.mesh.con_FacetCell.con.shape[0] ):
        bFacets[ifacet] = 1*(len([icell for icell in p.mesh.con_FacetCell.con[ifacet]
                                 if not(icell==-1)]) == 1)
     '''


    # write post
    postMesh = meshio.Mesh(
        p.mesh.pos,
        [ ("triangle", p.mesh.con_CellPoint.con), ],
        #[ ("line", p.mesh.con_FacetPoint.con), ],
        #point_data={"BoundaryNodes": bNodes},
        cell_data={"isPointInside":[isPointInside]},
    )

    postMesh.write(
        "tmp.vtu",  # str, os.PathLike, or buffer/open file
    )
