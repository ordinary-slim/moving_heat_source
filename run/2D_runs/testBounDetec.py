import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
import meshio
from wrapper import Problem
import pdb

def mesh(box):
    cell_type="quad4"
    meshDen = 1
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    print( "nels = ", cells.shape[0] )
    cells = cells.astype( int )
    return points, cells

if __name__=="__main__":
    inputFile = "input.txt"
    box = [-2, 2, -2, 2]

    p     = Problem("FRF")

    points, cells = mesh(box)
    p.input["points"] = points
    p.input["cells"] = cells
    p.input["cell_type"]="quad4"

    #read input
    p.parseInput( inputFile )

    p.initialize()
    owner = p.mesh.findOwnerElement( np.array( [0, 0, 0] ) )

    #build boun nodes
    bNodes = np.zeros( p.mesh.nnodes )
    #pdb.set_trace()
    bFacets = np.zeros(p.mesh.con_FacetCell.con.shape[0])
    for ifacet in range( p.mesh.con_FacetCell.con.shape[0] ):
        bFacets[ifacet] = 1*(len([icell for icell in p.mesh.con_FacetCell.con[ifacet]
                                 if not(icell==-1)]) == 1)


    # write post
    postMesh = meshio.Mesh(
        p.mesh.pos,
        #[ ("quad", p.mesh.con_CellPoint.con), ],
        [ ("line", p.mesh.con_FacetPoint.con), ],
        #point_data={"BoundaryNodes": bNodes},
        cell_data={"bEdges":[bFacets]},
    )

    postMesh.write(
        "tmp.vtu",  # str, os.PathLike, or buffer/open file
    )
