import MovingHeatSource as mhs
import numpy as np
import meshzoo

def mesh2d(box, meshDen=4):
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def mesh3d(box, meshDen=4):
    cell_type="hexa8"
    nelsX = int((box[1] - box[0])*meshDen)
    nelsY = int((box[3] - box[2])*meshDen)
    nelsZ = int((box[5] - box[4])*meshDen)

    points, cells = meshzoo.cube_hexa(
        np.linspace( box[0], box[1], nelsX+1),
        np.linspace( box[3], box[2], nelsY+1),
        np.linspace( box[5], box[4], nelsZ+1),
    )
    cells = cells.astype( np.uint32 )
    return points, cells, cell_type


def run2d():
    inputFile = "input.yaml"
    box = [-1, 1, 0, 1]

    # read input
    problemInput = mhs.readInput( inputFile )

    problemInput = dict( problemInput )

    # Mesh
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = mesh2d(box)

    myMesh = mhs.Mesh(meshInput)

    p = mhs.Problem(myMesh, problemInput, caseName="collisions_2d")

    origin = np.array( [-1.0, 0.0, 0.0] )
    destin = np.array( [+1.0,+1.0, 0.0] )
    mdwidth = 0.1
    mdheight = 0.1

    printer = mhs.Printer( p, mdwidth, mdheight )
    collidedEls = mhs.MeshTag( myMesh, myMesh.dim )
    printer.deposit( origin, destin, collidedEls )

    p.writepos(
            cellMeshTags={"collidedEls" : collidedEls},
            )

def run3d():
    inputFile = "input.yaml"
    box = [-1, 1, 0, 1, -1, 1]

    # read input
    problemInput = mhs.readInput( inputFile )

    problemInput = dict( problemInput )

    # Mesh
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = mesh3d(box, meshDen=4)

    myMesh = mhs.Mesh(meshInput)

    p = mhs.Problem(myMesh, problemInput, caseName="collisions_3d")

    origin = np.array( [-1.0, 0.0, -1.0] )
    destin = np.array( [+1.0, 1.0, +1.0] )
    mdwidth = 0.01
    mdheight = 0.01

    printer = mhs.Printer( p, mdwidth, mdheight )
    collidedEls = mhs.MeshTag( myMesh, myMesh.dim )
    printer.deposit( origin, destin, collidedEls )

    p.writepos(
            cellMeshTags={"collidedEls" : collidedEls},
            )

def test():
    run2d()
    run3d()

    refDSs = ["post_collisions_2d_reference/collisions_2d_0.vtu",
              "post_collisions_3d_reference/collisions_3d_0.vtu"]
    newDSs = ["post_collisions_2d/collisions_2d_0.vtu",
              "post_collisions_3d/collisions_3d_0.vtu"]

    tests = []
    for refds, newds in zip(refDSs, newDSs):
        tests.append( mhs.meshio_comparison( refds, newds,
                                            psets=[],
                                            csets=["collidedEls"] ) )

    assert all(tests)
