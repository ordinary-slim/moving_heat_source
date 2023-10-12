import MovingHeatSource as mhs
import numpy as np
from MovingHeatSource.adaptiveStepper import meshRectangle, meshBox

def run2d():
    inputFile = "input.yaml"
    box = np.array([-1, 1, 0, 1])

    # read input
    problemInput = mhs.readInput( inputFile )

    problemInput = dict( problemInput )

    # Mesh
    myMesh = meshRectangle( box, elSize=0.25 )

    p = mhs.Problem(myMesh, problemInput, caseName="collisions_2d")

    origin = np.array( [-1.0, 0.0, 0.0] )
    destin = np.array( [+1.0,+1.0, 0.0] )
    mdwidth = 0.1
    mdheight = 0.1

    printer = mhs.Printer( p, mdwidth, mdheight/2, mdheight/2 )
    collidedEls = mhs.MeshTag( myMesh, myMesh.dim )
    printer.deposit( origin, destin, collidedEls )

    p.writepos(
            cellMeshTags={"collidedEls" : collidedEls},
            )

def run3d():
    inputFile = "input.yaml"
    box = np.array([-1, 1, 0, 1, -1, 1])

    # read input
    problemInput = mhs.readInput( inputFile )

    problemInput = dict( problemInput )

    # Mesh
    myMesh = meshBox( box, elSize=0.25 )

    p = mhs.Problem(myMesh, problemInput, caseName="collisions_3d")

    origin = np.array( [-1.0, 0.0, -1.0] )
    destin = np.array( [+1.0, 1.0, +1.0] )
    mdwidth = 0.01
    mdheight = 0.01

    printer = mhs.Printer( p, mdwidth, mdheight/2, mdheight/2 )
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

if __name__=="__main__":
    test()
