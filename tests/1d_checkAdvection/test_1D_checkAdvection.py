import MovingHeatSource as mhs
import numpy as np

def getMesh(nels, L):
    cell_type="line2"
    points = np.linspace( -L/2, +L/2, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    cells = np.transpose( np.vstack( (np.arange(0, nels), np.arange(1, nels+1) ), dtype=int ) )
    numberOfGaussPointsCells = 3
    return mhs.Mesh( {
                    "points" : points,
                    "cells" : cells,
                    "cell_type" : cell_type,
                    "numberOfGaussPointsCells" : numberOfGaussPointsCells,
                    "dimension" : 1,
                    })

def run1():
    inputFile = "input1.yaml"
    stableProblemInput = mhs.readInput( inputFile )
    unstableProblemInput = dict( stableProblemInput )

    # MESHING
    L = stableProblemInput["L"]
    m = getMesh(stableProblemInput["nels"], L)

    stableProblemInput["isStabilized"] = 1
    unstableProblemInput["isStabilized"] = 0

    pStable = mhs.Problem(m, stableProblemInput, caseName="stable")
    pUnstable = mhs.Problem(m, unstableProblemInput, caseName="unstable")

    # Manufactured Initial Condition
    f = lambda pos : 2*(pos[0] >= 0.0)
    for p in [pUnstable, pStable]:
        p.forceState( f )

    while ( pStable.time < pStable.input["Tfinal"] ):
        for p in [pStable, pUnstable]:
            p.iterate()#assembly + solve

    # Write post
    for p in [pStable, pUnstable]:
        p.writepos()

def run2():
    inputFile = "input2.yaml"

    # Read input
    problemInput = mhs.readInput( inputFile )
    L = problemInput["L"]
    Tfinal = problemInput["Tfinal"]
    nels = problemInput["nels"]
    maxIter = problemInput["maxIter"]

    mesh = getMesh(nels, L)

    # Initialize problems
    p  = mhs.Problem(mesh, problemInput, caseName="nonUnitParams" )

    for _ in range(maxIter):
        p.iterate()
    p.writepos()

def test1():
    run1()
    refds = "post_stable_reference/stable_5.vtu"
    newds = "post_stable/stable_5.vtu"

    assert mhs.meshio_comparison( refds, newds, tol=1e-5 )

def test2():
    run2()
    newds = "post_nonUnitParams/nonUnitParams_20.vtu"
    refds = "post_nonUnitParams_reference/nonUnitParams_20.vtu"
    assert mhs.meshio_comparison( refds, newds, tol=1e-5 )

if __name__=="__main__":
    test1()
    test2()
