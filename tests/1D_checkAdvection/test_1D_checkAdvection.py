import MovingHeatSource as mhs
import numpy as np

def mesh(leftEnd, rightEnd, nels):
    cell_type="line2"
    points = np.linspace( leftEnd, rightEnd, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    auxCells = []
    for iel in range(nels):
        auxCells.append( [iel, iel+1] )
    cells = np.array( auxCells, dtype=int )
    numberOfGaussPointsCells = 3
    return points, cells, cell_type, numberOfGaussPointsCells

def run():
    inputFile = "input.yaml"
    stableProblemInput = mhs.readInput( inputFile )
    unstableProblemInput = dict( stableProblemInput )

    # MESHING
    nels = 100
    leftEnd = -50.0
    rightEnd = +50.0
    points, cells, cell_type, ngpoins = mesh(leftEnd, rightEnd, nels)

    meshInput = {}
    meshInput["points"] = points
    meshInput["cells"] = cells
    meshInput["cell_type"]=cell_type
    meshInput["numberOfGaussPointsCells"] = 3
    m = mhs.Mesh( meshInput )

    stableProblemInput["isStabilized"] = 1
    unstableProblemInput["isStabilized"] = 0

    pStable = mhs.Problem(m, stableProblemInput, caseName="stable")
    pUnstable = mhs.Problem(m, unstableProblemInput, caseName="unstable")

    # Manufactured Initial Condition
    L = rightEnd - leftEnd
    center = (rightEnd + leftEnd) / 2
    f = lambda pos : 2*(pos[0] >= center)
    for p in [pUnstable, pStable]:
        p.forceState( f )

    maxIter = pStable.input["maxIter"]

    while ( pStable.time < pStable.input["Tfinal"] ):
        for p in [pStable, pUnstable]:
            p.iterate()#assembly + solve

    # Write post
    for p in [pStable, pUnstable]:
        p.writepos()

def test():
    run()
    refds = "post_stable_reference/stable_5.vtu"
    newds = "post_stable/stable_5.vtu"

    assert mhs.meshio_comparison( refds, newds, 1e-5 )

if __name__=="__main__":
    test()

