import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release')
import MovingHeatSource as mhs
from wrapper import Problem, readInput, meshio_comparison
import numpy as np

def mesh(leftEnd, rightEnd, elDen):
    cell_type="line2"
    nels = int((rightEnd-leftEnd)*elDen)
    points = np.linspace( leftEnd, rightEnd, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    auxCells = []
    for iel in range(nels):
        auxCells.append( [iel, iel+1] )
    cells = np.array( auxCells, dtype=int )
    return points, cells, cell_type


def run(caseName):
    inputFile = "input.txt"

    problemInput = readInput( inputFile )
    # Mesh
    elDen = 1
    leftEnd = -20.0
    rightEnd = +20.0
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = mesh(leftEnd, rightEnd, elDen)
    m = mhs.Mesh( meshInput )

    # Dirichlet condition
    dirichletNodes = [0, m.nnodes-1]
    dirichletValues = [5, 5]

    problemInput["dirichletNodes"] = dirichletNodes
    problemInput["dirichletValues"] = dirichletValues
    problemInput["steadyState"] = 1

    # Initialize problems
    p  = Problem(m, problemInput, caseName=caseName)

    # Solve
    p.iterate()#assembly + solve
    p.writepos()

def test_1d_dirichlet():
    run("1d_dirichlet")
    referenceDs = "post_1d_dirichlet_reference/1d_dirichlet_reference_1.vtu"
    trialDs =  "post_1d_dirichlet/1d_dirichlet_1.vtu"

    # COMPARISON
    assert meshio_comparison(referenceDs, trialDs)
