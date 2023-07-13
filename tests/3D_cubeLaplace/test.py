import MovingHeatSource as mhs
import numpy as np
import meshzoo

def exactSol( x ):
    return (1 - np.power(x[0], 2) - np.power(x[1],2) - np.power(x[2],2))

def setDirichlet( p ):
    #set Dirichlet BC. boundary nodes to 0
    dirichletNodes = []
    dirichletValues = []
    tol = 1e-7
    for inode in range(p.domain.mesh.nnodes):
        pos = p.domain.mesh.pos[inode, :]
        if (abs(abs(pos[0]) - 1) < tol) or (abs(abs(pos[1]) - 1) < tol) or (abs(abs(pos[2]) - 1) < tol):
            dirichletNodes.append( inode )
            dirichletValues.append( exactSol(pos) )
    p.setDirichlet( dirichletNodes, dirichletValues )

def mesh(box, meshDen=4):
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

def run():
    inputFile = "input.txt"
    problemInput = mhs.readInput( inputFile )

    box = [-1, 1, -1, 1, -1, 1]
    meshDen = 4
    p, c, cell_type = mesh(box, meshDen)
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = p, c, cell_type
    myMesh = mhs.Mesh( meshInput )

    p = mhs.Problem( myMesh, problemInput )

    print("Setting BCs...")
    setDirichlet( p )

    # Pre-assembly, updating free dofs
    p.preAssemble(False)

    p.assemble()

    p.ls.solve()

    p.gather()

    p.postIterate()

    fexact = p.project( exactSol )
    p.writepos(
            functions={
                "fexact":fexact,
                },
            )
def test():
    run()
    reference = "post_case_reference/case_0.vtu"
    new = "post_case/case_0.vtu"
    assert mhs.meshio_comparison(reference, new)

if __name__=="__main__":
    test()
