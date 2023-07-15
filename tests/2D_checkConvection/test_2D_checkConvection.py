import MovingHeatSource as mhs
import numpy as np
import meshzoo

maxIter = 40
dt = 1
Tfinal = 5
initialT = 10

def mesh(box):
    cell_type="quad4"
    meshDen = 1
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( np.uint32 )
    return points, cells, cell_type

def run():
    inputFile = "input.yaml"
    box = [-10, 10, -5, 5]

    # Read input
    problemInput = mhs.readInput( inputFile )

    # Mesh
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = mesh(box)
    m = mhs.Mesh(meshInput)

    # Initialize problems
    p  = mhs.Problem(m, problemInput, caseName="2d_checkConvection")


    f = lambda pos : initialT
    p.forceState( f )

    it = 0
    time = 0.0
    while (( time < Tfinal) and ( it < maxIter )):
        time += dt
        it += 1
        p.iterate()
    p.writepos()


def test():
    run()
    referenceDs = "post_2d_checkConvection_reference/2d_checkConvection_5.vtu"
    trialDs =  "post_2d_checkConvection/2d_checkConvection_5.vtu"

    # COMPARISON
    assert mhs.meshio_comparison(referenceDs, trialDs)

if __name__=="__main__":
    test()
