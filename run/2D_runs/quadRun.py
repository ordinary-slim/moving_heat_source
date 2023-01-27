import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem

def mesh():
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(-10, 10, 41),
        np.linspace(-5, 5, 21),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    return points, cells

def isInsideBox( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        xmin = min(el.pos[:, 0])
        xmax = min(el.pos[:, 0])
        ymin = min(el.pos[:, 1])
        ymax = min(el.pos[:, 1])

        activeElements.append(1*(xmin>=box[0] and xmax <= box[1] and ymin >= box[2] and ymax <= box[3]))
    return activeElements

if __name__=="__main__":
    inputFile = "input.txt"
    postFolder = "postTrial"

    problem = Problem("firstRun")
    problem.parseInput( inputFile )
    points, cells = mesh()
    problem.input["points"] = points
    problem.input["cells"] = cells
    problem.input["cell_type"]="quad4"
    problem.input["isAdvection"] = 1
    problem.input["advectionSpeedX"] = -problem.input["speedX"]
    problem.input["speedX"] = 0.0

    problem.initialize()

    activeElements = isInsideBox( problem.mesh, [-5, 5, -5, 5] )
    problem.activate( activeElements )

    for iteration in range( problem.input["maxIter"] ):
        problem.iterate()
        problem.writepos(postFolder)
