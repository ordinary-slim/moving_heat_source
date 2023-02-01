import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem

def mesh(box):
    cell_type="quad4"
    meshDen = 4
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    return points, cells

def isInsideBox( mesh, box ):
    activeElements = []
    for ielem in range( mesh.nels ):
        el = mesh.getElement( ielem )
        pos = np.array(el.pos)
        for inode in range(el.nnodes):
            pos[inode] += mesh.x0
        xmin = min(pos[:, 0])
        xmax = max(pos[:, 0])
        ymin = min(pos[:, 1])
        ymax = max(pos[:, 1])

        activeElements.append(1*(xmin>=box[0] and xmax <= box[1] and ymin >= box[2] and ymax <= box[3]))
    return activeElements

if __name__=="__main__":
    inputFile = "input.txt"
    boxRef = [-5, 5, -5, 5]
    boxInac = [-10, 10, -10, 10]

    problemRef = Problem("reference")
    problemTest = Problem("test")

    points, cells = mesh(boxRef)
    problemRef.input["points"] = points
    problemRef.input["cells"] = cells
    points, cells = mesh(boxInac)
    problemTest.input["points"] = points
    problemTest.input["cells"] = cells

    for p in [problemRef, problemTest]:
        p.parseInput( inputFile )
        p.input["cell_type"]="quad4"
        p.input["isAdvection"] = 1
        p.input["advectionSpeedX"] = -p.input["speedX"]
        p.input["speedX"] = 0.0
        p.initialize()

        for iteration in range(5):
            activeElements = isInsideBox( p.mesh, boxRef )
            p.activate( activeElements )
            p.iterate()
            p.writepos()
