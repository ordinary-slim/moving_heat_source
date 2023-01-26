import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
import pdb
from prepos import PrePosProcessor

def mesh():
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(-5, 5, 21),
        np.linspace(-5, 5, 21),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    return points, cells

if __name__=="__main__":
    inputFile = "input.txt"
    postFolder = "post"

    wrapper = PrePosProcessor("firstRun", mhs.Problem())
    wrapper.parseInput( inputFile )
    points, cells = mesh()
    wrapper.input["points"] = points
    wrapper.input["cells"] = cells
    wrapper.input["cell_type"]="quad4"

    wrapper.initialize()

    for iteration in range( wrapper.input["maxIter"] ):
        wrapper.iterate()
        wrapper.writepos(postFolder)
