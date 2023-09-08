import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
import numpy as np
import meshzoo
from CustomStepper import CustomStepper
import sys
import pdb

inputFile = "input.yaml"
tol = 1e-7
adimR_domain = 10 
problemInput = mhs.readInput( inputFile )
adimR = problemInput["radius"] / np.linalg.norm(problemInput["HeatSourceSpeed"])
boxPhys = [-10, 10, -10, 10]
gcodeFile = problemInput["path"]
radiusHs = problemInput["radius"]
fineElSize = radiusHs/4

def meshBox(box, elementSize=0.25):
    cell_type="quad4"
    nelsX = int((box[1]-box[0]) / elementSize) +1
    nelsY = int((box[3]-box[2]) / elementSize) +1
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], nelsX),
        np.linspace(box[2], box[3], nelsY),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def getMesh(elementSize=fineElSize):
    # Mesh
    meshInput, meshInputMoving = {}, {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = meshBox(boxPhys, elementSize=elementSize)
    return mhs.Mesh( meshInput )

def runReference():
    mesh = getMesh()
    pFRF     = mhs.Problem(mesh, problemInput, caseName="FRF")
    # Set path
    pFRF.setPath( gcodeFile )

    pFRF.setDt( 0.5*adimR )

    while not(pFRF.mhs.path.isOver(pFRF.time)):
        pFRF.iterate()
        pFRF.writepos()

def runCoupled():
    meshFixed  = getMesh()
    pFixed   = mhs.Problem(meshFixed, problemInput, caseName="fixed")
    myDriver = CustomStepper( pFixed, maxAdimtDt=adimR_domain-2, elementSize=fineElSize, adimMinRadius=3, threshold=0.015 )

    while not(pFixed.mhs.path.isOver(pFixed.time)):
        myDriver.iterate()


if __name__=="__main__":

    isRunReference = ("--run-reference" in sys.argv)
    onlyRunReference = ("--only-reference" in sys.argv)
    if isRunReference or onlyRunReference:
        runReference()

    if not(onlyRunReference):
        runCoupled()
