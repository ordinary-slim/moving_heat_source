import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
import numpy as np
import meshzoo
from CustomStepper import CustomStepper
import argparse

inputFile = "input.yaml"
tol = 1e-7
adimR_domain = 10 
problemInput = mhs.readInput( inputFile )
adimR = problemInput["radius"] / np.linalg.norm(problemInput["HeatSourceSpeed"])
boxPhys = [-7, 7, -7, 7]
gcodeFile = problemInput["path"]
radiusHs = problemInput["radius"]
fineElSize = radiusHs/2/problemInput["elSizeFactor"]

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

def writeGcode():
    gcodeLines = []
    '''
    lenX = boxPhys[1] - boxPhys[0]
    lenY = boxPhys[3] - boxPhys[2]
    lenPathX = lenX * 0.75
    lenPathY = lenY * 0.75
    '''
    x0 = -5
    x1 = +5
    y0 = -5
    y1 = +5
    speed = np.array(problemInput["HeatSourceSpeed"])
    V = np.linalg.norm( speed )
    gcodeLines.append(
            "G0 F{} X{} Y{} Z0".format( V, x0, y0 ))
    gcodeLines.append(
            "G1 X{} E0.1".format( x1 ))
    gcodeLines.append(
            "G1 Y{} E0.2".format( y1 ))
    gcodeLines.append(
            "G1 X{} E0.3".format( x0 ))
    gcodeLines.append(
            "G1 Y{} E0.4".format( y0 ))
    with open(problemInput["path"], 'w') as gf:
        gf.writelines([l+"\n" for l in gcodeLines])

def runReference(caseName="frf"):
    mesh = getMesh(elementSize=fineElSize)
    pFRF     = mhs.Problem(mesh, problemInput, caseName=caseName)
    # Set path
    pFRF.setPath( gcodeFile )

    pFRF.setDt( 0.5*adimR/problemInput["frfTstepFactor"] )

    while not(pFRF.mhs.path.isOver(pFRF.time)):
        pFRF.iterate()
        pFRF.writepos()

def runCoupled(caseName="coupled"):
    meshFixed = getMesh(elementSize=fineElSize)
    pFixed   = mhs.Problem(meshFixed, problemInput, caseName=caseName)
    myDriver = CustomStepper( pFixed,
                              elementSize=fineElSize,
                              adimMinRadius=4,
                              slowDown=False,
                              threshold=0.02,
                              alwaysCoupled=False,
                              adimFineDt=0.5,
                             )
                              #adimPosZLen=0.5,
                              #adimNegZLen=2.0,
                              #adimSideRadius=2.0,)

    while not(pFixed.mhs.path.isOver(pFixed.time)):
        myDriver.iterate()


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run-reference', action='store_true')
    parser.add_argument('-c', '--run-coupled', action='store_true')
    parser.add_argument('--layers', default=-1, type=int)
    parser.add_argument('--case-name', default='case')
    args = parser.parse_args()
    writeGcode()
    if args.run_reference:
        runReference(caseName=args.case_name)
    if args.run_coupled:
        runCoupled(caseName=args.case_name)
