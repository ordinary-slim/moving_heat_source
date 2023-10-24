import MovingHeatSource as mhs
import numpy as np
import sys
from meshing import getMeshPhysical, fineElSize
from CustomStepper import CustomStepper, DriverReference
from scanningPath import writeGcode
from MyLogger import MyLogger
import pickle

inputFile = "input.yaml"
problemInput = mhs.readInput( inputFile )
gcodeFile = problemInput["path"]
tol = 1e-7
fineElSizeMoving = fineElSize/problemInput["fineElFactorMovingSubdomain"]

# read input
problemInput = mhs.readInput( inputFile )


def deactivateBelowSurface(p, surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < surfaceZ):
            activeEls.append( ielem )
    substrateEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( substrateEls )

def setAdimR( adimR, input ):
    r = input["radius"]
    speed  = np.linalg.norm( input["HeatSourceSpeed"] )
    return (adimR * r / speed)

def runReference(caseName="reference"):
    mesh = getMeshPhysical()
    pReference = mhs.Problem( mesh, problemInput, caseName=caseName)
    deactivateBelowSurface( pReference ) 

    driver = DriverReference( pReference )

    logger = MyLogger()

    while not(driver.problem.mhs.path.isOver(pReference.time)):
        logger.iterate( driver )

    with open("{}.log".format(caseName), "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)

def runCoupled(caseName="fixed"):
    adimR_tstep = 2
    fixedProblemInput = dict( problemInput )

    # Mesh
    meshFixed = getMeshPhysical()

    # mhs.Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput]:
        input["dt"] = dt

    pFixed = mhs.Problem(meshFixed, fixedProblemInput, caseName=caseName)

    deactivateBelowSurface( pFixed )

    driver = CustomStepper( pFixed,
                           adimFineDt=0.5 / problemInput["fineTStepFactor"],
                           maxAdimtDt=1,
                           elementSize=fineElSizeMoving,
                           threshold=0.3,
                           adimMinRadius=2,
                           #adimZRadius=1.0,
                           )
    
    logger = MyLogger()
    while not(driver.pFixed.mhs.path.isOver( driver.getTime() ) ) :
        logger.iterate( driver )

    with open("{}.log".format(caseName), "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)


if __name__=="__main__":
    isRunReference = ("--run-reference" in sys.argv)
    isRunCoupled = ("--run-coupled" in sys.argv)
    nLayers = None
    for arg in sys.argv:
        if "--layers" in arg:
            nLayers = int( arg.split("=")[-1] )
    caseName = "case"
    for arg in sys.argv:
        if "--case-name" in arg:
            caseName = str( arg.split("=")[-1] )
    writeGcode( nLayers=nLayers )
    if isRunReference:
        runReference(caseName=caseName)
    if isRunCoupled:
        runCoupled(caseName=caseName)
