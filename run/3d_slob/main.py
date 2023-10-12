import MovingHeatSource as mhs
import numpy as np
import sys
from meshing import getMeshPhysical, fineElSize
from CustomStepper import CustomStepper, DriverReference, DriverAnalytical
from scanningPath import writeGcode
from MyLogger import MyLogger
import pickle

inputFile = "input.yaml"
problemInput = mhs.readInput( inputFile )
gcodeFile = problemInput["path"]
tol = 1e-7
maxIter = np.inf

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

def runReference():
    mesh = getMeshPhysical()
    pReference = mhs.Problem( mesh, problemInput, caseName="reference" )
    deactivateBelowSurface( pReference ) 

    driver = DriverReference( pReference )

    logger = MyLogger()

    iteration = 0
    while not(driver.problem.mhs.path.isOver(pReference.time)) and (iteration < maxIter):
        logger.iterate( driver )
        iteration += 1

    with open("reference.log", "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)

def runCoupled():
    adimR_tstep = 2
    fixedProblemInput = dict( problemInput )

    # Mesh
    meshFixed = getMeshPhysical()

    # mhs.Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput]:
        input["dt"] = dt

    pFixed         = mhs.Problem(meshFixed, fixedProblemInput, caseName="new")

    deactivateBelowSurface( pFixed )

    driver = CustomStepper( pFixed, maxAdimtDt=3, elementSize=fineElSize/2, threshold=0.2, adimMinRadius=1.25, adimZRadius=1.0 )
    
    logger = MyLogger()
    iteration = 0
    while (not(driver.pFixed.mhs.path.isOver( driver.getTime() ) ) and (iteration < maxIter)):
        logger.iterate( driver )
        iteration += 1

    with open("coupled.log", "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)

def runAnalytical():
    mesh = getMeshPhysical()
    p = mhs.Problem( mesh, problemInput, caseName="analytical" )
    deactivateBelowSurface( p ) 

    driver = DriverAnalytical( p )

    iteration = 0
    while not(driver.problem.mhs.path.isOver(p.time)) and (iteration < maxIter):
        driver.iterate()
        iteration += 1

if __name__=="__main__":
    isRunCoupled = False
    isRunReference = False
    isRunAnalytical = False

    isRunCoupled = ("--run-coupled" in sys.argv)
    isRunReference = ("--run-reference" in sys.argv)
    isRunAnalytical = ("--run-analytical" in sys.argv)

    nLayers = None
    for arg in sys.argv:
        if "--layers" in arg:
            nLayers = int( arg.split("=")[-1] )
        if "--maxIter" in arg:
            maxIter = int( arg.split("=")[-1] )
    writeGcode( nLayers=nLayers )

    if isRunReference:
        runReference()
    if isRunCoupled:
        runCoupled()
    if isRunAnalytical:
        runAnalytical()
