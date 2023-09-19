import MovingHeatSource as mhs
from CustomStepper import MyAdaptiveStepper, DriverReference
import thinWall
from scanningPath import writeGcode
import sys
import re
from MyLogger import MyLogger
import pickle

nLayers=None
problemInput = mhs.readInput("input.yaml")
gcodeFile = "Path.gcode"
layerThickness = problemInput["layerThickness"]
radius = problemInput["radius"]
maxTsteps = None

def deactivateBelowSurface(p,
                           surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeels = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < surfaceZ):
            activeels.append( ielem )
    substrateels = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeels )
    p.domain.setActivation( substrateels )

def runReference():
    mesh = thinWall.generateMesh()
    pReference = mhs.Problem( mesh, problemInput, caseName="reference" )
    deactivateBelowSurface( pReference ) 

    driver = DriverReference( pReference )

    logger = MyLogger()

    while not(driver.problem.mhs.path.isOver(pReference.time)):
        logger.iterate( driver )
        if maxTsteps:
            if pReference.iter >= maxTsteps:
                break

    with open("reference.log", "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)

def runCoupled():
    mesh = thinWall.generateMesh()
    pFixed = mhs.Problem( mesh, problemInput, caseName="fixed" )

    deactivateBelowSurface( pFixed ) 
    
    elementSize = thinWall.fineElementSize
    driver = MyAdaptiveStepper( pFixed, factor=2, maxAdimtDt=2,
                 threshold= 0.3, elementSize=elementSize, isCoupled=True, adimMinRadius=1.5 )

    logger = MyLogger()
    while not(pFixed.mhs.path.isOver(driver.getTime())):
        logger.iterate( driver )
        if maxTsteps:
            if driver.pFixed.iter >= maxTsteps:
                break

    with open("coupled.log", "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)

if __name__=="__main__":
    for arg in sys.argv:
        matchTsteps = re.search( r"--max-timesteps=(\d+)", arg )
        if matchTsteps:
            maxTsteps = int( matchTsteps.group(1) )
        matchLayers = re.search( r"--layers=(\d+)", arg )
        if matchLayers:
            nLayers = int( matchLayers.group(1) )

    writeGcode( gcodeFile, nLayers=nLayers )

    if ("--run-reference" in sys.argv) or ("--only-reference" in sys.argv):
        runReference()
    if not("--only-reference") in sys.argv:
        runCoupled()
