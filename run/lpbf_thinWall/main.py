import MovingHeatSource as mhs
from customStepper import MyAdaptiveStepper
import thinWall
from scanningPath import writeGcode
import numpy as np
import sys
import re
import pdb

nLayers=None
problemInput = mhs.readInput("input.yaml")
gcodeFile = "Path.gcode"
layerThickness = problemInput["layerThickness"]
radius = problemInput["radius"]
maxTsteps = None

def setDtFromAdimR( problem, adimR, maxDt=None):
    r = problem.input["radius"]
    speed = max( np.linalg.norm( problem.input["HeatSourceSpeed"] ), np.linalg.norm( problem.input["advectionSpeed"] ) )
    dt =  adimR * r / speed 
    if maxDt:
        dt = min(dt, maxDt)
    problem.setDt( dt )

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
    for p in [pReference]:
        p.setPath( gcodeFile )
        deactivateBelowSurface( p ) 
    dt2trackEnd = None
    printerRef = None
    if "printer" in problemInput:
        printerRef = mhs.Printer( pReference,
                                 problemInput["printer"]["width"],
                                 problemInput["printer"]["height"],
                                 problemInput["printer"]["depth"]
                                 )
    while not(pReference.mhs.path.isOver(pReference.time)):
        setDtFromAdimR( pReference, 0.5, dt2trackEnd )
        tnp1 = pReference.time + pReference.dt
        track = pReference.mhs.path.interpolateTrack( tnp1 )
        if (track.hasDeposition) and (printerRef is not None):
            printerRef.deposit( pReference.mhs.position,
                                pReference.mhs.path.interpolatePosition(tnp1)
                               )
        pReference.iterate()
        dt2trackEnd = pReference.mhs.currentTrack.endTime - pReference.time
        pReference.writepos()
        if maxTsteps:
            if pReference.iter >= maxTsteps:
                break

def runCoupled():
    mesh = thinWall.generateMesh()
    pFixed = mhs.Problem( mesh, problemInput, caseName="fixed" )

    deactivateBelowSurface( pFixed ) 
    
    elementSize = thinWall.fineElementSize
    driver = MyAdaptiveStepper( pFixed, factor=2, adimMaxSubdomainSize=10,
                 threshold= 0.3, elementSize=elementSize, isCoupled=True )

    while not(pFixed.mhs.path.isOver(driver.getTime())):
        driver.iterate()
        if maxTsteps:
            if driver.pFixed.iter >= maxTsteps:
                break

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
