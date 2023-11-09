import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import LpbfAdaptiveStepper, TrackType, activatePowderLayer, deactivateBelowSurface
import numpy as np
import meshzoo
import pdb

class CustomStepper(LpbfAdaptiveStepper):
    def setBCs( self ):
        self.pFixed.setConvection( resetBcs = False )
        if self.isCoupled:
            self.pMoving.setConvection( resetBcs = False )

    def increaseDt( self ):
        if self.adimDt <= 1.5 + 1e-7:
            self.adimDt += self.adimFineDt
        else:
            self.adimDt += 2*self.adimFineDt

    def computeSizeSubdomain( self, adimDt = None ):
        if adimDt is None:
            adimDt = self.adimDt
        return min(adimDt + 4*self.adimFineDt, adimDt * 2 )

    def getIsPrinting( self ):
        return (self.pFixed.mhs.currentTrack.type == TrackType.printing)

class DriverReference:
    def __init__(self, problem):
        self.problem = problem
        if "path" in self.problem.input:
            self.problem.setPath( self.problem.input["path"] )
        self.dt2trackEnd = None
        self.printer = None
        if "printer" in self.problem.input:
            self.printer = mhs.Printer( self.problem,
                                     self.problem.input["printer"]["width"],
                                     self.problem.input["printer"]["height"],
                                     self.problem.input["printer"]["depth"]
                                     )
        self.isCoupled = False
        self.tol = 1e-7
        self.nextTrack = self.problem.mhs.currentTrack
        self.onNewTrack = True
        self.adimDt = 0.5 / self.problem.input["fineTStepFactor"]

    def setDtFromAdimR( self, adimR, maxDt=None):
        if (self.nextTrack.type == TrackType.printing) or (self.nextTrack.type == TrackType.cooling):
            r =self.problem.input["radius"]
            speed = max( np.linalg.norm(self.problem.input["HeatSourceSpeed"] ), np.linalg.norm(self.problem.input["advectionSpeed"] ) )
            dt =  adimR * r / speed 
        else:
            dt = self.problem.input["cooling_dt"]
        if maxDt:
            dt = min(dt, maxDt)
        self.problem.setDt( dt )

    def print(self):
        self.printer.melt( self.problem.mhs.path.interpolatePosition(self.tnp),
                           self.problem.mhs.position,
                           idxTargetSet = 0,
                           )

    def iterate( self ):
        self.setDtFromAdimR( self.adimDt, self.dt2trackEnd )
        self.tnp = float(self.problem.time)
        self.tnp1 = self.problem.time + self.problem.dt
        self.track = self.problem.mhs.path.interpolateTrack( self.tnp1 )
        self.problem.setConvection( resetBcs = True )

        if not(self.problem.hasPreIterated):
            self.problem.preIterate(canPreassemble=False)
        # DEPOSIT NEW LAYER
        if (self.onNewTrack) and (self.problem.mhs.currentTrack.type == TrackType.recoating):
            activatePowderLayer( self.problem, self.printer )
        # PRINT
        if (self.track.type == TrackType.printing) and (self.printer is not None):
            self.print()

        self.problem.preAssemble()
        self.problem.assemble()
        self.problem.ls.solve()
        self.problem.gather()
        self.problem.postIterate()
        print( "{} iter# {}, time={}".format(
            self.problem.caseName,
            self.problem.iter,
            self.problem.time) )

        self.onNewTrack = False
        self.dt2trackEnd = self.problem.mhs.currentTrack.endTime - self.problem.time
        if self.dt2trackEnd < self.tol:
            self.onNewTrack = True
            self.nextTrack = self.problem.mhs.getNextTrack()
            self.dt2trackEnd = self.nextTrack.endTime - self.problem.time

        self.writepos()

    def writepos(self):
        self.problem.writepos(
                cellMeshTags={
                    "material":self.problem.domain.materialTag,
                    },
                )

    def getIsPrinting( self ):
        return (self.problem.mhs.currentTrack.type == TrackType.printing)

