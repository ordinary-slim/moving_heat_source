import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper
import numpy as np
import meshzoo
import pdb

#TODO: Store hatch and not adim + nonAdim
#TODO: Build pMoving here in order to reduce risk of mistake

class CustomStepper(AdaptiveStepper):

    valuesOfMetric = [np.inf]
    lastVals4mean = 2
    maximumTs = [np.inf]*lastVals4mean

    def setBCs( self ):
        self.pFixed.setConvection( resetBcs = False )
        if self.isCoupled:
            self.pMoving.setConvection( resetBcs = False )

    def computeSteadinessMetric( self, verbose=True ):
        # TODO: Change this
        delta = self.pMoving.unknown - self.pMoving.previousValues[0]
        self.metric = delta.getL2Norm() / self.pMoving.unknown.getL2Norm()
        self.maximumTs.append( max(self.pMoving.unknown.values) )
        movingAv = np.mean( np.array( self.maximumTs[(-1-self.lastVals4mean):-1] ) )
        self.maxTChange = abs(self.maximumTs[-1] - movingAv ) / movingAv
        self.valuesOfMetric.append( self.metric )
        self.relChangeMetric = abs( (self.valuesOfMetric[-1] - self.valuesOfMetric[-2]) / self.valuesOfMetric[-1] )
        if verbose:
            print( " || fn+1 - fn || / || fn+1 || = {} ".format( self.valuesOfMetric[-1] ) )
            print( " threshold = {}".format( self.threshold ) )
            print( "Relative change metric = {}%".format( 100*self.relChangeMetric ) )
            try:
                print( "Previous max T = {}, current max T = {}".format( self.maximumTs[-2], self.maximumTs[-1] ) )
            except IndexError:
                pass
            print( "Relative change max T change over last iters = {}%".format( 100*self.maxTChange ) )

    def checkSteadinessCriterion(self):
        return (((self.valuesOfMetric[-1] < self.threshold) or (self.relChangeMetric < 0.01)) \
                and (self.maxTChange < 0.05))

    def increaseDt( self ):
        self.adimDt = min( self.factor * self.adimDt, self.adimDt + self.adimFineDt )
        self.maximumTs.extend( [np.inf]*self.lastVals4mean )

    def computeSizeSubdomain( self, adimDt = None ):
        if adimDt is None:
            adimDt = self.adimDt
        return 4

    def onNewTrackOperations(self):
        self.rotateSubdomain()
        self.isCoupled = True
        speed = self.nextTrack.getSpeed()
        self.pMoving.setAdvectionSpeed( -speed )
        self.pMoving.domain.setSpeed( speed )
        self.pMoving.mhs.setPower( self.nextTrack.power )

    def deposit( self ):
        '''
        Do deposition for interval [tn, tn1] at tn+1
        '''
        if (self.pFixed.mhs.currentTrack.hasDeposition):
            origin = self.pFixed.mhs.path.interpolatePosition( self.getTime() - self.dt )
            destination = self.pFixed.mhs.position

            self.printerFixed.deposit( origin, destination,
                                self.physicalDomain,
                                )
            self.printerMoving.deposit( origin, destination,
                                self.pMoving.domain.activeElements,
                                modifyValues=True,
                                )

    def writepos( self ):
        self.pFixed.writepos(
            nodeMeshTags={
                "gammaNodes":self.pFixed.gammaNodes,
                "forcedDofs":self.pFixed.forcedDofs,
                },
            cellMeshTags={
                "physicalDomain":self.physicalDomain,
                },
                          )
        self.pMoving.writepos(
            functions={"prevVal":self.pMoving.previousValues[0]},
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )



class DriverReference:
    def __init__(self, problem, adimFineDt = 0.5):
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
        self.adimDt = adimFineDt / self.problem.input["fineTStepFactor"]
        self.isCoupled = False
        self.tol = 1e-7
        self.nextTrack = None

    def setDtFromAdimR( self, adimR, maxDt=None):
        r =self.problem.input["radius"]
        speed = max( np.linalg.norm(self.problem.input["HeatSourceSpeed"] ), np.linalg.norm(self.problem.input["advectionSpeed"] ) )
        dt =  adimR * r / speed 
        if maxDt:
            dt = min(dt, maxDt)
        self.problem.setDt( dt )

    def iterate( self ):
        self.setDtFromAdimR( self.adimDt, self.dt2trackEnd )
        tnp1 = self.problem.time + self.problem.dt
        track = self.problem.mhs.path.interpolateTrack( tnp1 )
        if (track.hasDeposition) and (self.printer is not None):
            self.printer.deposit( self.problem.mhs.position,
                                self.problem.mhs.path.interpolatePosition(tnp1)
                               )
        self.problem.setConvection( resetBcs = True )
        self.problem.iterate()
        self.dt2trackEnd = self.problem.mhs.currentTrack.endTime - self.problem.time
        if self.dt2trackEnd < self.tol:
            self.nextTrack = self.problem.mhs.getNextTrack()
            self.dt2trackEnd = self.nextTrack.endTime - self.problem.time
        self.problem.writepos()
