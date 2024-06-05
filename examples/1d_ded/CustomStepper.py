import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper, TrackType
import numpy as np
import meshzoo
import pdb
import line_profiler

#TODO: Store hatch and not adim + nonAdim
#TODO: Build pMoving here in order to reduce risk of mistake

class CustomStepper(AdaptiveStepper):

    valuesOfMetric = [1e9]
    lastVals4mean = 2
    maximumTs = [np.inf]*(lastVals4mean+1)

    def onNewTrackOperations(self):
        self.rotateSubdomain()
        self.isCoupled = True
        speed = self.nextTrack.getSpeed()
        self.pMoving.setAdvectionSpeed( -speed )
        self.pMoving.domain.setSpeed( speed )
        self.pMoving.mhs.setPower( self.nextTrack.power )

    def increaseDt( self ):
        self.adimDt += self.adimFineDt
        for _ in range( self.lastVals4mean ):
            self.maximumTs.append( np.inf )

    def computeSizeSubdomain( self, adimDt = None ):
        if adimDt is None:
            adimDt = self.adimDt
        return adimDt + 4

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
        return (((self.valuesOfMetric[-1] < self.threshold) or (self.relChangeMetric < 0.02)) \
                and (self.maxTChange < 0.1))

    def update(self):
        '''
        Called at end of iteration
        Computes size of subdomain and dt
        '''
        # TODO: fix this for when track changes from t to tnp1
        
        time = self.pFixed.time
        self.onNewTrack = False

        # If path is over
        if (self.pFixed.mhs.path.isOver(time)):
            return
        tUnit = self.pFixed.mhs.radius / self.pFixed.mhs.currentTrack.speed

        dt2trackEnd = self.pFixed.mhs.currentTrack.endTime - time
        adimMaxDt = self.adimMaxDt

        if ( dt2trackEnd < 1e-7) or ( time== 0.0 ):#new track
            self.onNewTrack = True
            self.adimDt = self.adimFineDt

            if time != 0.0:
                self.nextTrack = self.pFixed.mhs.getNextTrack()

            # Make sure we don't skip over new track
            adimMaxDt = min( (self.nextTrack.endTime - time)/ tUnit, adimMaxDt )
        else:
            adimDt2TrackEnd = dt2trackEnd/tUnit
            maxDt2TrackEnd = adimDt2TrackEnd
            #if (maxDt2TrackEnd > self.adimMinRadius + 1e-7):
                #maxDt2TrackEnd -= self.adimMinRadius
            #else:
                #maxDt2TrackEnd = min( self.adimFineDt, adimDt2TrackEnd )

            adimMaxDt = min( maxDt2TrackEnd, self.adimMaxDt )
            self.computeSteadinessMetric(verbose=True)
            if (self.checkSteadinessCriterion()):
                self.increaseDt()

        # Cap dt
        self.adimDt = min( adimMaxDt, self.adimDt )
        self.setSizeSubdomain()
        self.dt = tUnit * self.adimDt


    def iterate( self ):
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.setActivation(self.physicalDomain)

        self.setDt()
        if self.onNewTrack:
            self.onNewTrackOperations()

        self.shapeSubdomain()

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#tn intersect

        # Motion, other operations
        self.pMoving.preIterate(canPreassemble=False)
        self.pFixed.preIterate(canPreassemble=False)

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#physical domain intersect

        if (self.hasPrinter and self.pFixed.input["print"]):
            self.deposit()

        if self.isCoupled:
            self.pFixed.substractExternal(self.pMoving, updateGamma=False)
            self.pFixed.updateInterface( self.pMoving )
            self.pMoving.updateInterface( self.pFixed )
            # Set interface boundary conditions
            self.pFixed.setGamma2Dirichlet()
            self.pMoving.setGamma2Neumann()
        self.setBCs()

        self.pMoving.preAssemble(allocateLs=False)

        if self.isCoupled:
            # Pre-assembly, updating free dofs
            self.pFixed.preAssemble(allocateLs=False)
        
            ls = mhs.LinearSystem.Create( self.pMoving, self.pFixed )
            # Assembly
            self.pMoving.assemble( self.pFixed )
            self.pFixed.assemble( self.pMoving )
            # Build ls
            ls.assemble()
            # Solve ls
            ls.solve()
            # Recover solution
            self.pFixed.gather()
            self.pMoving.gather()

            self.pFixed.unknown.interpolateInactive( self.pMoving.unknown, ignoreOutside = False )
            self.pMoving.unknown.interpolateInactive( self.pFixed.unknown, ignoreOutside = False )
            # Post iteration
            self.pFixed.postIterate()
            self.pMoving.postIterate()
        else:
            self.pFixed.clearGamma()
            self.pFixed.preAssemble(allocateLs=True)
            self.pFixed.iterate()
            try:
                self.pMoving.unknown.interpolate( self.pFixed.unknown, ignoreOutside = False )
            except ValueError:
                self.writepos()
                raise
            self.pMoving.postIterate()

        # Compute next dt && domain size
        try:
            self.update()
        except ZeroDivisionError:
            self.writepos()
            raise
        # Write vtk files
        self.writepos()

    def writepos( self ):
        self.pFixed.writepos(
            functions={
                "Tn":(self.pFixed.previousValues[0]),
                "DTdt":((self.pFixed.unknown - self.pFixed.previousValues[0]) / self.pFixed.dt),
                },
            nodeMeshTags={
                "gammaNodes":self.pFixed.gammaNodes,
                "forcedDofs":self.pFixed.forcedDofs,
                },
            cellMeshTags={
                "physicalDomain":self.physicalDomain,
                },
                          )
        self.pMoving.writepos(
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )




class DriverReference:
    def __init__(self, problem, adimDt=0.5):
        self.problem = problem
        if "path" in self.problem.input:
            self.problem.setPath( self.problem.input["path"] )
        self.dt2trackEnd = None
        self.printer = None
        self.adimDt  = adimDt
        if "printer" in self.problem.input:
            self.printer = mhs.Printer( self.problem,
                                     self.problem.input["printer"]["width"],
                                     self.problem.input["printer"]["height"],
                                     self.problem.input["printer"]["depth"]
                                     )
        self.isCoupled = False
        self.tol = 1e-7
        self.nextTrack = None

    def setDtFromAdimR( self, adimR, maxDt=None):
        r =self.problem.input["radius"]
        speed = np.linalg.norm(self.problem.input["HeatSourceSpeed"] )
        dt =  adimR * r / speed 
        if maxDt:
            dt = min(dt, maxDt)
        self.problem.setDt( dt )

    def solve(self):
        self.problem.iterate()
    def writepos(self):
        self.problem.writepos(
                functions={
                    "Tn":(self.problem.previousValues[0]),
                    "DTdt":((self.problem.unknown - self.problem.previousValues[0]) / self.problem.dt),
                    }
                )

    def iterate( self ):
        self.setDtFromAdimR( self.adimDt, self.dt2trackEnd )
        tnp1 = self.problem.time + self.problem.dt
        track = self.problem.mhs.path.interpolateTrack( tnp1 )
        if ((track.type==TrackType.printing) and (self.printer is not None) and (self.problem.input["print"])):
            self.printer.deposit( self.problem.mhs.position,
                                self.problem.mhs.path.interpolatePosition(tnp1)
                               )
        #self.problem.setConvection( resetBcs = True )
        self.solve()
        self.dt2trackEnd = self.problem.mhs.currentTrack.endTime - self.problem.time
        if self.dt2trackEnd < self.tol:
            self.nextTrack = self.problem.mhs.getNextTrack()
            self.dt2trackEnd = self.nextTrack.endTime - self.problem.time
        self.writepos()
