import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper, TrackType
import numpy as np
import meshzoo
import pdb

class CustomStepper(AdaptiveStepper):

    def setBCs( self ):
        self.pFixed.setConvection( resetBcs = False )
        if self.isCoupled:
            self.pMoving.setConvection( resetBcs = False )

    #def computeSizeSubdomain( self, adimDt = None ):
        #if adimDt is None:
            #adimDt = self.adimDt
        #return adimDt + 4

    def deposit( self ):
        '''
        Do deposition for interval [tn, tn1] at tn+1
        '''
        if (self.pFixed.mhs.currentTrack.type == TrackType.printing):
            origin = self.pFixed.mhs.path.interpolatePosition( self.getTime() - self.dt )
            destination = self.pFixed.mhs.position

            self.printerFixed.melt( origin, destination,
                                    idxTargetSet = 0,
                                )
            self.printerMoving.melt( origin, destination,
                                    idxTargetSet = 0,
                                )

    def writepos( self ):
        self.pFixed.writepos(
            nodeMeshTags={
                "gammaNodes":self.pFixed.gammaNodes,
                "forcedDofs":self.pFixed.forcedDofs,
                },
            cellMeshTags={
                "material":self.pFixed.domain.materialTag,
                "physicalDomain":self.physicalDomain,
                },
                          )
        self.pMoving.writepos(
            cellMeshTags={
                "material":self.pMoving.domain.materialTag,
                },
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )

    def onNewTrackOperations(self):
        self.isCoupled = False
        speed = self.nextTrack.getSpeed()
        self.pMoving.setAdvectionSpeed( -speed )
        self.pMoving.domain.setSpeed( speed )
        self.pMoving.mhs.setPower( self.nextTrack.power )

        if self.nextTrack.type == TrackType.printing:
            self.rotateSubdomain()
        if (self.nextTrack.type == TrackType.recoating):
            activatePowderLayer3d( self.pFixed, self.printerFixed )
            activatePowderLayer3d( self.pMoving, self.printerMoving )

    def iterate( self ):
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.setActivation(self.physicalDomain)

        self.setDt()
        self.setCoupling()

        # RECOAT
        if self.isNewLayer():
            self.onNewLayerOperations()
        # SET ADVECTION IN MOVING SUBDOMAIN
        if self.onNewTrack:
            self.onNewTrackOperations()


        self.shapeSubdomain()

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#tn intersect

        # Motion, other operations
        self.pMoving.preIterate(canPreassemble=False)
        self.pFixed.preIterate(canPreassemble=False)

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#physical domain intersect

        self.pMoving.domain.setMaterialSets(
                self.pMoving.domain.projectCellTag(
                    self.pFixed.domain.materialTag,
                    self.pFixed.domain,
                    )
                )

        if self.hasPrinter:
            self.deposit()

        # At this point, activation in fixed problem corresponds
        # to physical domain
        self.physicalDomain = mhs.MeshTag( self.pFixed.domain.activeElements )

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

        self.prettyPrintIteration()

        # Compute next dt && domain size
        if self.isAdaptive:
            try:
                self.update()
            except ZeroDivisionError:
                self.writepos()
                raise
        # Write vtk files
        self.writepos()



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
        self.setDtFromAdimR( 0.5, self.dt2trackEnd )
        self.tnp = float(self.problem.time)
        self.tnp1 = self.problem.time + self.problem.dt
        self.track = self.problem.mhs.path.interpolateTrack( self.tnp1 )
        self.problem.setConvection( resetBcs = True )

        if not(self.problem.hasPreIterated):
            self.problem.preIterate(canPreassemble=False)
        # DEPOSIT NEW LAYER
        if (self.onNewTrack) and (self.problem.mhs.currentTrack.type == TrackType.recoating):
            activatePowderLayer3d( self.problem, self.printer )
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

def deactivateBelowSurface3d(p, surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeEls = []
    powderEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < surfaceZ):
            activeEls.append( ielem )
        else:
            powderEls.append( ielem )
    # DEACTIVE
    substrateEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( substrateEls )
    # SET 2 POWDER
    powderEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, powderEls )
    p.domain.setMaterialSets( powderEls )

def activatePowderLayer3d( p, printer, powderSet = 1 ):
    '''
    If new layer of problem, activate powder layer
    '''
    currentZ = p.mhs.position[2]
    thresholdZ = currentZ + p.input["printer"]["height"]

    nels = p.domain.mesh.nels
    newPowderEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < thresholdZ) and \
        not(p.domain.activeElements[ielem]):
            newPowderEls.append( ielem )
    p.domain.activeElements.setIndices( newPowderEls, 1 )
    p.domain.setActivation( p.domain.activeElements )
    printer.setDepositionTemperature()
    p.domain.materialTag.setIndices( newPowderEls, powderSet )
