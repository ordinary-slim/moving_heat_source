import numpy as np
import gmsh
import MovingHeatSource as mhs
from MovingHeatSource.gcode import TrackType
import pdb

#TODO: Store hatch and not adim + nonAdim
#TODO: Build pMoving here in order to reduce risk of mistake

class AdaptiveStepper:

    tol = 1e-7

    def __init__(self,
                 pFixed,
                 isAdaptive=True,
                 factor=2,
                 maxAdimtDt=8,
                 threshold= 0.01,
                 elementSize=[0.25]*3,
                 shift=None,
                 adimFineDt = 0.5,
                 adimMinRadius=2,
                 adimPosZLen=None,
                 adimNegZLen=None,
                 adimSideRadius=None,
                 alwaysCoupled=False,
                 slowDown=True,
                 slowAdimDt=None,
                 ):

        self.pFixed = pFixed
        self.isAdaptive = isAdaptive
        self.adimFineDt = adimFineDt
        # Set up laser path
        if "path" in self.pFixed.input:
            self.pFixed.setPath( self.pFixed.input["path"] )
        else:
            raise Exception("AdapativeStepper requires gcode.")
        self.nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time )

        self.alwaysCoupled = alwaysCoupled
        self.isCoupled = alwaysCoupled
        self.slowDown = slowDown
        self.adimMaxDt = maxAdimtDt
        self.adimMaxSubdomainSize = self.computeSizeSubdomain(self.adimMaxDt)
        self.adimMinRadius = adimMinRadius
        self.adimNegZLen = adimMinRadius
        self.adimSideRadius = adimMinRadius
        if not(adimNegZLen is None):
            self.adimNegZLen = adimNegZLen
        self.adimPosZLen = self.adimNegZLen
        if not(adimPosZLen is None):
            self.adimPosZLen = adimPosZLen
        if not(adimSideRadius is None):
            self.adimSideRadius = adimSideRadius
        self.slowAdimDt  = adimFineDt
        if not(slowAdimDt is None):
            self.slowAdimDt = slowAdimDt
        self.buildMovingProblem(elementSize=elementSize, shift=shift)
        # TODO: better initialization
        self.threshold = threshold
        self.factor = factor

        # Cooling
        self.cooling_dt = self.adimFineDt
        if "cooling_dt" in self.pFixed.input:
            self.cooling_dt = self.pFixed.input["cooling_dt"]
        self.traveledDistance = 0.0
        self.update()
        self.onNewTrack = True
        self.hasPrinter = False
        self.physicalDomain = mhs.MeshTag( self.pFixed.domain.activeElements )

        # Match heat source positions at t = 0.0
        self.pMoving.mhs.setPosition( self.pFixed.mhs.position )

        # Set up printer
        if "printer" in self.pFixed.input:
            self.setupPrinter( pFixed.input["printer"]["width"], pFixed.input["printer"]["height"], pFixed.input["printer"]["depth"] )

        # Initialize steadiness criterion
        self.metric = 1e99
        self.idxSolverUncoupledIter = 1
        if "idxSolverUncoupledIter" in self.pFixed.input:
            self.idxSolverUncoupledIter = self.pFixed.input["idxSolverUncoupledIter" ]
        self.idxSolverCoupledIter = self.pFixed.idxSolver
        if "idxSolverCoupledIter" in self.pFixed.input:
            self.idxSolverCoupledIter = self.pFixed.input["idxSolverCoupledIter" ]

    def buildMovingProblem(self, elementSize, shift):
        meshInputMoving = {}
        try:
            len(elementSize)
        except TypeError:
            elementSize = [elementSize]*3
        self.meshMoving = self.meshAroundHS(elementSizes=elementSize, shift=shift)
        self.currentOrientation = np.array([1.0, 0.0, 0.0])
        movingProblemInput = dict( self.pFixed.input )
        movingProblemInput["isStabilized"] = 1
        movingProblemInput["isAdvection"] = 1
        movingProblemInput["advectionSpeed"] = -self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["speedDomain"]      = self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["HeatSourceSpeed"] = np.zeros(3)
        movingProblemInput["initialPosition"] = self.pFixed.mhs.position
        self.pMoving  = mhs.Problem(self.meshMoving, movingProblemInput, caseName=self.pFixed.caseName + "_moving")
        # Quick-fix: Mesh is built oriented around X
        if (self.nextTrack.type == TrackType.printing):
            self.rotateSubdomain()

    def getTime(self):
        return self.pFixed.time

    def setupPrinter( self, mdwidth, mdheight, mddepth ):
        self.printerFixed = mhs.Printer( self.pFixed, mdwidth, mdheight, mddepth )
        self.printerMoving = mhs.Printer( self.pMoving, mdwidth, mdheight, mddepth )
        self.hasPrinter = True

    def deposit( self ):
        '''
        Do deposition for interval [tn, tn1] at tn+1
        '''
        if (self.pFixed.mhs.currentTrack.type == TrackType.printing):
            origin = self.pFixed.mhs.path.interpolatePosition( self.getTime() - self.dt )
            destination = self.pFixed.mhs.position

            self.printerFixed.deposit( origin, destination,
                                self.physicalDomain,
                                )
            self.printerMoving.deposit( origin, destination,
                                self.pMoving.domain.activeElements,
                                )

    def computeSteadinessMetric( self, verbose=True ):
        delta = self.pMoving.unknown - self.pMoving.previousValues[0]
        self.metric = delta.getL2Norm() / self.pMoving.unknown.getL2Norm()
        if verbose:
            print( "metric = {}".format( self.metric ) )

    def checkSteadinessCriterion( self ):
        return (self.metric < self.threshold)

    def increaseDt( self ):
        self.adimDt = min( self.factor * self.adimDt, self.adimDt + 2*self.adimFineDt )

    def computeSizeSubdomain( self, adimDt = None ):
        if adimDt is None:
            adimDt = self.adimDt
        return min(adimDt + 4*self.adimFineDt, adimDt * 2 )

    def setSizeSubdomain( self, adimDt = None ):
        self.adimSubdomainSize = self.computeSizeSubdomain( adimDt = adimDt )

    def update(self):
        '''
        Called at end of iteration
        Computes size of subdomain and dt
        '''
        # If path is over
        if (self.pFixed.mhs.path.isOver(self.getTime())):
            return

        self.onNewTrack = False

        self.dt2trackEnd = self.pFixed.mhs.currentTrack.endTime - self.getTime()

        # Check if next step starts new track
        if ( self.dt2trackEnd < 1e-7) or ( self.getTime()== 0.0 ):#new track
            self.onNewTrack = True
            if self.getTime() != 0.0:
                self.nextTrack = self.pFixed.mhs.getNextTrack()
            # Store next dt 2 track end
            self.dt2trackEnd = (self.nextTrack.endTime - self.getTime())

        if   self.nextTrack.type == TrackType.printing:
            self.printingUpdate()
        else:
            self.coolingUpdate()

    def printingUpdate( self ):
        speed = np.linalg.norm( self.nextTrack.getSpeed() )
        tUnit = self.pFixed.mhs.radius / speed
        self.maxAdimDtPrinting = self.adimMaxDt
        if (self.onNewTrack):
            self.adimDt = self.adimFineDt
            # Make sure we don't skip over new track
            self.maxAdimDtPrinting = min( self.dt2trackEnd / tUnit, self.maxAdimDtPrinting )
        else:
            adimDt2TrackEnd = self.dt2trackEnd/tUnit
            maxDt2TrackEnd = adimDt2TrackEnd
            if self.slowDown:
                if (maxDt2TrackEnd > self.adimMinRadius + 1e-7):
                    maxDt2TrackEnd -= self.adimMinRadius
                else:
                    maxDt2TrackEnd = min( self.slowAdimDt, adimDt2TrackEnd )

            self.maxAdimDtPrinting = min( maxDt2TrackEnd, self.adimMaxDt )
            self.computeSteadinessMetric(verbose=True)
            if (self.checkSteadinessCriterion()):
                self.isCoupled = True
                self.increaseDt()

        # Cap dt
        self.adimDt = min( self.maxAdimDtPrinting, self.adimDt )
        self.setSizeSubdomain()
        self.dt = tUnit * self.adimDt
        self.traveledDistance += self.adimDt

    def coolingUpdate(self):
        self.dt = min( self.cooling_dt, self.dt2trackEnd )
        self.setSizeSubdomain( adimDt=self.adimFineDt )
        self.traveledDistance = 0.0

    def rotateSubdomain( self ):
        '''
        Align mesh of moving subproblem with track
        Called at tn
        '''
        nextOrientation    = self.nextTrack.getSpeed() / np.linalg.norm( self.nextTrack.getSpeed() )
        center = self.pMoving.mhs.position
        angle = np.arccos( np.dot( nextOrientation, self.currentOrientation ) / np.linalg.norm( nextOrientation ) / np.linalg.norm( self.currentOrientation ) )
        if (angle > 1e-5):
            self.pMoving.domain.inPlaneRotate( center, angle )
            self.pMoving.unknown.interpolate( self.pFixed.unknown, ignoreOutside = True )
            self.currentOrientation = nextOrientation

    def shapeSubdomain( self ):
        '''
        At t^n, do things
        '''
        if not(self.nextTrack.type == TrackType.printing):
            return
        # OBB
        radius = self.pFixed.mhs.radius

        # compute front and sides
        sideRadius = self.adimSideRadius * radius
        adimBackRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize )
        backRadiusAabb = max( adimBackRadius, self.adimMinRadius ) * radius
        frontRadiusAabb = self.adimMinRadius*radius
        zRadius    = max(self.adimNegZLen, self.adimPosZLen) * radius
        xAxis      = self.nextTrack.getSpeed() / self.nextTrack.speed

        #backRadiusObb = max(backRadius - radius, 0.0)
        pback  = self.pMoving.mhs.position - backRadiusAabb*xAxis
        pfront = self.pMoving.mhs.position + frontRadiusAabb*xAxis
        paabb  = self.pMoving.mhs.position + (frontRadiusAabb - backRadiusAabb)*xAxis/2
        hLensAabb = np.array( [ (frontRadiusAabb + backRadiusAabb)/2.0,
                            sideRadius,
                            zRadius,] )
        aabb = mhs.MyAABB( paabb, hLensAabb, True )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( aabb )
        #collidingElsBackSphere = self.pMoving.domain.mesh.findCollidingElements( p0, self.adimMinRadius*radius )
        #collidingElsFrontSphere = self.pMoving.domain.mesh.findCollidingElements( self.pMoving.mhs.position, self.adimMinRadius*radius )
        #subdomainEls += collidingElsBackSphere
        #subdomainEls += collidingElsFrontSphere
        subdomain = mhs.MeshTag( self.pMoving.domain.mesh, self.pMoving.domain.mesh.dim, subdomainEls )
        self.pMoving.domain.intersect( subdomain )


    def setDt( self ):
        self.pMoving.setDt( self.dt )
        self.pFixed.setDt( self.dt )

    '''
    def setCoupling( self ):
        # Set coupling
        if ((not(self.nextTrack.type == TrackType.printing) or \
            ((self.adimDt <= self.adimFineDt+1e-7) and not(self.isCoupled))) \
            and not(self.alwaysCoupled)):
            self.isCoupled = False
        else:
            self.isCoupled = True
    '''

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
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )


    def onNewTrackOperations(self):
        if not(self.alwaysCoupled):
            self.isCoupled = False
        speed = self.nextTrack.getSpeed()
        self.pMoving.setAdvectionSpeed( -speed )
        self.pMoving.domain.setSpeed( speed )
        self.pMoving.mhs.setPower( self.nextTrack.power )

        if self.nextTrack.type == TrackType.printing:
            self.rotateSubdomain()
        if (self.nextTrack.type == TrackType.recoating):
            activatePowderLayer( self.pFixed, self.printerFixed )
            activatePowderLayer( self.pMoving, self.printerMoving )

    def onNewLayerOperations(self):
        pass

    def setBCs(self):
        pass

    def isNewLayer( self ):
        isNewLayer = (self.nextTrack.isNewY) and ((self.onNewTrack) and self.nextTrack.type == TrackType.printing)
        return isNewLayer

    def prettyPrintIteration( self ):
        if self.isCoupled:
            print( "{} iter# {}, time={}, travelled distance = {}R".format(
                self.pFixed.caseName,
                self.pFixed.iter,
                self.pFixed.time,
                self.traveledDistance) )
        # Called at tn
        if self.pFixed.mhs.currentTrack.type == TrackType.printing :
            print(" dt = {}R, domainSize = {}R, is coupled = {}".format( self.adimDt,
                                                                        self.adimSubdomainSize,
                                                                        self.isCoupled ) )

    def iterate( self ):
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.setActivation(self.physicalDomain)

        self.setDt()
        #self.setCoupling()

        if self.isNewLayer():
            self.onNewLayerOperations()
        if self.onNewTrack:
            self.onNewTrackOperations()

        self.shapeSubdomain()

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#tn intersect

        # Motion, other operations
        self.pMoving.preIterate(canPreassemble=False)
        self.pFixed.preIterate(canPreassemble=False)

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#physical domain intersect

        if self.hasPrinter:
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
            ls.setSolver( self.idxSolverCoupledIter )
            ls.setInitialGuess( self.pMoving, self.pFixed )
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
            self.pFixed.idxSolver = self.idxSolverUncoupledIter
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

    def meshAroundHS( self, elementSizes, shift=None ):
        # Compute lengths of box
        radius = self.pFixed.mhs.radius
        minRadius = self.adimMinRadius  * radius
        YRadius   = self.adimSideRadius * radius
        negZRadius   = self.adimNegZLen * radius
        posZRadius   = self.adimPosZLen * radius
        meshCenter = np.array(self.pFixed.mhs.position)
        if not(shift is None):
            meshCenter += shift
        lengths = {}
        lengths["trailLength"] = self.adimMaxSubdomainSize * radius
        lengths["capotLength"] = min( lengths["trailLength"], minRadius )
        lengths["halfLengthYMinus"] = min( lengths["trailLength"], lengths["capotLength"], YRadius )
        lengths["halfLengthYPlus"] = min( lengths["trailLength"], lengths["capotLength"], YRadius )
        lengths["halfLengthZMinus"] = min( lengths["trailLength"], negZRadius )
        lengths["halfLengthZPlus"] = min( lengths["trailLength"], posZRadius )
        # Round to element size
        for key in ["trailLength", "capotLength"]:
            lengths[key] = np.ceil( lengths[key] / elementSizes[0] ) * elementSizes[0]
        for key in ["halfLengthYMinus", "halfLengthYPlus"]:
            lengths[key] = np.ceil( lengths[key] / elementSizes[1] ) * elementSizes[1]
        for key in ["halfLengthZMinus", "halfLengthZPlus"]:
            lengths[key] = np.ceil( lengths[key] / elementSizes[2] ) * elementSizes[2]
        # Define box
        bounds = np.array( [[meshCenter[0] - lengths["trailLength"], meshCenter[0] + lengths["capotLength"]],
               [meshCenter[1] - lengths["halfLengthYMinus"], meshCenter[1] + lengths["halfLengthYPlus"]],
               [meshCenter[2] - lengths["halfLengthZMinus"], meshCenter[2] + lengths["halfLengthZPlus"]]] )
        # Mesh
        tolSearch = 1e-10
        if "toleranceSearches" in self.pFixed.input:
            tolSearch = self.pFixed.input["toleranceSearches"]
        if self.pFixed.domain.dim()==1:
            return meshLine(bounds[:1,:], elementSizes, tolSearch = tolSearch)
        elif self.pFixed.domain.dim()==2:
            return meshRectangle(bounds[:2,:], elementSizes, tolSearch = tolSearch)
        elif self.pFixed.domain.dim()==3:
            return meshBox(bounds, elementSizes, tolSearch = tolSearch)

class LpbfAdaptiveStepper(AdaptiveStepper):
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

    def iterate( self ):
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.setActivation(self.physicalDomain)

        self.setDt()
        #self.setCoupling()

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
            ls.setSolver( self.idxSolverCoupledIter )
            ls.setInitialGuess( self.pMoving, self.pFixed )
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
            self.pFixed.idxSolver = self.idxSolverUncoupledIter
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



def meshLine(box, elSize=[0.25]*1, tolSearch = 1e-10, popup=False):

    gmsh.initialize()
    box = box.reshape(2)
    xMin, xMax = box
    xLen  = xMax - xMin
    nelsX = np.round( xLen / elSize[0] ).astype(int)

    # negativeXFace 
    gmsh.model.geo.addPoint( xMin, 0.0, 0.0, tag = 1 )
    gmsh.model.geo.addPoint( xMax, 0.0, 0.0, tag = 2 )

    line = gmsh.model.geo.addLine(  1, 2, tag=1 )

    gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsX+1)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(1)

    gmsh.model.addPhysicalGroup(1, [line], tag=1, name="Domain")

    if popup:
        gmsh.fltk.run()

    mesh = mhs.gmshModelToMesh( gmsh.model, tolSearch=tolSearch )
    gmsh.finalize()
    return mesh

def meshRectangle(box, elSize=[0.25]*2, recombine = True, tolSearch = 1e-10, popup=False):

    gmsh.initialize()
    box = box.reshape(4)
    xMin, xMax, yMin, yMax = box
    xLen  = xMax - xMin
    yLen  = yMax - yMin
    nelsX = np.round( xLen / elSize[0] ).astype(int)
    nelsY = np.round( yLen / elSize[1] ).astype(int)

    # negativeXFace 
    gmsh.model.geo.addPoint( xMin, yMin, 0.0, tag = 1 )
    gmsh.model.geo.addPoint( xMin, yMax, 0.0, tag = 2 )

    line = gmsh.model.geo.addLine(  1, 2, tag=1 )

    gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsY+1)

    # Extrusion
    _, surface, _, _ = gmsh.model.geo.extrude([(1,line)], xLen, 0.0, 0.0,
                           numElements=[nelsX], recombine=recombine)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    gmsh.model.addPhysicalGroup(2, [surface[-1]], tag=1, name="Domain")

    if popup:
        gmsh.fltk.run()

    mesh = mhs.gmshModelToMesh( gmsh.model, tolSearch=tolSearch )
    gmsh.finalize()
    return mesh

def meshBox(box, elSize=[0.25]*3, recombine = True, tolSearch = 1e-10, popup=False):

    gmsh.initialize()
    box = box.reshape(6)
    xMin, xMax, yMin, yMax, zMin, zMax = box
    xLen  = xMax - xMin
    yLen  = yMax - yMin
    zLen  = zMax - zMin
    nelsX = np.round( xLen / elSize[0] ).astype(int)
    nelsY = np.round( yLen / elSize[1] ).astype(int)
    nelsZ = np.round( zLen / elSize[2] ).astype(int)

    # negativeXFace 
    gmsh.model.geo.addPoint( xMin, yMin, zMin, tag = 1 )
    gmsh.model.geo.addPoint( xMin, yMax, zMin, tag = 2 )
    gmsh.model.geo.addPoint( xMin, yMax, zMax, tag = 3 )
    gmsh.model.geo.addPoint( xMin, yMin, zMax, tag = 4 )

    lines = []
    lines.append( gmsh.model.geo.addLine(  1, 2, tag=1 ) )
    lines.append( gmsh.model.geo.addLine(  2, 3, tag=2 ) )
    lines.append( gmsh.model.geo.addLine(  3, 4, tag=3 ) )
    lines.append( gmsh.model.geo.addLine(  4, 1, tag=4 ) )

    for line in [1, 3]:
        gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsY+1)
    for line in [2, 4]:
        gmsh.model.geo.mesh.setTransfiniteCurve(line, nelsZ+1)

    negativeXContour = gmsh.model.geo.addCurveLoop( lines, 1 )
    negativeXFace = gmsh.model.geo.addPlaneSurface([negativeXContour], 1)

    gmsh.model.geo.mesh.setTransfiniteSurface(negativeXFace)
    gmsh.model.geo.mesh.setRecombine(2,negativeXFace) 

    # Extrusion
    gmsh.model.geo.extrude([(2,negativeXFace)], xLen, 0.0, 0.0,
                           numElements=[nelsX], recombine=recombine)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.model.addPhysicalGroup(3, [1], tag=1, name="Domain")

    if popup:
        gmsh.fltk.run()

    mesh = mhs.gmshModelToMesh( gmsh.model, tolSearch=tolSearch )
    gmsh.finalize()
    return mesh

def deactivateBelowSurface(p, heightOSurface = 0):
    idxHeight = p.domain.dim() - 1
    nels = p.domain.mesh.nels
    activeEls = []
    powderEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElementGeometry( ielem )
        if (e.getCentroid()[idxHeight] < heightOSurface):
            activeEls.append( ielem )
        else:
            powderEls.append( ielem )
    # DEACTIVE
    substrateEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( substrateEls )
    # SET 2 POWDER
    powderEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, powderEls )
    p.domain.setMaterialSets( powderEls )

def activatePowderLayer( p, printer, powderSet = 1 ):
    '''
    If new layer of problem, activate powder layer
    '''
    idxHeight = p.domain.dim() - 1
    currentHeight = p.mhs.position[idxHeight]
    thresholdHeight = currentHeight + p.input["printer"]["height"]

    nels = p.domain.mesh.nels
    newPowderEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElementGeometry( ielem )
        if (e.getCentroid()[idxHeight] < thresholdHeight) and \
        not(p.domain.activeElements[ielem]):
            newPowderEls.append( ielem )
    p.domain.activeElements.setIndices( newPowderEls, 1 )
    p.domain.setActivation( p.domain.activeElements )
    printer.setDepositionTemperature()
    p.domain.materialTag.setIndices( newPowderEls, powderSet )
