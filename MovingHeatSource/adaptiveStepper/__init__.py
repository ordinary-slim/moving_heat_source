import numpy as np
import gmsh
import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
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
                 adimZRadius=None,
                 isCoupled=True,
                 ):

        self.pFixed = pFixed
        self.isAdaptive = isAdaptive
        self.adimFineDt = adimFineDt
        # Set up laser path
        if "path" in self.pFixed.input:
            self.pFixed.setPath( self.pFixed.input["path"] )
        self.isCoupled = isCoupled
        self.adimMaxDt = maxAdimtDt
        self.adimMaxSubdomainSize = self.computeSizeSubdomain(self.adimMaxDt)
        self.adimMinRadius = adimMinRadius
        self.adimZRadius = adimMinRadius
        if not(adimZRadius is None):
            self.adimZRadius = adimZRadius
        self.buildMovingProblem(elementSize=elementSize, shift=shift)
        # TODO: better initialization
        self.threshold = threshold
        self.factor = factor

        self.dt = pFixed.dt
        tscale = self.pFixed.mhs.radius / self.pFixed.mhs.currentTrack.speed
        self.adimDt = pFixed.dt / tscale
        self.nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time )
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


    def buildMovingProblem(self, elementSize, shift):
        meshInputMoving = {}
        try:
            len(elementSize)
        except TypeError:
            elementSize = [elementSize]*3
        self.meshMoving = self.meshAroundHS(elementSizes=elementSize, shift=shift)
        movingProblemInput = dict( self.pFixed.input )
        movingProblemInput["isStabilized"] = 1
        movingProblemInput["isAdvection"] = 1
        movingProblemInput["advectionSpeed"] = -self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["speedDomain"]      = self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["HeatSourceSpeed"] = np.zeros(3)
        movingProblemInput["initialPosition"] = self.pFixed.mhs.position
        self.pMoving  = mhs.Problem(self.meshMoving, movingProblemInput, caseName=self.pFixed.caseName + "_moving")
        # Quick-fix: Mesh is built oriented around X
        self.rotateSubdomain( currentOrientation=np.array([1.0, 0.0, 0.0]),
                               nextOrientation=self.pFixed.mhs.currentTrack.getSpeed() )

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
        if (self.pFixed.mhs.currentTrack.hasDeposition):
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

    def setSizeSubdomain( self ):
        self.adimSubdomainSize = self.computeSizeSubdomain()

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
            if (maxDt2TrackEnd > self.adimMinRadius + 1e-7):
                maxDt2TrackEnd -= self.adimMinRadius
            else:
                maxDt2TrackEnd = min( self.adimFineDt, adimDt2TrackEnd )

            adimMaxDt = min( maxDt2TrackEnd, self.adimMaxDt )
            self.computeSteadinessMetric(verbose=True)
            if (self.checkSteadinessCriterion()):
                self.increaseDt()

        # Cap dt
        self.adimDt = min( adimMaxDt, self.adimDt )
        self.setSizeSubdomain()
        self.dt = tUnit * self.adimDt

    def rotateSubdomain( self, currentOrientation=None, nextOrientation=None ):
        '''
        Align mesh of moving subproblem with track
        Called at tn
        '''
        if currentOrientation is None:
            currentOrientation = self.pFixed.mhs.currentTrack.getSpeed()
        if nextOrientation is None:
            nextOrientation    = self.nextTrack.getSpeed()
        # If z-motion, do nothing
        if (np.linalg.norm( np.cross(nextOrientation, np.array([0.0, 0.0, 1.0])) ) < 1e-7):
            return
        center = self.pMoving.mhs.position
        angle = np.arccos( np.dot( nextOrientation, currentOrientation ) / np.linalg.norm( nextOrientation ) / np.linalg.norm( currentOrientation ) )
        if (angle > 1e-5):
            self.pMoving.domain.inPlaneRotate( center, angle )
            self.pMoving.unknown.interpolate( self.pFixed.unknown, ignoreOutside = True )

    def shapeSubdomain( self ):
        '''
        At t^n, do things
        '''
        if not(self.nextTrack.hasDeposition):
            return
        # OBB
        radius = self.pFixed.mhs.radius

        # compute front and sides
        sideRadius = self.adimMinRadius * radius
        adimBackRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize )
        backRadiusAabb = max( adimBackRadius, self.adimMinRadius ) * radius
        frontRadiusAabb = self.adimMinRadius*radius
        zRadius    = self.adimZRadius * radius
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
        # Called at tn
        print(" dt = {}R, domainSize = {}R".format( self.adimDt, self.adimSubdomainSize ) )
        self.pMoving.setDt( self.dt )
        self.pFixed.setDt( self.dt )

    def setCoupling( self ):
        # Set coupling
        if ( ((self.adimDt <= self.adimFineDt+1e-7) and not(self.isCoupled)) \
            or not(self.nextTrack.hasDeposition) ):
            self.isCoupled = False
        else:
            self.isCoupled = True

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
        self.rotateSubdomain()
        self.isCoupled = False
        speed = self.nextTrack.getSpeed()
        self.pMoving.setAdvectionSpeed( -speed )
        self.pMoving.domain.setSpeed( speed )
        self.pMoving.mhs.setPower( self.nextTrack.power )


    def setBCs(self):
        pass

    def iterate( self ):
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.setActivation(self.physicalDomain)

        self.setDt()
        self.setCoupling()
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
        minRadius = self.adimMinRadius * radius
        zRadius   = self.adimZRadius * radius
        meshCenter = np.array(self.pFixed.mhs.position)
        if not(shift is None):
            meshCenter += shift
        lengths = {}
        lengths["trailLength"] = self.adimMaxSubdomainSize * radius
        lengths["capotLength"] = min( lengths["trailLength"], minRadius )
        lengths["halfLengthY"] = min( lengths["trailLength"], lengths["capotLength"] )
        lengths["halfLengthZ"] = min( lengths["trailLength"], zRadius )
        # Round to element size
        for key in ["trailLength", "capotLength"]:
            lengths[key] = np.ceil( lengths[key] / elementSizes[0] ) * elementSizes[0]
        lengths["halfLengthY"] = np.ceil( lengths["halfLengthY"] / elementSizes[1] ) * elementSizes[1]
        lengths["halfLengthZ"] = np.ceil( lengths["halfLengthZ"] / elementSizes[2] ) * elementSizes[2]
        # Define box
        bounds = np.array( [[meshCenter[0] - lengths["trailLength"], meshCenter[0] + lengths["capotLength"]],
               [meshCenter[1] - lengths["halfLengthY"], meshCenter[1] + lengths["halfLengthY"]],
               [meshCenter[2] - lengths["halfLengthZ"], meshCenter[2] + lengths["halfLengthZ"]]] )
        # Mesh
        tolSearch = 1e-10
        if "toleranceSearches" in self.pFixed.input:
            tolSearch = self.pFixed.input["toleranceSearches"]
        if self.pFixed.domain.dim()==2:
            return meshRectangle(bounds[:2,:], elementSizes, tolSearch = tolSearch)
        elif self.pFixed.domain.dim()==3:
            return meshBox(bounds, elementSizes, tolSearch = tolSearch)

def meshRectangle(box, elSize=[0.25]*2, tolSearch = 1e-10, popup=False):

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
                           numElements=[nelsX], recombine=True)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    gmsh.model.addPhysicalGroup(2, [surface[-1]], tag=1, name="Domain")

    if popup:
        gmsh.fltk.run()

    mesh = mhs.gmshModelToMesh( gmsh.model, tolSearch=tolSearch )
    gmsh.finalize()
    return mesh

def meshBox(box, elSize=[0.25]*3, tolSearch = 1e-10, popup=False):

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
                           numElements=[nelsX], recombine=True)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.model.addPhysicalGroup(3, [1], tag=1, name="Domain")

    if popup:
        gmsh.fltk.run()

    mesh = mhs.gmshModelToMesh( gmsh.model, tolSearch=tolSearch )
    gmsh.finalize()
    return mesh

