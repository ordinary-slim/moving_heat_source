import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper
import numpy as np
import meshzoo

#TODO: Store hatch and not adim + nonAdim
#TODO: Build pMoving here in order to reduce risk of mistake

class CustomStepper(AdaptiveStepper):

    tol = 1e-7

    def __init__(self,
                 pFixed,
                 factor=2,
                 maxAdimtDt=8,
                 threshold= 0.01,
                 elementSize=0.25,
                 adimMinRadius=2,
                 adimZRadius=None,
                 isCoupled=True, ):

        self.pFixed = pFixed
        # Set up laser path
        if "path" in self.pFixed.input:
            self.pFixed.setPath( self.pFixed.input["path"] )
        self.isCoupled = isCoupled
        self.adimMaxDt = maxAdimtDt
        self.adimMaxSubdomainSize = self.adimMaxDt + 1 # ad-hoc
        self.buildMovingProblem(elementSize=elementSize)
        # TODO: better initialization
        self.adimFineDt = 0.5
        self.adimFineSubdomainSize = 1
        self.threshold = threshold
        self.factor = factor
        self.adimMinRadius = adimMinRadius
        if not(adimZRadius is None):
            self.adimZRadius = adimZRadius
        else:
            self.adimZRadius = adimMinRadius

        self.dt = pFixed.dt
        tscale = self.pFixed.mhs.radius / self.pFixed.mhs.currentTrack.speed
        self.adimDt = pFixed.dt / tscale
        self.adimSubdomainSize = self.adimFineSubdomainSize
        self.update()
        self.onNewTrack = True
        self.hasPrinter = False
        self.physicalDomain = mhs.MeshTag( self.pFixed.domain.activeElements )

        # Match heat source positions at t = 0.0
        self.pMoving.mhs.setPosition( self.pFixed.mhs.position )

        self.nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time + self.dt )
        # Set up printer
        if "printer" in self.pFixed.input:
            self.setupPrinter( pFixed.input["printer"]["width"], pFixed.input["printer"]["height"], pFixed.input["printer"]["depth"] )

        # Initialize steadiness criterion
        self.metric = 1e99


    def buildMovingProblem(self, elementSize=0.25):
        meshInputMoving = {}
        meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(self.adimMaxSubdomainSize, self.pFixed, elementSize=elementSize)
        self.meshMoving = mhs.Mesh(meshInputMoving)
        movingProblemInput = dict( self.pFixed.input )
        movingProblemInput["isStabilized"] = 1
        movingProblemInput["isAdvection"] = 1
        movingProblemInput["advectionSpeed"] = -self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["speedDomain"]      = self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["HeatSourceSpeed"] = np.zeros(3)
        self.pMoving  = mhs.Problem(self.meshMoving, movingProblemInput, caseName="moving")
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

    def setSizeSubdomain( self ):
        self.adimSubdomainSize = min(self.adimDt + 2.0, self.adimDt * 2 )
        return self.adimSubdomainSize

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
        if ( dt2trackEnd < 1e-7) or ( time == 0.0 ):#new track
            self.onNewTrack = True
            self.adimDt = self.adimFineDt
            self.adimSubdomainSize = self.adimFineSubdomainSize
        else:
            adimMaxDt = min( dt2trackEnd/tUnit, self.adimMaxDt )
            self.computeSteadinessMetric(verbose=True)
            if (self.checkSteadinessCriterion()):
                self.adimDt = min( self.factor * self.adimDt, self.adimDt + 1 )

        # Cap dt
        self.adimDt = min( adimMaxDt, self.adimDt )
        self.setSizeSubdomain()
        self.dt = tUnit * self.adimDt

        if self.onNewTrack:
            self.nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time + self.dt )

    def rotateSubdomain( self, currentOrientation=None, nextOrientation=None ):
        '''
        Align mesh of moving subproblem with track
        Called at tn
        '''
        if currentOrientation is None:
            currentOrientation = self.pFixed.mhs.currentTrack.getSpeed()
        if nextOrientation is None:
            nextOrientation    = self.nextTrack.getSpeed()
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
        backRadius = max( adimBackRadius, self.adimMinRadius ) * radius
        zRadius    = self.adimZRadius * radius
        xAxis      = self.nextTrack.getSpeed() / self.nextTrack.speed

        #backRadiusObb = max(backRadius - radius, 0.0)
        p0 = self.pMoving.mhs.position - backRadius*xAxis
        p1 = self.pMoving.mhs.position + self.adimMinRadius*xAxis
        obb = mhs.myOBB( p0, p1, 2*sideRadius, 2*zRadius )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( obb )
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
        # Set coupling
        if ( ((self.adimDt <= 0.5+1e-7) and not(self.isCoupled)) \
            or not(self.nextTrack.hasDeposition) ):
            self.isCoupled = False
        else:
            self.isCoupled = True

    def writepos( self ):
        activeInExternal = self.pFixed.getActiveInExternal( self.pMoving, 1e-7 )
        self.pFixed.writepos(
            nodeMeshTags={
                "gammaNodes":self.pFixed.gammaNodes,
                "forcedDofs":self.pFixed.forcedDofs,
                "activeInExternal":activeInExternal,
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
        if self.onNewTrack:
            self.onNewTrackOperations()

        self.shapeSubdomain()

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#tn intersect

        # Motion, other operations
        self.pMoving.preiterate(canPreassemble=False)
        self.pFixed.preiterate(canPreassemble=False)

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
        self.update()
        # Write vtk files
        self.writepos()

def meshRectangle(box, elementSize=0.25):
    cell_type="quad4"
    nelsX = round( (box[0][1]-box[0][0]) / elementSize )
    nelsY = round( (box[1][1]-box[1][0]) / elementSize )
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0][0], box[0][1], nelsX+1),
        np.linspace(box[1][0], box[1][1], nelsY+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def meshBox(box, elementSize=0.25):
    cell_type="hexa8"
    nelsX = round( (box[0][1]-box[0][0]) / elementSize )
    nelsY = round( (box[1][1]-box[1][0]) / elementSize )
    nelsZ = round( (box[2][1]-box[2][0]) / elementSize )

    points, cells = meshzoo.cube_hexa(
        np.linspace( box[0][0], box[0][1], nelsX+1),
        np.linspace( box[1][1], box[1][0], nelsY+1),
        np.linspace( box[2][1], box[2][0], nelsZ+1),
    )
    cells = cells.astype( np.uint32 )
    return points, cells, cell_type

def meshAroundHS( adimR, pFixed, elementSize=0.25 ):
    radius = pFixed.mhs.radius
    initialPosition = pFixed.mhs.position
    lengths = {}
    lengths["trailLength"] = adimR * radius
    lengths["capotLength"] = min( lengths["trailLength"], 2*radius )
    lengths["halfLengthY"] = min( lengths["trailLength"], lengths["capotLength"] )
    lengths["halfLengthZ"] = min( lengths["trailLength"], lengths["capotLength"] )
    # Round to element size
    for key in lengths.keys():
        lengths[key] = np.ceil( lengths[key] / elementSize ) * elementSize
    bounds = np.array( [[initialPosition[0] - lengths["trailLength"], initialPosition[0] + lengths["capotLength"]],
           [initialPosition[1] - lengths["halfLengthY"], initialPosition[1] + lengths["halfLengthY"]],
           [initialPosition[2] - lengths["halfLengthZ"], initialPosition[2] + lengths["halfLengthZ"]]] )
    if pFixed.domain.dim()==2:
        return meshRectangle(bounds, elementSize)
    elif pFixed.domain.dim()==3:
        return meshBox(bounds, elementSize)

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

    def setDtFromAdimR( self, adimR, maxDt=None):
        r =self.problem.input["radius"]
        speed = max( np.linalg.norm(self.problem.input["HeatSourceSpeed"] ), np.linalg.norm(self.problem.input["advectionSpeed"] ) )
        dt =  adimR * r / speed 
        if maxDt:
            dt = min(dt, maxDt)
        self.problem.setDt( dt )

    def iterate( self ):
        self.setDtFromAdimR( 0.5, self.dt2trackEnd )
        tnp1 = self.problem.time + self.problem.dt
        track = self.problem.mhs.path.interpolateTrack( tnp1 )
        if (track.hasDeposition) and (self.printer is not None):
            self.printer.deposit( self.problem.mhs.position,
                                self.problem.mhs.path.interpolatePosition(tnp1)
                               )
        self.problem.iterate()
        self.dt2trackEnd = self.problem.mhs.currentTrack.endTime - self.problem.time
        self.problem.writepos()
