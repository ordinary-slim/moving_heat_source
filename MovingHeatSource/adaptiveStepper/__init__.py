import numpy as np
import meshzoo
import MovingHeatSource as mhs

#TODO: Store hatch and not adim + nonAdim
#TODO: Build pMoving here in order to reduce risk of mistake

class AdaptiveStepper:

    tol = 1e-7

    def __init__(self,
                 pFixed,
                 factor=2,
                 adimMaxSubdomainSize=10,
                 threshold= 0.01,
                 isCoupled=True, ):

        self.isCoupled = isCoupled
        self.pFixed = pFixed
        self.adimMaxSubdomainSize = adimMaxSubdomainSize
        self.buildMovingProblem()
        # TODO: better initialization
        self.adimFineDt = 0.5
        self.adimFineSubdomainSize = 1
        self.adimMaxDt = adimMaxSubdomainSize / factor
        self.threshold = threshold
        self.factor = factor
        self.adimMinRadius = 2

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

    def buildMovingProblem(self):
        adimR = self.pFixed.input["radius"] / np.linalg.norm(self.pFixed.input["HeatSourceSpeed"])
        meshInputMoving = {}
        meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(self.adimMaxSubdomainSize, self.pFixed)
        self.meshMoving = mhs.Mesh(meshInputMoving)
        movingProblemInput = dict( self.pFixed.input )
        movingProblemInput["isAdvection"] = 1
        movingProblemInput["advectionSpeed"] = -self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["speedFRF"]      = self.pFixed.input["HeatSourceSpeed"]
        movingProblemInput["HeatSourceSpeed"] = np.zeros(3)
        self.pMoving  = mhs.Problem(self.meshMoving, movingProblemInput, caseName="moving")

    def getTime(self):
        return self.pFixed.time

    def setupPrinter( self, mdwidth, mdheight ):
        self.printerFixed = mhs.Printer( self.pFixed, mdwidth, mdheight )
        self.printerMoving = mhs.Printer( self.pMoving, mdwidth, mdheight )
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
            delta = self.pMoving.unknown - self.pMoving.previousValues[0]
            metric = delta.getL2Norm() / self.pMoving.unknown.getL2Norm()
            print( " || fn+1 - fn || / || fn+1 || = {} ".format( metric ) )
            if (metric < self.threshold):
                self.adimDt = min( self.factor * self.adimDt, self.adimDt + 1 )

        # Cap dt
        self.adimDt = min( adimMaxDt, self.adimDt )
        self.adimSubdomainSize = self.adimDt + 1.0
        self.dt = tUnit * self.adimDt

        if self.onNewTrack:
            self.nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time + self.dt )

    def rotateSubdomain( self ):
        '''
        Align mesh of moving subproblem with track
        Called at tn
        '''
        currentTrack = self.pFixed.mhs.currentTrack
        center = self.pMoving.mhs.position
        angle = np.arccos( np.dot( self.nextTrack.getSpeed(), currentTrack.getSpeed() ) / self.nextTrack.speed / currentTrack.speed )
        if (angle > 1e-5):
            self.pMoving.domain.inPlaneRotate( center, angle )
            self.pMoving.unknown.interpolate( self.pFixed.unknown, ignoreOutside = True )

    def shapeSubdomain( self ):
        '''
        At t^n, do things
        '''
        # OBB
        radius = self.pFixed.mhs.radius

        # compute front and sides
        sideRadius = self.adimMinRadius * radius
        adimBackRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize )
        backRadius = max( adimBackRadius, self.adimMinRadius ) * radius
        zRadius    = 1
        xAxis      = self.nextTrack.getSpeed() / self.nextTrack.speed

        backRadiusObb = max(backRadius - radius, 0.0)
        p0 = self.pMoving.mhs.position - backRadiusObb*xAxis
        obb = mhs.myOBB( p0, self.pMoving.mhs.position, 2*sideRadius, 2*zRadius )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( obb )
        collidingElsBackSphere = self.pMoving.domain.mesh.findCollidingElements( p0, self.adimMinRadius*radius )
        collidingElsFrontSphere = self.pMoving.domain.mesh.findCollidingElements( self.pMoving.mhs.position, self.adimMinRadius*radius )
        subdomainEls += collidingElsBackSphere
        subdomainEls += collidingElsFrontSphere
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
        # New track operations
        if (self.onNewTrack):
            #self.rotateSubdomain()
            self.isCoupled = False
            speed = self.nextTrack.getSpeed()
            print("speed = {}".format( speed ) )
            self.pMoving.setAdvectionSpeed( -speed )
            self.pMoving.domain.setSpeed( speed )
            self.pMoving.mhs.setPower( self.nextTrack.power )

    def iterate( self ):
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.setActivation(self.physicalDomain)

        self.setDt()
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
            #Dirichet gamma
            self.pFixed.setGamma2Dirichlet()
            self.pMoving.setGamma2Neumann()

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
        else:
            self.pFixed.clearGamma()
            self.pFixed.preAssemble(allocateLs=True)
            self.pFixed.iterate()
            self.pMoving.unknown.interpolate( self.pFixed.unknown, ignoreOutside = False )

        # Post iteration
        self.pFixed.postIterate()
        self.pMoving.postIterate()

        activeInExternal = self.pFixed.getActiveInExternal( self.pMoving, 1e-7 )

        # Compute next dt && domain size
        self.update()

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

def meshBox(box, meshDen=4):
    cell_type="quad4"
    nelsX = int(meshDen*(box[1]-box[0])) +1
    nelsY = int(meshDen*(box[3]-box[2])) +1
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], nelsX),
        np.linspace(box[2], box[3], nelsY),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def meshAroundHS( adimR, pFixed, meshDen=4 ):
    radius = pFixed.mhs.radius
    initialPosition = pFixed.mhs.position
    trailLength = adimR * radius
    capotLength = min( trailLength, 2*radius )
    halfLengthY = min( trailLength, capotLength )
    box = [initialPosition[0] - trailLength, initialPosition[0] + capotLength,
           initialPosition[1] - halfLengthY, initialPosition[1] + halfLengthY]
    return meshBox(box, meshDen)

