import numpy as np
import pdb
import MovingHeatSource as mhs

#TODO: Store hatch and not adim + nonAdim

class AdaptiveStepper:

    tol = 1e-7

    def __init__(self, pFixed, pMoving, factor=2, adimMaxSubdomainSize=10, threshold= 0.01, isCoupled=True, rotateSubdomain=False ):
        self.isCoupled = isCoupled
        self.pFixed = pFixed
        self.pMoving = pMoving
        # TODO: better initialization
        self.adimFineDt = 0.5
        self.adimFineSubdomainSize = 1
        self.adimMaxSubdomainSize = adimMaxSubdomainSize
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
        # Initialize physical domain
        self.physicalDomain = mhs.MeshTag( self.pFixed.domain.activeElements )

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
            self.adimDt = self.adimFineDt
            self.adimSubdomainSize = self.adimFineSubdomainSize
            self.onNewTrack = True
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

    def rotateSubdomain( self ):
        '''
        Align mesh of moving subproblem with track
        '''
        nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time + self.dt )
        currentTrack = self.pFixed.mhs.currentTrack
        center = self.pMoving.mhs.position
        angle = np.arccos( np.dot( nextTrack.getSpeed(), currentTrack.getSpeed() ) / nextTrack.speed / currentTrack.speed )
        if (angle > 1e-5):
            self.pMoving.domain.inPlaneRotate( center, angle )
            self.pMoving.unknown.interpolate( self.pFixed.unknown, ignoreOutside = True )

    def shapeSubdomain( self ):
        '''
        At t^n, do things
        '''
        nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time + self.dt )
        # OBB
        radius = self.pFixed.mhs.radius

        # compute front and sides
        sideRadius = self.adimMinRadius * radius
        adimBackRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize )
        backRadius = max( adimBackRadius, self.adimMinRadius ) * radius
        zRadius    = 1
        xAxis      = nextTrack.getSpeed() / nextTrack.speed

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
        print(" dt = {}R, domainSize = {}R".format( self.adimDt, self.adimSubdomainSize ) )
        self.pMoving.setDt( self.dt )
        self.pFixed.setDt( self.dt )
        # New track operations
        if (self.onNewTrack):
            self.rotateSubdomain()
            self.isCoupled = False
            self.pMoving.setAdvectionSpeed( -self.pFixed.mhs.speed )
            self.pMoving.domain.setSpeed( self.pFixed.mhs.speed )

        # Set coupling
        if (self.adimDt <= 0.5+1e-7) and not(self.isCoupled):
            self.isCoupled = False
        else:
            self.isCoupled = True

    def iterate( self ):
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.setActivation(self.physicalDomain)

        self.setDt()
        self.shapeSubdomain()

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#tn intersect

        # Motion, other operations
        self.pMoving.preiterate( canPreassemble=False )
        self.pFixed.preiterate( canPreassemble=False )

        self.pMoving.intersectExternal(self.pFixed, updateGamma=False)#physical domain intersect

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
            self.pMoving.unknown.interpolateInactive( self.pFixed.unknown, ignoreOutside = True )
        else:
            self.pFixed.preAssemble(allocateLs=True)
            self.pFixed.iterate()
            self.pMoving.unknown.interpolate( self.pFixed.unknown, ignoreOutside = True )

        # Post iteration
        self.pFixed.postIterate()
        self.pMoving.postIterate()

        activeInExternal = self.pFixed.getActiveInExternal( self.pMoving, 1e-7 )

        # Compute next dt && domain size
        self.update()

        self.pFixed.writepos(
            nodeMeshTags={
                "gammaNodes":self.pFixed.gammaNodes,
                "activeInExternal":activeInExternal,
                },
            cellMeshTags={
                "physicalDomain":self.physicalDomain,
                },
                          )
        self.pMoving.writepos(
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )
