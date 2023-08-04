import numpy as np
import pdb
import MovingHeatSource as mhs

#TODO: Store hatch and not adim + nonAdim

class AdaptiveStepper:

    tol = 1e-7

    def __init__(self, pFixed, pMoving, factor=2, adimMaxSubdomainSize=10, threshold= 0.01, isCoupled=True ):
        self.isCoupled = isCoupled
        self.pFixed = pFixed
        self.pMoving = pMoving
        # TODO: better initialization
        self.adimFineDt = 0.25
        self.adimFineSubdomainSize = 0.75 
        self.adimMaxSubdomainSize = adimMaxSubdomainSize
        self.adimMaxDt = adimMaxSubdomainSize / factor
        self.threshold = threshold
        self.factor = factor
        self.adimMinRadius = 1.5

        self.dt = pFixed.dt
        tscale = self.pFixed.mhs.radius / self.pFixed.mhs.path.currentTrack.speed
        self.adimDt = pFixed.dt / tscale
        self.adimSubdomainSize = self.adimFineSubdomainSize
        self.update()

    def update(self):
        '''
        Called at end of iteration
        Computes size of subdomain and dt
        '''
        # TODO: fix this for when track changes from t to tnp1
        
        time = self.pFixed.time

        # If path is over
        if (self.pFixed.mhs.path.isOver(time)):
            return
        tUnit = self.pFixed.mhs.radius / self.pFixed.mhs.path.currentTrack.speed

        dt2trackEnd = self.pFixed.mhs.path.currentTrack.endTime - time
        adimMaxDt = self.adimMaxDt
        if ( dt2trackEnd < 1e-5) or ( time == 0.0 ):#end of track
            self.adimDt = self.adimFineDt
            self.adimSubdomainSize = self.adimFineSubdomainSize
        else:
            adimMaxDt = min( dt2trackEnd/tUnit, self.adimMaxDt )
            delta = self.pMoving.unknown - self.pMoving.previousValues[0]
            metric = delta.getL2Norm() / self.pMoving.unknown.getL2Norm()
            print( " || fn+1 - fn || / || fn+1 || = {} ".format( metric ) )
            if (metric < self.threshold):
                self.adimDt = self.factor * self.adimDt
                self.adimSubdomainSize = self.factor * self.adimSubdomainSize

        # Cap dt
        self.adimDt = min( adimMaxDt, self.adimDt )
        self.dt = tUnit * self.adimDt

    def shapeSubdomain( self ):
        # OBB
        radius = self.pFixed.mhs.radius

        # compute front and sides
        sideRadius = self.adimMinRadius * radius
        adimBackRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize )
        backRadius = max( adimBackRadius, self.adimMinRadius ) * radius
        zRadius    = 1
        xAxis      = self.pFixed.mhs.path.currentTrack.getSpeed() / self.pFixed.mhs.path.currentTrack.speed

        backRadiusObb = max(backRadius - radius, 0.0)
        p0 = self.pMoving.mhs.currentPosition - backRadiusObb*xAxis
        obb = mhs.myOBB( p0, self.pMoving.mhs.currentPosition, 2*sideRadius, 2*zRadius )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( obb )
        collidingElsBackSphere = self.pMoving.domain.mesh.findCollidingElements( p0, self.adimMinRadius*radius )
        collidingElsFrontSphere = self.pMoving.domain.mesh.findCollidingElements( self.pMoving.mhs.currentPosition, self.adimMinRadius*radius )
        subdomainEls += collidingElsBackSphere
        subdomainEls += collidingElsFrontSphere
        subdomain = mhs.MeshTag( self.pMoving.domain.mesh, self.pMoving.domain.mesh.dim, subdomainEls )
        self.pMoving.domain.intersect( subdomain )


    def setDt( self ):
        print(" dt = {}R, domainSize = {}R".format( self.adimDt, self.adimSubdomainSize ) )
        self.pMoving.setDt( self.dt )
        self.pFixed.setDt( self.dt )
        '''
        if (self.adimDt <= 0.25+1e-7):
            self.isCoupled = False
        else:
            self.isCoupled = True
        '''

    def iterate( self ):
        self.setDt()

        # MY SCHEME ITERATE
        self.pFixed.setAssembling2External( True )
        self.pMoving.setAssembling2External( True )
        # PRE-ITERATE AND DOMAIN OPERATIONS
        self.pMoving.domain.resetActivation()
        self.pFixed.domain.resetActivation()

        self.pMoving.intersectExternal( self.pFixed, False )
        self.pFixed.preiterate(False)

        # Update speeds moving domain
        self.pMoving.setAdvectionSpeed( -self.pFixed.mhs.speed )
        self.pMoving.domain.mesh.setSpeedFRF( self.pFixed.mhs.speed )
        self.pMoving.preiterate(False)

        self.pMoving.intersectExternal( self.pFixed, False )
        self.shapeSubdomain()

        if self.isCoupled:
            self.pFixed.substractExternal( self.pMoving, False )
        self.pFixed.updateInterface( self.pMoving )
        self.pMoving.updateInterface( self.pFixed )
        #Dirichet gamma
        self.pFixed.setGamma2Dirichlet()
        # Pre-assembly, updating free dofs
        self.pMoving.preAssemble(True)
        self.pFixed.preAssemble(True)
        ls = mhs.LinearSystem( self.pMoving, self.pFixed )
        ls.cleanup()
        # Assembly
        self.pMoving.assemble()
        self.pFixed.assemble()
        # Assembly Gamma
        self.pFixed.assembleDirichletGamma( self.pMoving )
        self.pMoving.assembleNeumannGamma( self.pFixed )
        # Build ls
        ls.assemble()
        # Solve ls
        ls.solve()
        # Recover solution
        self.pFixed.gather()
        self.pMoving.gather()
        self.pFixed.unknown.interpolateInactive( self.pMoving.unknown, False )
        if self.isCoupled:
            self.pMoving.unknown.interpolateInactive( self.pFixed.unknown, True )
        else:
            self.pMoving.unknown.interpolate( self.pFixed.unknown )

        # Post iteration
        self.pFixed.postIterate()
        self.pMoving.postIterate()

        # Compute next dt && domain size
        self.update()

        self.pFixed.writepos(
            nodeMeshTags={ "gammaNodes":self.pFixed.gammaNodes, },)
        self.pMoving.writepos(
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )
