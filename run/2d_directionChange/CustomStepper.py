import numpy as np
import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper

#TODO: Store hatch and not adim + nonAdim

class CustomStepper( AdaptiveStepper ):

    def setBCs( self ):
        self.pFixed.setConvection( resetBcs = False )
        if self.isCoupled:
            self.pMoving.setConvection( resetBcs = False )

    def computeSizeSubdomain( self, adimDt = None ):
        if adimDt is None:
            adimDt = self.adimDt
        adimSize = np.round(4*(1 + adimDt + (adimDt**2.1)*self.adimFineDt)) / 4
        adimSize = min( adimSize, self.pFixed.input["maxAdimSize"] )
        return adimSize

    def increaseDt( self ):
        self.adimDt += self.adimFineDt

    '''
    def onNewTrackOperations(self):
        self.rotateSubdomain()
        self.isCoupled = False
        speed = self.nextTrack.getSpeed()
        self.pMoving.setAdvectionSpeed( -speed )
        self.pMoving.domain.setSpeed( speed )
        self.pMoving.mhs.setPower( self.nextTrack.power )
    '''

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
        obb = mhs.MyOBB( p0, self.pMoving.mhs.position, 2*sideRadius, 2*zRadius )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( obb )
        collidingElsBackSphere = self.pMoving.domain.mesh.findCollidingElements( p0, self.adimMinRadius*radius )
        collidingElsFrontSphere = self.pMoving.domain.mesh.findCollidingElements( self.pMoving.mhs.position, self.adimMinRadius*radius )
        subdomainEls += collidingElsBackSphere
        subdomainEls += collidingElsFrontSphere
        subdomain = mhs.MeshTag( self.pMoving.domain.mesh, self.pMoving.domain.mesh.dim, subdomainEls )
        self.pMoving.domain.intersect( subdomain )
