import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper

class CustomStepper(AdaptiveStepper):
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
            self.computeSteadinessMetric()
            if (self.checkSteadinessCriterion()):
                self.adimDt = min( self.factor * self.adimDt, self.adimDt + 1 )

        # Cap dt
        self.adimDt = min( adimMaxDt, self.adimDt )
        self.adimSubdomainSize = self.adimDt + 1.0
        self.dt = tUnit * self.adimDt

        if self.onNewTrack:
            self.nextTrack = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time + self.dt )


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
        zRadius    = 1
        xAxis      = self.nextTrack.getSpeed() / self.nextTrack.speed

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

