import numpy as np
import pdb
import MovingHeatSource as mhs

#TODO: Store hatch and not adim + nonAdim

class AdaptiveStepper:
    tol = 1e-7
    def __init__(self, pFixed, pMoving, factor=2, adimMaxSubdomainSize=10, threshold= 0.01 ):
        self.pFixed = pFixed
        self.pMoving = pMoving
        # TODO: better initialization
        self.adimFineDt = 0.5
        self.adimFineSubdomainSize = 1
        self.adimMaxSubdomainSize = adimMaxSubdomainSize
        self.adimMaxDt = adimMaxSubdomainSize / factor
        self.threshold = threshold
        self.factor = factor

        self.dt = pFixed.dt
        self.currentTrack = pFixed.mhs.path.currentTrack
        tscale = self.pFixed.mhs.radius / self.currentTrack.speed
        self.adimDt = pFixed.dt / tscale
        self.adimSubdomainSize = self.adimFineSubdomainSize
        self.update()

    def update(self):
        # TODO: fix this for when track changes from t to tnp1
        
        self.currentTrack = self.pFixed.mhs.path.currentTrack
        time = self.pFixed.time

        if (self.pFixed.mhs.path.isOver(time)):
            return
        tUnit = self.pFixed.mhs.radius / self.currentTrack.speed

        dt2trackEnd = self.currentTrack.endTime - time
        print( "dt2trackEnd = ", dt2trackEnd )
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

        if not(self.pFixed.mhs.path.isOver(time)):
            self.getObbSubdomain()

    def getObbSubdomain( self ):
        # Given currentSubdomainSize and currentDt
        # compute OBB of subdomain
        track = self.pFixed.mhs.path.interpolateTrack( self.pFixed.time + self.dt )
        sUnit = self.pFixed.mhs.radius
        tnp1  = self.pFixed.time + self.dt
        heatSourcePosition = self.pFixed.mhs.path.interpolatePosition( tnp1 ) \
                + self.pFixed.domain.mesh.shiftFRF - self.pMoving.domain.mesh.shiftFRF

        # compute front and sides
        frontRadius = min( 1, self.adimSubdomainSize ) * sUnit
        sideRadius = min( 1, self.adimSubdomainSize ) * sUnit
        backRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize ) * sUnit
        zRadius    = 1
        xAxis      = track.getSpeed() / track.speed
        #TODO: Fix change ref frame here
        self.OBB   = mhs.myOBB( heatSourcePosition - backRadius*sUnit*xAxis,
                                heatSourcePosition + frontRadius*sUnit*xAxis,
                                2*sideRadius,
                                2*zRadius )

    def setDt( self ):
        self.pMoving.setDt( self.dt )
        self.pFixed.setDt( self.dt )
