import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper
from MovingHeatSource.cpp import TrackType
import numpy as np
import meshzoo
import pdb
import line_profiler

#TODO: Store hatch and not adim + nonAdim
#TODO: Build pMoving here in order to reduce risk of mistake

class CustomStepper(AdaptiveStepper):

    valuesOfMetric = [1e9]
    lastVals4mean = 2
    maximumTs = [np.inf]*(lastVals4mean+1)

    def shapeSubdomain( self ):
        '''
        At t^n, do things
        '''
        if not(self.nextTrack.type != TrackType.printing):
            return
        # OBB
        radius = self.pFixed.mhs.radius

        # compute front and sides
        sideRadius = self.adimMinRadius * radius
        adimBackRadius = min( self.adimMaxSubdomainSize, self.adimSubdomainSize )
        backRadius = max( adimBackRadius, self.adimMinRadius ) * radius
        zRadius    = self.adimZRadius * radius
        xAxis      = self.nextTrack.getSpeed() / self.nextTrack.speed

        backRadiusAabb = max(backRadius - radius, 0.0)
        frontRadiusAabb = self.adimMinRadius*radius
        pback = self.pMoving.mhs.position - backRadiusAabb*xAxis
        pfront = self.pMoving.mhs.position + self.adimMinRadius*xAxis
        paabb  = self.pMoving.mhs.position + (frontRadiusAabb - backRadiusAabb)*xAxis/2
        hLens = np.array( [ (frontRadiusAabb + backRadiusAabb)/2.0,
                            sideRadius,
                            zRadius,] )
        aabb = mhs.MyAABB( paabb, hLens, True )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( aabb )
        collidingElsBackSphere = self.pMoving.domain.mesh.findCollidingElements( pback, self.adimMinRadius*radius )
        #collidingElsFrontSphere = self.pMoving.domain.mesh.findCollidingElements( self.pMoving.mhs.position, self.adimMinRadius*radius )
        subdomainEls += collidingElsBackSphere
        #subdomainEls += collidingElsFrontSphere
        subdomain = mhs.MeshTag( self.pMoving.domain.mesh, self.pMoving.domain.mesh.dim, subdomainEls )
        self.pMoving.domain.intersect( subdomain )


    def increaseDt( self ):
        self.adimDt = self.adimDt
        for _ in range( self.lastVals4mean ):
            self.maximumTs.append( np.inf )

    def setBCs( self ):
        pass

    def computeSizeSubdomain( self, adimDt = None ):
        if adimDt is None:
            adimDt = self.adimDt
        return 10

    def getIsPrinting( self ):
        return (self.pFixed.mhs.currentTrack.type == TrackType.printing)

    def getNdofs( self ):
        return self.pFixed.ls.ndofs()

    def computeSteadinessMetric( self, verbose=True ):
        # TODO: Change this
        delta = self.pMoving.unknown - self.pMoving.previousValues[0]
        self.metric = delta.getL2Norm() / self.pMoving.unknown.getL2Norm()
        self.maximumTs.append( max(self.pMoving.unknown.values) )
        movingAv = np.mean( np.array( self.maximumTs[(-1-self.lastVals4mean):-1] ) )
        self.maxTChange = abs(self.maximumTs[-1] - movingAv ) / movingAv
        self.valuesOfMetric.append( self.metric )
        self.relChangeMetric = abs( (self.valuesOfMetric[-1] - self.valuesOfMetric[-2]) / self.valuesOfMetric[-1] )
        if verbose:
            print( " || fn+1 - fn || / || fn+1 || = {} ".format( self.valuesOfMetric[-1] ) )
            print( " threshold = {}".format( self.threshold ) )
            print( "Relative change metric = {}%".format( 100*self.relChangeMetric ) )
            try:
                print( "Previous max T = {}, current max T = {}".format( self.maximumTs[-2], self.maximumTs[-1] ) )
            except IndexError:
                pass
            print( "Relative change max T change over last iters = {}%".format( 100*self.maxTChange ) )

    def checkSteadinessCriterion(self):
        return (((self.valuesOfMetric[-1] < self.threshold) or (self.relChangeMetric < 0.01)) \
                and (self.maxTChange < 0.05))

    def writepos( self ):
        self.pFixed.writepos(
            shift=-self.pFixed.mhs.position,
            nodeMeshTags={
                "gammaNodes":self.pFixed.gammaNodes,
                "forcedDofs":self.pFixed.forcedDofs,
                },
            cellMeshTags={
                "physicalDomain":self.physicalDomain,
                },
                          )
        self.pMoving.writepos(
            shift=-self.pFixed.mhs.position,
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )




class DriverReference:
    def __init__(self, problem, adimDt=0.5, matlab=None):
        self.problem = problem
        if "path" in self.problem.input:
            self.problem.setPath( self.problem.input["path"] )
        self.dt2trackEnd = None
        self.printer = None
        self.adimDt  = adimDt
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
        r =self.problem.input["radius"]
        speed = max( np.linalg.norm(self.problem.input["HeatSourceSpeed"] ), np.linalg.norm(self.problem.input["advectionSpeed"] ) )
        dt =  adimR * r / speed 
        if maxDt:
            dt = min(dt, maxDt)
        self.problem.setDt( dt )

    def solve(self):
        self.problem.iterate()
    def writepos(self):
        self.problem.writepos(
                shift=-self.problem.mhs.position,
                )

    def iterate( self ):
        self.setDtFromAdimR( self.adimDt, self.dt2trackEnd )
        tnp1 = self.problem.time + self.problem.dt
        track = self.problem.mhs.path.interpolateTrack( tnp1 )
        #self.problem.setConvection( resetBcs = True )
        self.solve()
        self.dt2trackEnd = self.problem.mhs.currentTrack.endTime - self.problem.time
        if self.dt2trackEnd < self.tol:
            self.nextTrack = self.problem.mhs.getNextTrack()
            self.dt2trackEnd = self.nextTrack.endTime - self.problem.time
        self.writepos()


    def getIsPrinting( self ):
        return (self.problem.mhs.currentTrack.type == TrackType.printing)

    def getNdofs( self ):
        return self.problem.ls.ndofs()


class DriverAnalytical( DriverReference ):
    def solve(self):
        self.problem.fakeIter()
        hs = self.problem.mhs
        radius = hs.radius
        rho = self.problem.material.density
        k = self.problem.material.conductivity
        cp = self.problem.material.specificHeat
        Tenv = self.problem.input["environmentTemperature"]
        speed = np.linalg.norm( hs.speed )
        kappa = k / rho / cp
        TRosenthal = np.zeros( self.problem.domain.mesh.nnodes )
        TNGuyen = np.zeros( self.problem.domain.mesh.nnodes )
        # Params NGuyen
        time = self.problem.time
        N = 3*np.sqrt(3) * hs.power / (rho*cp*np.power(np.pi, 1.5))
        taoMax = time
        taoSamples = np.linspace( 0, taoMax )
        for inode in range(self.problem.domain.mesh.nnodes):
            pos = self.problem.domain.mesh.pos[ inode, : ] - hs.position
            R = np.linalg.norm( pos )
            # Rosenthal solution
            TRosenthal[inode] = Tenv + hs.power / 4 / np.pi / k / R * \
                    np.exp( -speed*(R + pos[0]) / kappa )
            # NGuyen solution
            integrand = np.zeros( taoSamples.size )
            for idx, tao in enumerate(taoSamples):
                A1 = -3*((pos[0] - hs.speed[0]*(tao - time))**2 + pos[1]**2 + pos[2]**2)
                A1 /= (12*kappa*(time-tao) + radius**2)
                A1 = np.exp( A1 )
                integrand[idx] = A1 / np.power( (12*kappa*(time-tao) + radius**2), 1.5)
            integral = np.trapz( integrand, taoSamples )
            TNGuyen[inode] = Tenv + N * integral

        self.rosenthal = mhs.Function( self.problem.domain, TRosenthal )
        self.nguyen = mhs.Function( self.problem.domain, TNGuyen )

    def writepos(self):
        self.problem.writepos(
                functions={
                    "TRosenthal":self.rosenthal,
                    "TNGuyen":self.nguyen,
                    },
                )
