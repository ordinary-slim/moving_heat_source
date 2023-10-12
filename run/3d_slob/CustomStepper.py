import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import AdaptiveStepper
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

        backRadiusObb = max(backRadius - radius, 0.0)
        p0 = self.pMoving.mhs.position - backRadiusObb*xAxis
        p1 = self.pMoving.mhs.position + self.adimMinRadius*xAxis
        obb = mhs.MyOBB( p0, p1, 2*sideRadius, 2*zRadius )
        subdomainEls = self.pMoving.domain.mesh.findCollidingElements( obb )
        collidingElsBackSphere = self.pMoving.domain.mesh.findCollidingElements( p0, self.adimMinRadius*radius )
        #collidingElsFrontSphere = self.pMoving.domain.mesh.findCollidingElements( self.pMoving.mhs.position, self.adimMinRadius*radius )
        subdomainEls += collidingElsBackSphere
        #subdomainEls += collidingElsFrontSphere
        subdomain = mhs.MeshTag( self.pMoving.domain.mesh, self.pMoving.domain.mesh.dim, subdomainEls )
        self.pMoving.domain.intersect( subdomain )


    def increaseDt( self ):
        self.adimDt = self.adimDt + 0.5
        for _ in range( self.lastVals4mean ):
            self.maximumTs.append( np.inf )

    def setBCs( self ):
        pass

    def computeSizeSubdomain( self, adimDt = None ):
        if adimDt is None:
            adimDt = self.adimDt
        #adimSubdomainSize = max(self.adimSubdomainSize, 2.75)
        #adimSubdomainSize = min(self.adimDt + 2.0, self.adimDt * 2 )
        return adimDt + 4

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
        try:
            self.update()
        except ZeroDivisionError:
            self.writepos()
            raise
        # Write vtk files
        self.writepos()

    def writepos( self ):
        activeInExternal = self.pFixed.getActiveInExternal( self.pMoving, 1e-7 )
        self.pFixed.writepos(
            shift=-self.pFixed.mhs.position,
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
            shift=-self.pFixed.mhs.position,
            nodeMeshTags={ "gammaNodes":self.pMoving.gammaNodes, },
            )




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
        self.tol = 1e-7
        self.nextTrack = None

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
        self.setDtFromAdimR( 0.5, self.dt2trackEnd )
        tnp1 = self.problem.time + self.problem.dt
        track = self.problem.mhs.path.interpolateTrack( tnp1 )
        if (track.hasDeposition) and (self.printer is not None):
            self.printer.deposit( self.problem.mhs.position,
                                self.problem.mhs.path.interpolatePosition(tnp1)
                               )
        #self.problem.setConvection( resetBcs = True )
        self.solve()
        self.dt2trackEnd = self.problem.mhs.currentTrack.endTime - self.problem.time
        if self.dt2trackEnd < self.tol:
            self.nextTrack = self.problem.mhs.getNextTrack()
            self.dt2trackEnd = self.nextTrack.endTime - self.problem.time
        self.writepos()

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
