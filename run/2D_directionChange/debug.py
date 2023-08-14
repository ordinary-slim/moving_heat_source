import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
import numpy as np
import meshzoo
from AdaptiveStepper import AdaptiveStepper
import pdb

Tfinal = 2
inputFile = "input.yaml"
tol = 1e-7
adimR_domain = 3
problemInput = mhs.readInput( inputFile )
adimR = problemInput["radius"] / problemInput["HeatSourceSpeedX"]

def meshBox(box, meshDen=2):
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

def meshAroundHS( adimR, problemInput, meshDen=2 ):
    radius = problemInput["radius"]
    initialPositionX = problemInput["initialPositionX"]
    initialPositionY = problemInput["initialPositionY"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 3*radius )
    halfLengthY = min( trailLength, capotLength )
    box = [initialPositionX - trailLength, initialPositionX + capotLength,
           initialPositionY - halfLengthY, initialPositionY + halfLengthY]
    return meshBox(box, meshDen)

class MyIterator:
    def __init__( self, pFixed, pMoving, cutoffRadius=1, isCoupled=True, isSteady=False  ):
        self.isCoupled = isCoupled
        self.adaptiveStepper = AdaptiveStepper( pFixed, pMoving, adimR_domain )

    def iterate( self, pFixed, pMoving ):
        #self.adaptiveStepper.setDt()

        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.resetActivation()
        #TODO: Replace by intersect with OBB
        #pMoving.domain.intersectObb( self.adaptiveStepper.OBB )

        pMoving.intersectExternal( pFixed, False )
        pFixed.preiterate(False)

        # Update speeds moving domain
        pMoving.setAdvectionSpeed( -pFixed.mhs.speed )
        pMoving.domain.mesh.setSpeedFRF( pFixed.mhs.speed )
        pMoving.preiterate(False)

        pMoving.intersectExternal( pFixed, False )
        if self.isCoupled:
            pFixed.substractExternal( pMoving, False )
        pFixed.updateInterface( pMoving )
        pMoving.updateInterface( pFixed )
        #Dirichet gamma
        pFixed.setGamma2Dirichlet()
        # Pre-assembly, updating free dofs
        pMoving.preAssemble(allocateLs=True)
        pFixed.preAssemble(allocateLs=True)
        ls = mhs.LinearSystem.Create( self.pMoving, self.pFixed )
        # Assembly
        pMoving.assemble()
        pFixed.assemble()
        # Assembly Gamma
        pFixed.assembleDirichletGamma( pMoving )
        pMoving.assembleNeumannGamma( pFixed )
        # Build ls
        ls.assemble()
        # Solve ls
        ls.solve()
        # Recover solution
        pFixed.gather()
        pMoving.gather()
        pFixed.unknown.interpolateInactive( pMoving.unknown, False )
        if self.isCoupled:
            pMoving.unknown.interpolateInactive( pFixed.unknown, True )
        else:
            pMoving.unknown.interpolate( pFixed.unknown )

        # Post iteration
        pFixed.postIterate()
        pMoving.postIterate()

        # Compute next dt && domain size
        #self.adaptiveStepper.update()

        pFixed.writepos(
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },)
        pMoving.writepos(
            nodeMeshTags={ "gammaNodes":pMoving.gammaNodes, },
            )

if __name__=="__main__":
    boxPhys = [-10, 10, -10, 10]

    # read input
    fixedProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshInputPhys, meshInputMoving = {}, {}
    meshInputPhys["points"], meshInputPhys["cells"], meshInputPhys["cell_type"] = meshBox(boxPhys)
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput)

    meshFixed           = mhs.Mesh(meshInputPhys)
    meshMoving           = mhs.Mesh(meshInputMoving)

    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeedX"] = -fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["speedFRF_X"]      = fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["HeatSourceSpeedX"] = 0.0

    pFixed   = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pMoving  = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")


    maxIter = pFixed.input["maxIter"]

    pFixed.setDt( 1*adimR )
    pMoving.setDt( 1*adimR )

    # Set path
    pFixed.mhs.setPath( *gcode2laserPath( "Path.gcode" ) )
    myIterator = MyIterator( pFixed, pMoving )

    while not(pFixed.mhs.path.isOver):
        myIterator.iterate( pFixed, pMoving )
