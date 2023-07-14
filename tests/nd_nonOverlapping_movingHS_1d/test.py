import MovingHeatSource as mhs
import numpy as np
import matplotlib.pyplot as plt

def mesh(leftEnd, rightEnd, meshDen=4):
    cell_type="line2"
    nels = int((rightEnd-leftEnd)*meshDen)
    points = np.linspace( leftEnd, rightEnd, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    auxCells = []
    for iel in range(nels):
        auxCells.append( [iel, iel+1] )
    cells = np.array( auxCells, dtype=int )
    return points, cells, cell_type

def meshAroundHS( adimR, problemInput, meshDen=4 ):
    radius = problemInput["radius"]
    initialPositionX = problemInput["initialPositionX"]
    initialPositionY = problemInput["initialPositionY"]
    halfLength = adimR * radius
    box = [initialPositionX - halfLength, initialPositionX + halfLength,
           initialPositionY - halfLength, initialPositionY + halfLength]
    return mesh(box[0], box[1], meshDen=meshDen)

def plotProblem( p ):
    plt.plot( p.domain.mesh.posFRF[:, 0],
             p.unknown.values)
    gammaNodesIndices = p.gammaNodes.getIndices()
    for idx in gammaNodesIndices:
        x = p.domain.mesh.posFRF[idx, 0]
        plt.axvline(x=x)

def setAdimR( adimR, input ):
    r = input["radius"]
    HeatSourceSpeedX = max( abs(input["HeatSourceSpeedX"]), abs(input["advectionSpeedX"]))
    HeatSourceSpeedY = max( abs(input["HeatSourceSpeedY"]), abs(input["advectionSpeedY"]))
    HeatSourceSpeedZ = max( abs(input["HeatSourceSpeedZ"]), abs(input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [HeatSourceSpeedX, HeatSourceSpeedY, HeatSourceSpeedZ] ) )
    return (adimR * r / speed)

def run():
    inputFile = "input.py"
    intervalDomain = [-25, 25]
    adimR_tstep = 10
    adimR_domain = 11

    # read input
    problemInput = mhs.readInput( inputFile )

    fixedProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshInputFixed, meshInputMoving = {}, {}
    meshDen = 1
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(intervalDomain[0], intervalDomain[1], meshDen=meshDen)
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput, meshDen=meshDen)

    meshFixed  = mhs.Mesh(meshInputFixed)
    meshMoving = mhs.Mesh(meshInputMoving)

    # mhs.Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput, movingProblemInput,]:
        input["dt"] = dt

    #set MRF business NO TRANSPORT
    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeedX"] = -fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["speedFRF_X"]      = fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["HeatSourceSpeedX"] = 0.0

    pFixed         = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pMoving        = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")

    maxIter = pFixed.input["maxIter"]
    Tfinal = pFixed.input["Tfinal"]

    while (pMoving.time < Tfinal - 1e-7) :

        pFixed.setAssembling2External( True )
        pMoving.setAssembling2External( True )

        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.resetActivation()
        pMoving.intersectExternal( pFixed, False, False )
        pFixed.preiterate(False)
        pMoving.preiterate(False)
        pMoving.intersectExternal( pFixed, False, False )
        pFixed.substractExternal( pMoving, False, True )
        pMoving.updateInterface( pFixed )

        #Dirichet gamma
        pFixed.setGamma2Dirichlet()

        # Pre-assembly, updating free dofs
        pMoving.preAssemble(True)
        pFixed.preAssemble(True)
        ls = mhs.LinearSystem( pMoving, pFixed )
        ls.cleanup()

        pMoving.assemble()
        pFixed.assemble()

        pFixed.assembleDirichletGamma( pMoving )
        pMoving.assembleNeumannGamma( pFixed )

        ls.assemble()

        ls.solve()

        pFixed.gather()
        pMoving.gather()

        # Get inactive points information from other
        pFixed.unknown.interpolateInactive( pMoving.unknown, False )
        pMoving.unknown.interpolateInactive( pFixed.unknown, True )

        pFixed.postIterate()
        pMoving.postIterate()

    pFixed.writepos(
        nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },
            )

def test():
    run()
    refds = "post_fixed_reference/fixed_4.vtu"
    newds = "post_fixed/fixed_4.vtu"

    # COMPARISON
    assert mhs.meshio_comparison(refds, newds)

if __name__=="__main__":
    test()
