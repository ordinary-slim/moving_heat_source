import MovingHeatSource as mhs
import numpy as np

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
    initialPosition = problemInput["initialPosition"]
    halfLength = adimR * radius
    box = [initialPosition[0] - halfLength, initialPosition[0] + halfLength,
           initialPosition[1] - halfLength, initialPosition[1] + halfLength]
    return mesh(box[0], box[1], meshDen=meshDen)

'''
import matplotlib.pyplot as plt
def plotProblem( p ):
    plt.plot( p.domain.posLab[:, 0],
             p.unknown.values)
    gammaNodesIndices = p.gammaNodes.getIndices()
    for idx in gammaNodesIndices:
        x = p.domain.posLab[idx, 0]
        plt.axvline(x=x)
'''

def setAdimR( adimR, input ):
    r = input["radius"]
    speed = max( np.linalg.norm(input["HeatSourceSpeed"]), np.linalg.norm(input["advectionSpeed"]) )
    return (adimR * r / speed)

def run():
    inputFile = "input.yaml"
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
    movingProblemInput["advectionSpeed"] = -fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["speedDomain"]      = fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["HeatSourceSpeed"] = np.zeros(3)

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
        pMoving.intersectExternal( pFixed, False )
        pFixed.preiterate(False)
        pMoving.preiterate(False)
        pMoving.intersectExternal( pFixed, False )
        pFixed.substractExternal( pMoving, True )
        pMoving.updateInterface( pFixed )

        #Dirichet gamma
        pFixed.setGamma2Dirichlet()

        # Pre-assembly, updating free dofs
        pMoving.preAssemble(False)
        pFixed.preAssemble(False)
        ls = mhs.LinearSystem.Create( pMoving, pFixed )
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

        pMoving.writepos(
            nodeMeshTags={ "gammaNodes":pMoving.gammaNodes, },
                )
        pFixed.writepos(
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },
                )

def test():
    run()
    refdatasetsfixed = ["post_fixed_reference/fixed_{}.vtu".format(i+1) for i
             in range(4)]
    newdatasetsfixed = ["post_fixed/fixed_{}.vtu".format(i+1) for i in range(4)]

    refdatasetsmoving = ["post_moving_reference/moving_{}.vtu".format(i+1) for i
             in range(4)]
    newdatasetsmoving = ["post_moving/moving_{}.vtu".format(i+1) for i in range(4)]

    # COMPARISON
    comparisonsFixed = []
    for refds, newds in zip(refdatasetsfixed, newdatasetsfixed):
        comparisonsFixed.append( mhs.meshio_comparison(refds, newds) )
    comparisonsMoving = []
    for refds, newds in zip(refdatasetsmoving, newdatasetsmoving):
        comparisonsMoving.append( mhs.meshio_comparison(refds, newds) )

    assert all(comparisonsFixed + comparisonsMoving)

if __name__=="__main__":
    test()
