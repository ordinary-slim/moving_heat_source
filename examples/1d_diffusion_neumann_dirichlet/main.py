import MovingHeatSource as mhs
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse

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

def plotProblems( pLeft, pRight ):
    # Default settings
    mpl.rcParams.update(mpl.rcParamsDefault)

    plt.style.use("bmh")
    mpl.rcParams.update({'font.size':22})

    #sns.set_palette( sns.color_palette("bright") )
    plt.plot( pLeft.domain.mesh.posFRF[:, 0],
             pLeft.unknown.values,
             label="Left subproblem")
    plt.plot( pRight.domain.mesh.posFRF[:, 0],
             pRight.unknown.values,
             label="Right subproblem")
    gammaNodesIndices = pLeft.gammaNodes.getIndices()

    xGamma = pLeft.domain.mesh.posFRF[gammaNodesIndices[0], 0]

    gammaSettings = {"linestyle" : "--", "linewidth" : 1, "color" : "k"}
    gammaSettings["label"] = "Gamma"
    plt.axvline(x=xGamma, **gammaSettings)

    plt.gca().set_xlim( [pLeft.domain.mesh.posFRF[0, 0],
                         pRight.domain.mesh.posFRF[-1,0],] )
    #plt.gca().set_ylim( [0, 10] )
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    #plt.pause(0.5)


def main(plotting):
    if plotting:
        fig = plt.figure( figsize=[ 12, 8 ] )

    inputFile = "input.yaml"
    boxLeft = [-25, 0]
    boxRight = [0, 25]

    # Read input
    problemInput = mhs.readInput( inputFile )

    # Mesh
    leftMeshInput, rightMeshInput = {}, {}
    meshDen = 1
    leftMeshInput["points"], leftMeshInput["cells"], leftMeshInput["cell_type"] = mesh(*boxLeft, meshDen=meshDen)
    rightMeshInput["points"], rightMeshInput["cells"], rightMeshInput["cell_type"] = mesh(*boxRight, meshDen=meshDen)

    # open integration
    leftMeshInput["numberOfGaussPointsCells"] =  3
    rightMeshInput["numberOfGaussPointsCells"] = 3

    meshLeft = mhs.Mesh(leftMeshInput)
    meshRight = mhs.Mesh(rightMeshInput)

    # Initialize problems
    pLeft  = mhs.Problem(meshLeft, problemInput, caseName="left")
    #problemInput["advectionSpeed"] = np.array([-10, 0.0, 0.0])
    pRight  = mhs.Problem(meshRight, problemInput, caseName="right")

    f = lambda pos : max( 0.0, 10.0 - abs( pos[0] ) )
    for p in [pLeft, pRight]:
        p.forceState( f )

    tol = 1e-7
    Tfinal = problemInput["Tfinal"]
    while (pLeft.time < Tfinal - tol):
        print( "iter# {}, time={}".format(
            pLeft.iter,
            pLeft.time) )


        pLeft.preIterate( canPreassemble=False )
        pRight.preIterate( canPreassemble=False )

        pLeft.updateInterface( pRight )
        pRight.updateInterface( pLeft )

        # Dirichlet interface left
        pRight.setGamma2Dirichlet()
        pLeft.setGamma2Neumann()

        # Pre-assembly, updating free dofs
        pLeft.preAssemble(allocateLs=True)
        pRight.preAssemble(allocateLs=True)
        ls = mhs.LinearSystem.Create( pLeft, pRight )

        pLeft.assemble( pRight )
        pRight.assemble( pLeft )


        ls.assemble()

        ls.solve()

        pRight.gather()
        pLeft.gather()

        pRight.postIterate()
        pLeft.postIterate()

        pRight.writepos()
        pLeft.writepos()
        if plotting:
            plt.clf()
            plotProblems( pLeft, pRight )

    if plotting:
        plt.show()

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    plotting = args.plot
    main(plotting)
