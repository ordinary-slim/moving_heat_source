import MovingHeatSource as mhs
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from CustomStepper import DriverReference, CustomStepper
import argparse

inputFile = "input.yaml"
problemInput = mhs.readInput( inputFile )
startPoint = -problemInput["L"]/4
endPoint = +problemInput["L"]/4

def getMesh(nels, L):
    cell_type="line2"
    points = np.linspace( -L/2, +L/2, nels+1, dtype=np.float64)
    points = points.reshape( (nels+1, 1) )
    cells = np.transpose( np.vstack( (np.arange(0, nels), np.arange(1, nels+1) ), dtype=int ) )
    return mhs.Mesh( {
                    "points" : points,
                    "cells" : cells,
                    "cell_type" : cell_type,
                    "dimension" : 1,
                    })

def writeGcode():
    speed = np.linalg.norm( problemInput["HeatSourceSpeed"] )
    lines = []
    lines.append(
            "G0 F{} X{} Y0.0 Z0.0".format( speed, startPoint )
            )
    lines.append(
            "G1 X{} E0.1".format( endPoint )
            )
    with open(problemInput["path"], 'w') as f:
        f.writelines( [l+"\n" for l in lines] )


def plotProblem( p ):
    # Default settings
    mpl.rcParams.update(mpl.rcParamsDefault)

    plt.style.use("bmh")
    mpl.rcParams.update({'font.size':22})

    #sns.set_palette( sns.color_palette("bright") )
    plt.plot( p.domain.posLab[:, 0],
             p.unknown.values,
             #label="T",
             )

    plt.gca().set_xlim( [p.domain.posLab[0, 0],
                         p.domain.posLab[-1,0],] )
    plt.xlabel("x")
    plt.ylabel("u")
    #plt.legend()
    plt.pause(0.2)

def set2powder(p, xInterface = None ):
    if xInterface is None:
        xInterface = startPoint
    nels = p.domain.mesh.nels
    matSets = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[0] >= xInterface ):
            matSets.append( ielem )
    p.domain.setMaterialSets( mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, matSets ) )

def runReference(caseName="frf", plot=False):
    mesh = getMesh(nels=problemInput["nels"], L = problemInput["L"])

    # Initialize problems
    p  = mhs.Problem(mesh, problemInput, caseName=caseName)
    set2powder( p, xInterface=0.0 )
    if problemInput["print"]:
        set2powder( p )
    driver = DriverReference( p, problemInput["adimDtRef"] )

    while not(driver.problem.mhs.path.isOver( driver.problem.time ) ) :
        driver.iterate()
        driver.writepos()
        if plot:
            plt.clf()
            plotProblem( driver.problem )
    if plot:
        plt.show()


def runCoupled(caseName="coupled", plot=False):
    fixedProblemInput = dict( problemInput )

    # Mesh
    mesh = getMesh(nels=problemInput["nels"], L = problemInput["L"])

    # mhs.Problem params

    pFixed = mhs.Problem(mesh, fixedProblemInput, caseName=caseName)
    if fixedProblemInput["print"]:
        set2powder( pFixed )

    adimFineDt = 2
    adimMaxDt  = 2
    fineElSizeMoving = float(problemInput["L" ]) / problemInput["nels"]
    driver = CustomStepper( pFixed,
                           adimFineDt=adimFineDt,
                           maxAdimtDt=adimMaxDt,
                           elementSize=fineElSizeMoving,
                           threshold=0.3,
                           adimMinRadius=2,
                           #shift=np.array([0.0, 0.1, 0.0]),
                           #adimZRadius=1.0,
                           )
    
    while not(driver.pFixed.mhs.path.isOver( driver.getTime() ) ) :
        driver.iterate()
        if plot:
            plt.clf()
            plotProblem( driver.pFixed )
    if plot:
        plt.show()


if __name__=="__main__":
    writeGcode()
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run-reference', action='store_true')
    parser.add_argument('-c', '--run-coupled', action='store_true')
    parser.add_argument('--case-name', default='case')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.run_reference:
        runReference(caseName=args.case_name, plot=args.plot)
    if args.run_coupled:
        runCoupled(caseName=args.case_name, plot=args.plot)
