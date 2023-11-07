import MovingHeatSource as mhs
from meshing import getMeshPhysical, fineElSize
from CustomStepper import CustomStepper, DriverReference
from MovingHeatSource.adaptiveStepper import LpbfAdaptiveStepper, deactivateBelowSurface
from scanningPath import writeGcode
from MyLogger import MyLogger
import pickle
import argparse
from line_profiler import LineProfiler

inputFile = "input.yaml"
problemInput = mhs.readInput( inputFile )
gcodeFile = problemInput["path"]
tol = 1e-7
fineElSizeMoving = fineElSize/problemInput["fineElFactorMovingSubdomain"]

# read input
problemInput = mhs.readInput( inputFile )

def runReference(caseName="reference"):
    mesh = getMeshPhysical()
    pReference = mhs.Problem( mesh, problemInput, caseName=caseName)
    deactivateBelowSurface( pReference ) 

    driver = DriverReference( pReference )

    logger = MyLogger()

    while not(driver.problem.mhs.path.isOver(pReference.time)):
        logger.iterate( driver )

    with open("{}.log".format(caseName), "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)

def runCoupled(caseName="fixed"):
    fixedProblemInput = dict( problemInput )
    meshFixed = getMeshPhysical()
    pFixed = mhs.Problem(meshFixed, fixedProblemInput, caseName=caseName)
    deactivateBelowSurface( pFixed )

    driver = CustomStepper( pFixed,
                           adimFineDt=0.5 / problemInput["fineTStepFactorMoving"],
                           maxAdimtDt=5,
                           elementSize=fineElSizeMoving,
                           threshold=0.15,
                           adimMinRadius=1.8,
                           slowAdimDt=1.0,
                           adimPosZLen=0.5,
                           )
    
    logger = MyLogger()
    while not(driver.pFixed.mhs.path.isOver( driver.getTime() ) ) :
        logger.iterate( driver )

    with open("{}.log".format(caseName), "wb") as reflog:
        pickle.dump( logger, reflog, pickle.HIGHEST_PROTOCOL)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run-reference', action='store_true')
    parser.add_argument('-c', '--run-coupled', action='store_true')
    parser.add_argument('-cp', '--run-coupled-profiling', action='store_true')
    parser.add_argument('--layers', default=-1, type=int)
    parser.add_argument('--case-name', default='case')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    writeGcode( nLayers=args.layers )
    if args.run_reference:
        runReference(caseName=args.case_name)
    if args.run_coupled:
        runCoupled(caseName=args.case_name)
    if args.run_coupled_profiling:
        lp = LineProfiler()
        lp.add_module(mhs)
        lp.add_module(CustomStepper)
        lp.add_module(LpbfAdaptiveStepper)
        lp_wrapper = lp(runCoupled)
        lp_wrapper(caseName=args.case_name)
        lp.print_stats()
