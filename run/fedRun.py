import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import matplotlib.pyplot as plt
from myPlotHandler import myPlotHandler
from prepos import *

labelsMapping = {0: "FE", 1: "BE", 2:"BDF2", 3:"BDF3", 4:"BDF4"}

if __name__=="__main__":
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    dFine = dict(d)
    dFE = dict(d)
    dBE = dict(d)
    dBDF2 = dict(d)
    dBDF3 = dict(d)
    dBDF4 = dict(d)
    dFine["timeIntegration"] = 1
    fineStepsPerStep = 32
    dFine["dt"] = dFine["dt"] / float( fineStepsPerStep )
    dFE["timeIntegration"] = 0
    dBE["timeIntegration"] = 1
    dBDF2["timeIntegration"] = 2
    dBDF3["timeIntegration"] = 3
    dBDF4["timeIntegration"] = 4
    #problemsDics = [dFE, dBE, dBDF2]
    problemsDics = [dBE, dBDF2, dBDF3, dBDF4]
    problems = []

    pFine = mhs.Problem()
    pFine.initialize( dFine )

    #get maxNstepsRequired
    maxNstepsRequired = 4
    #do a few fine tsteps
    postFiles = []
    times = []
    for istep in range( maxNstepsRequired ):
        for fineIstep in range(fineStepsPerStep):
            pFine.iterate()
        postFiles.append( writePost( pFine, postFolder="postFine" ))
        times.append( pFine.time )

    #feed each time integrator and set current position
    for d in problemsDics:
        d["currentPositionX"] = pFine.mhs.currentPosition[0]
        p = mhs.Problem()
        p.initialize( d )
        problems.append( p )
    problems = [ (problems[i], problemsDics[i]) for i in range(len(problems)) ]
    for p, _ in problems:
        p.currentPosition = pFine.mhs.currentPosition
        p.time            = pFine.time;
        initializeTimeIntegrator( p, postFiles )

    #show IC
    plotHandler = myPlotHandler(problems[0][0].mesh,
                                pauseTime=0.25,
                                shortDescription="fedRun")
    plt.figure(dpi=200)
    plotHandler.plotProblem( pFine, linestyle="--",  label="Fine dt" )
    for p, d in problems:
        plotHandler.plotProblem( p, label=labelsMapping[d["timeIntegration"]] )
    plotHandler.pause()
    plotHandler.save()

    nsteps=int(d["maxIter"])
    #tstepping
    for istep in range(nsteps):
        plotHandler.clf( problems[0][0].mesh )
        for fineIstep in range(fineStepsPerStep):
            pFine.iterate()
        plotHandler.plotProblem( pFine, linestyle="--",  label="Fine dt" )
        for p, d in problems:
            p.iterate()
            plotHandler.plotProblem( p, label=labelsMapping[d["timeIntegration"]] )
        plotHandler.pause()
        plotHandler.save()
