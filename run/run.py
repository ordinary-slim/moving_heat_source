import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import matplotlib.pyplot as plt
from myPlotHandler import myPlotHandler

labelsMapping = {0: "FE", 1: "BE", 2:"BDF2", 3:"BDF3", 4:"BDF4"}

if __name__=="__main__":
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    dFE = dict(d)
    dBE = dict(d)
    dBDF2 = dict(d)
    dBDF3 = dict(d)
    dBDF4 = dict(d)
    dFE["timeIntegration"] = 0
    dBE["timeIntegration"] = 1
    dBDF2["timeIntegration"] = 2
    dBDF3["timeIntegration"] = 3
    dBDF4["timeIntegration"] = 4
    #problemsDics = [dFE, dBE, dBDF2]
    problemsDics = [dBE, dBDF2, dBDF3, dBDF4]
    problems = []
    for d in problemsDics:
        p = mhs.Problem()
        p.initialize( d )
        problems.append( p )
    problems = [ (problems[i], problemsDics[i]) for i in range(len(problems)) ]

    nsteps=int(d["maxIter"])
    #show IC
    plotHandler = myPlotHandler(problems[0][0].mesh, pauseTime=2)
    plt.figure(dpi=200)
    for p, d in problems:
        plotHandler.plotProblem( p, label=labelsMapping[d["timeIntegration"]] )
    plt.pause( plotHandler.pauseTime );
    #tstepping
    for istep in range(nsteps):
        plotHandler.clf( problems[0][0].mesh )
        for p, d in problems:
            p.iterate()
            plotHandler.plotProblem( p, label=labelsMapping[d["timeIntegration"]] )
        plt.pause( plotHandler.pauseTime );
