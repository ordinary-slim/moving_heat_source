import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import matplotlib.pyplot as plt
from myPlotHandler import myPlotHandler

labelsMapping = {0: "FE", 1: "BE", 2:"BDFS2"}

if __name__=="__main__":
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    dFE = dict(d)
    dBE = dict(d)
    dBDFS2 = dict(d)
    dFE["timeIntegration"] = 0
    dBE["timeIntegration"] = 1
    dBDFS2["timeIntegration"] = 2
    #problemsDics = [dFE, dBE, dBDFS2]
    problemsDics = [dBE, dBDFS2]
    problems = []
    for d in problemsDics:
        p = mhs.Problem()
        p.initialize( d )
        problems.append( p )
    problems = [ (problems[i], problemsDics[i]) for i in range(len(problems)) ]

    nsteps=int(d["maxIter"])
    #show IC
    plotHandler = myPlotHandler(problems[0][0].mesh)
    plt.figure(dpi=200)
    for p, d in problems:
        plotHandler.plotProblem( p, label=labelsMapping[d["timeIntegration"]] )
    plt.pause( 0.5 );
    #tstepping
    for istep in range(nsteps):
        plotHandler.clf( problems[0][0].mesh )
        for p, d in problems:
            p.iterate()
            plotHandler.plotProblem( p, label=labelsMapping[d["timeIntegration"]] )
        plt.pause( 0.5 );
