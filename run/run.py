import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import matplotlib.pyplot as plt
import numpy as np

class plotter:
    def __init__(self, m):
        self.left = min(m.pos)
        self.right = max(m.pos)
        self.Tmin = -1
        self.Tmax = 1
        self.range = 0
        self.figCleared = False

    def clf(self, mesh):
        plt.clf()
        self.figCleared = True
        plt.plot(mesh.pos, np.zeros( len(mesh.pos) ), '-o')

    def plotProblem( self, p, label="solution" ):
        plt.plot( p.mesh.pos, p.solution, label=label );
        plt.xlim( self.left, self.right );
        #get max min temperatures
        Tmin = min( p.solution )
        Tmax = max( p.solution )
        if (self.Tmin > Tmin or self.Tmax < Tmax):
            self.Tmin = Tmin
            self.Tmax = Tmax
            self.range = self.Tmax - self.Tmin
        plt.ylim( self.Tmin - self.range*0.1, self.Tmax + self.range*0.1)

        plt.annotate("time = " + str(round(p.time, 2)), (0.05, 0.95));
        plt.legend();

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
    pFE = mhs.Problem()
    pFE.initialize( dFE )
    pBE = mhs.Problem()
    pBE.initialize( dBE )
    pBDFS2 = mhs.Problem()
    pBDFS2.initialize( dBDFS2 )

    nsteps=int(d["maxIter"])

    #show IC
    plotHandler = plotter(pFE.mesh)
    plt.figure(dpi=200)
    #plotHandler.plotProblem( pFE, label="FE" )
    plotHandler.plotProblem( pBE, label="BE" )
    plotHandler.plotProblem( pBE, label="BDFS2" )
    plt.pause( 0.25 );
    #tstepping
    for istep in range(nsteps):
        pFE.iterate()
        pBE.iterate()
        pBDFS2.iterate()
        plotHandler.clf( pFE.mesh )
        #plotHandler.plotProblem( pFE, label="FE" )
        plotHandler.plotProblem( pBE, label="BE" )
        plotHandler.plotProblem( pBDFS2, label="BDFS2" )
        plt.pause( 0.25 );
