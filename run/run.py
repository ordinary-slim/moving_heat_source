import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import matplotlib.pyplot as plt

class plotter:
    def __init__(self):
        self.Tmin = 1
        self.Tmax = -1

    def plot_problem( self, p ):
        plt.clf();
        plt.plot( p.mesh.pos, p.solution, color="red", label="sol");
        plt.xlim( 0.0, 10 );
        #get max min temperatures
        Tmin = min( p.solution )
        Tmax = max( p.solution )
        if (self.Tmin > Tmin or self.Tmax < Tmax):
            self.Tmin = Tmin
            self.Tmax = Tmax
        plt.ylim( self.Tmin, self.Tmax )

        plt.annotate("time = " + str(round(p.time, 2)), (0.05, 0.95));
        plt.legend();
        plt.pause( 0.05 );

if __name__=="__main__":
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    p = mhs.Problem()
    p.initialize( d )

    nsteps=100

    #show IC
    plotHandler = plotter()
    plotHandler.plot_problem( p )
    #tstepping
    for istep in range(nsteps):
        p.iterate()
        plotHandler.plot_problem( p )
