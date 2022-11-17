import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import matplotlib.pyplot as plt

if __name__=="__main__":
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    p = mhs.Problem()
    p.initialize( d )

    nsteps=5
    for istep in range(nsteps):
        p.iterate()

        plt.clf();
        plt.plot( p.mesh.pos, p.solution, color="red", label="sol");
        plt.xlim( 0.0, 10 );
        plt.ylim( -1.0, 50.0 );

        #plt.annotate("t = " + to_string( t ), 0.05, 0.95);
        plt.legend();
        plt.pause( 2 );
