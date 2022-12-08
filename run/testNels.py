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
    nels = [2, 4, 8]
    problems = []
    for n in nels:
        d["nels"] = n
        p = mhs.Problem()
        p.initialize( d )
        problems.append( p )

    for p in problems:
        print("nels = {}, time = {}".format( p.mesh.nels, p.time ) )
        p.iterate()
        print("=================")
