import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import numpy as np
from prepos import writePost, initializeTimeIntegrator

def runWriteFine():
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    d["timeIntegration"] = 2
    problemsDics = [d]
    problems = []
    for d in problemsDics:
        p = mhs.Problem()
        p.initialize( d )
        problems.append( p )
    problems = [ (problems[i], problemsDics[i]) for i in range(len(problems)) ]

    nsteps=int(d["maxIter"])

    #tstepping
    writePost(p, postFolder="postFineReference")
    for istep in range(nsteps):
        for p, d in problems:
            p.iterate()
            writePost(p, postFolder="postFineReference")

if __name__=="__main__":
    runWriteFine()
