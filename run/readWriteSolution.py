import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import numpy as np
import pandas as pd
import os
import re

def writePost( p, fileName="", postFolder="" ):
    if not(fileName):
        fileName = "post{}.csv".format( str(round(p.time, 4)).replace(".", "_") )
    if not(postFolder):
        postFolder = "post"
        os.makedirs( postFolder, exist_ok=True )
    postFile = postFolder + "/" + fileName

    df = pd.DataFrame( {"X": p.mesh.pos, "T": p.solution}  )
    df.to_csv( postFile,
            index=False)
    return postFile

def initializeTimeIntegrator( p, postFiles ):
    prevSols = []
    counter = 1
    # sort post files by time
    times = []
    for postFile in postFiles:
        # get time
        t = re.search(r"(\d+_\d+)", postFile)
        if t:
            t = t.group(0).replace("_", ".")
        else:
            exit("Error!")
        times.append(float(t))
    # sort time list and retrieve permutation
    timesAndPermutation = [ (times[i], i) for i in range(len(times)) ]
    timesAndPermutation.sort(reverse=True)
    sortedTimes, permutation = zip(*timesAndPermutation)
    # sort postFiles list
    postFiles = [postFiles[i] for i in permutation]
    # time contained either in dataclass or postFileName or header
    for postFile in postFiles:
        prevSols.append( pd.read_csv( postFile, dtype=np.float64 )["T"].to_numpy() )
        counter += 1
    # build np array
    prevSols = np.array( prevSols )
    p.initializeIntegrator( np.transpose(prevSols) )

def test():
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    d["timeIntegration"] = 2
    p = mhs.Problem()
    p.initialize( d )

    # write postfiles
    postFiles = []
    postFiles.append( writePost( p ) )
    for i in range(3):
        p.iterate()
        postFiles.append( writePost( p ) )

    # read postfiles
    initializeTimeIntegrator( p, postFiles )


if __name__=="__main__":
    test()
