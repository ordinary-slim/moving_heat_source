import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
import numpy as np
import pandas as pd
import os
import shutil

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
    for postFile in postFiles:
        prevSols.append( pd.read_csv( postFile )["T"].to_numpy() )
        counter += 1
    # build np array
    prevSols = np.array( prevSols )
    p.initializeIntegrator( np.transpose(prevSols) )

def test():
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    p = mhs.Problem()
    p.initialize( d )

    # write postfiles
    postFile1 = writePost( p )
    for i in range(20):
        p.iterate()
    postFile2 = writePost( p )

    # read postfiles
    initializeTimeIntegrator( p, [postFile1, postFile2] )


if __name__=="__main__":
    test()
