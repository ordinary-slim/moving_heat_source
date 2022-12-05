import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
import MovingHeatSource as mhs
from readInput import *
from myPlotHandler import myPlotHandler
import numpy as np
import pandas as pd
import os
import re

def writePost( p, fileName="", postFolder="post" ):
    if not(fileName):
        fileName = "post{}.csv".format( str(round(p.time, 4)).replace(".", "_") )

    os.makedirs( postFolder, exist_ok=True )
    postFile = postFolder + "/" + fileName

    df = pd.DataFrame( {"X": p.mesh.pos, "T": p.solution}  )
    df.to_csv( postFile,
            index=False)
    return postFile

def getTimeFromFilename(postFile):
    # get time
    t = re.search(r"T(\d+_\d+)", postFile)
    if t:
        t = t.group(1).replace("_", ".")
    else:
        exit("Error!")
    return float(t)

def initializeTimeIntegrator( p, postFiles ):
    prevSols = []
    counter = 1
    # sort post files by time
    postFiles, times = sortPostFiles( postFiles, reverse=True )
    # time contained either in dataclass or postFileName or header
    for postFile in postFiles:
        prevSols.append( pd.read_csv( postFile, dtype=np.float64 )["T"].to_numpy() )
        counter += 1
    # build np array
    prevSols = np.array( prevSols )
    p.initializeIntegrator( np.transpose(prevSols) )

def sortPostFiles(postFiles,reverse=False):
    times = []
    for postFile in postFiles:
        # get time
        times.append(getTimeFromFilename(postFile))
    # sort time list and retrieve permutation
    timesAndPermutation = [ (times[i], i) for i in range(len(times)) ]
    timesAndPermutation.sort(reverse=reverse)
    sortedTimes, permutation = zip(*timesAndPermutation)
    # sort postFiles list
    postFiles = [postFiles[i] for i in permutation]
    return postFiles, sortedTimes

def plot1DPostFolder(postFolder="post"):
    import matplotlib.pyplot as plt
    times = []
    postFiles = list(os.listdir(postFolder))
    postFiles, times = sortPostFiles( postFiles )
    problems = []
    for pf in postFiles:
        problems.append( loadProblem(pf) )

    #grab mesh from first df
    mesh = p.mesh

    # get range for figure
    Tmin = +1e12
    Tmax = -1e12
    for df in dfs:
        tmp = max(df["T"])
        if tmp > Tmax:
            Tmax = tmp
        tmp = min(df["T"])
        if tmp < Tmin:
            Tmin = tmp

    plotHandler = myPlotHandler(mesh, Tmin=Tmin, Tmax=Tmax)
    plt.figure(dpi=200)
    for p, t in zip(problems, times):
        # plot
        plotHandler.clf( mesh )
        plotHandler.plotProblem( p, label=postFolder,
                updateLims=False,
                plotMhs=False
                )
        plt.pause( 0.25 );

def loadProblem( postFile ):
    class Bunch:
        def __init__(self, **kwds):
            self.__dict__.update(kwds)
    time = getTimeFromFilename(postFile)
    # build problem
    df = pd.read_csv( postFile )
    mesh = Bunch( pos=df["X"] )
    p = Bunch( mesh=mesh, solution=df["T"], time=time )
    return p

def test():
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    d["timeIntegration"] = 3
    d["nels"] = 5
    d["power"] = 0.0
    p = mhs.Problem()
    p.initialize( d )

    # write postfiles
    postFiles = []
    postFiles.append( writePost( p ) )

    p.iterate()
    p.solution = 26*np.ones( p.solution.size )
    postFiles.append( writePost( p ) )
    p.iterate()
    p.solution = 27*np.ones( p.solution.size )
    postFiles.append( writePost( p ) )
    p.iterate()
    p.solution = 28*np.ones( p.solution.size )
    postFiles.append( writePost( p ) )

    # read postfiles
    initializeTimeIntegrator( p, postFiles )
    p.iterate()
    print(p.solution)

if __name__=="__main__":
    test()
    #plot1DPostFolder("postFineReference")
