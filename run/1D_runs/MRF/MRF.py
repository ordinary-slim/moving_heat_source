import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
sys.path.insert(1, '../..')
import MovingHeatSource as mhs
from readInput import *
from prepos import *
import matplotlib.pyplot as plt
from myPlotHandler import myPlotHandler
import numpy as np
import os
import re
import pandas as pd

labelsMapping = {0: "FE", 1: "BE", 2:"BDF2", 3:"BDF3", 4:"BDF4"}
ordersMapping = {0: 1, 1: 1, 2:2, 3:3, 4:4}
stepsReqMapping = {0: 0, 1: 0, 2:1, 3:2, 4:3}
maxNstepsRequired = max(stepsReqMapping.values())#start off everyone @ same pos
maxNstepsRequired = maxNstepsRequired + 2

def formatFloat( fl ):
    return str(round(fl, 4)).replace(".", "_")

def getadimR( d ):
    L = d["Right"] - d["Left"]
    nels = d["nels"]
    h = L/float(nels)
    dt = d["dt"]
    speed = max(abs(d["speedX"]), abs(d["advectionSpeedX"]))
    return (speed * dt) / h

def set_dt( d, adimR ):
    L = d["Right"] - d["Left"]
    r = d["radius"]
    speed = max(abs(d["speedX"]), abs(d["advectionSpeedX"]))
    dt = (adimR * r) / float(speed)
    d["dt"] = dt
    print( "adimR = {}, dt = {}".format( adimR, dt ) )
    return dt

def main(figureFolder):
    fileName = "mrf.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    Tfinal = d["Tfinal"]
    adimR = 0.5
    dt = set_dt( d, adimR )
    if not figureFolder:
        figureFolder = "figures/adimR{}".format( formatFloat( adimR ) )
    os.makedirs( figureFolder, exist_ok=True )

    d["dt"] = dt
    # prepare fine problem
    ##determine fine tstep size
    approxFine_dt = pow(dt, 2)
    # treat case dt > 1
    approxFine_dt = min( approxFine_dt, dt / 16.0 )
    fineStepsPerStep = int( np.ceil( dt / approxFine_dt ) )
    fine_dt = d["dt"] / float( fineStepsPerStep )
    print( "Time step = {}, fine time step = {}".format( dt, fine_dt ) )

    problems = []
    #MRF
    pMRF = mhs.Problem()
    pMRF.initialize( d )
    pMRF.label = "MRF"
    problems.append( pMRF )
    #FRF
    pFRF = mhs.Problem()
    dFRF = dict(d)
    dFRF["sourceTerm"] = 0
    dFRF["speedX"] = -d["advectionSpeedX"]
    dFRF["isAdvection"] = 0
    pFRF.initialize( dFRF )
    pFRF.label = "FRF"
    problems.append( pFRF )
    #Fine
    dFine = dict( dFRF )
    dFine["dt"] = fine_dt
    dFine["timeIntegration"] = 2
    pFine = mhs.Problem()
    pFine.initialize( dFine )
    pFine.label = "Fine"

    plt.figure(dpi=200)
    plotHandler = myPlotHandler(problems[0].mesh,
                                pauseTime=0.1,
                                figureFolder=figureFolder,
                                Tmin=d["environmentTemperature"])
    plotHandler.left = -100.0
    promptFreq = 100
    counter = 0
    print( "Run: t0={}, dt={}, Tf={}".format( problems[0].time, problems[0].dt, Tfinal) )
    if plotHandler:
        plotHandler.plotProblem( pFine, label=pFine.label, linestyle='--')
        for p in problems:
            advectSolution = False
            if ("MRF" in p.label): advectSolution = True
            plotHandler.plotProblem( p, label=p.label, advectSolution=advectSolution)
        plotHandler.pause()
        plotHandler.save()
        plotHandler.clf( problems[0].mesh )

    iteration = 0
    maxIter = d["maxIter"]
    while( round(problems[0].time, 4) < Tfinal and iteration<maxIter ):
        iteration += 1
        for istep in range(fineStepsPerStep):
            pFine.iterate()
        plotHandler.plotProblem( pFine, linestyle="--", label=pFine.label, updateLims=True)

        for p in problems:
            counter += 1
            if counter%promptFreq==0:
                print( "Current time:", p.time )
            p.iterate()

            if plotHandler:
                advectSolution = False
                if p.label=="MRF": advectSolution = True
                plotHandler.plotProblem( p, label=p.label, advectSolution=advectSolution)

        if plotHandler:
            plt.title(r"$\mathcal{R}=" + str(adimR) + r"$")
            plotHandler.pause()
            plotHandler.save()
            plotHandler.clf( problems[0].mesh )

if __name__=="__main__":
    #compute data
    main("")
    #debug()
