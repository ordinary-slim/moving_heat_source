import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
from readInput import *
from prepos import *
import matplotlib.pyplot as plt
from myPlotHandler import myPlotHandler
import numpy as np
import os
import re
import pandas as pd
import pdb

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
    speed = d["speedX"]
    return (speed * dt) / h

def set_dt( d, adimR ):
    L = d["Right"] - d["Left"]
    r = d["radius"]
    speed = d["speedX"]
    dt = (adimR * r) / float(speed)
    d["dt"] = dt
    print( "adimR = {}, dt = {}".format( adimR, dt ) )
    return dt

def printFigures(figureFolder):
    fileName = "input.txt"
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
    # initialize fine iterator
    dFine = dict( d )
    dFine["dt"] = fine_dt
    dFine["timeIntegration"] = 2
    pFine = mhs.Problem()
    pFine.initialize( dFine )
    pFine.label = "Fine"
    ##do N blocks of fine time steps storing
    prevSols = [np.array(pFine.solution)]
    timesPrevSols  = [pFine.time]
    for istep in range( maxNstepsRequired ) :
        for fineIstep in range(fineStepsPerStep):
            pFine.iterate()
        prevSols.insert( 0, np.array(pFine.solution) )
        timesPrevSols.insert( 0, pFine.time )
    arrPrevSols = np.transpose(np.stack( prevSols ))

    problems = []
    for timeIntegration in [1]:
        d["timeIntegration"] = timeIntegration
        ##initialize time integrator for actual problem
        p = mhs.Problem()
        p.initialize( d )
        p.timeIntegration = timeIntegration
        p.label = "FRF, " + labelsMapping[timeIntegration]
        p.time = pFine.time
        p.mhs.currentPosition = pFine.mhs.currentPosition
        p.initializeIntegrator( arrPrevSols )
        problems.append( p )

    #MRF
    pMRF = mhs.Problem()
    dMRF = dict(d)
    dMRF["isAdvection"] = 1
    dMRF["advectionSpeedX"] = -d["speedX"]
    dMRF["speedX"] = 0.0
    dMRF["timeIntegration"] = 1
    pMRF.initialize( dMRF )
    pMRF.time = pFine.time
    pMRF.mhs.currentPosition = pFine.mhs.currentPosition
    pMRF.initializeIntegrator( arrPrevSols )
    pMRF.label = "MRF, " + labelsMapping[d["timeIntegration"]]
    problems.append( pMRF )

    plt.figure(dpi=200)
    plotHandler = myPlotHandler(problems[0].mesh,
                                pauseTime=0.1,
                                figureFolder=figureFolder,
                                Tmin=d["environmentTemperature"])

    plt.title(r"$\mathcal{R}=" + str(adimR) + r"$")
    promptFreq = 100
    counter = 0
    print( "Run: t0={}, dt={}, Tf={}".format( problems[0].time, problems[0].dt, Tfinal) )
    if plotHandler:
        plotHandler.plotProblem( pFine, label=pFine.label, linestyle='--');#debug
        for p in problems:
            plotHandler.plotProblem( p, label=p.label);#debug
        plotHandler.pause()
        plotHandler.save()
        plotHandler.clf( problems[0].mesh )

    while( round(problems[0].time, 4) < Tfinal ):
        for istep in range(fineStepsPerStep):
            pFine.iterate()
        plotHandler.plotProblem( pFine, linestyle="--", label=pFine.label, updateLims=True);#debug

        for p in problems:
            counter += 1
            if counter%promptFreq==0:
                print( "Current time:", p.time )
            p.iterate()

            if plotHandler:
                plotHandler.plotProblem( p, label=p.label, updateLims=True);#debug


        if plotHandler:
            plt.title(r"$\mathcal{R}=" + str(adimR) + r"$")
            plotHandler.pause()
            plotHandler.save()
            plotHandler.clf( problems[0].mesh )


if __name__=="__main__":
    #compute data
    printFigures("")
