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
import pdb

labelsMapping = {0: "FE", 1: "BE", 2:"BDF2", 3:"BDF3", 4:"BDF4"}
ordersMapping = {0: 1, 1: 1, 2:2, 3:3, 4:4}
stepsReqMapping = {0: 0, 1: 0, 2:1, 3:2, 4:3}
maxNstepsRequired = max(stepsReqMapping.values())#start off everyone @ same pos
maxNstepsRequired = maxNstepsRequired + 2

def formatFloat( fl ):
    return str(round(fl, 4)).replace(".", "_")

def computeL2Error( p1, p2 ):
    # safety check
    if (p1.time != p2.time):
        print("Something wrong!")
        exit()
    return np.sqrt( np.trapz( np.square( p1.solution - p2.solution ), p1.mesh.pos ) )

def advanceUntilTfinal( p, Tfinal, plotHandler=None ):
    promptFreq = 100
    counter = 0
    print( "Run: t0={}, dt={}, Tf={}".format( p.time, p.dt, Tfinal) )
    while( round(p.time, 4) < Tfinal ):
        counter += 1
        if counter%promptFreq==0:
            print( "Current time:", p.time )
        p.iterate()
        if plotHandler:
            plotHandler.plotProblem( p, label="Run"  );#debug
            plotHandler.pause();#debug
            plotHandler.clf( p.mesh )

def getRFL( d ):
    L = d["Right"] - d["Left"]
    nels = d["nels"]
    h = L/float(nels)
    dt = d["dt"]
    speed = d["speedX"]
    return (speed * dt) / h

def set_dt( d, RFL ):
    L = d["Right"] - d["Left"]
    r = d["radius"]
    speed = d["speedX"]
    dt = (RFL * r) / float(speed)
    d["dt"] = dt
    print( "RFL = {}, dt = {}".format( RFL, dt ) )
    return dt

def writeReferenceSolution(postFolder, referenceDT=1e-2):
    dt = referenceDT
    fileName = "l2.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    d["TimeIntegration"] = 2
    Tfinal = d["Tfinal"]
    fileName = "{}_DT{}_T{}.csv".format(
            labelsMapping[d["timeIntegration"]],
            formatFloat(dt),
            formatFloat( Tfinal ),
            )
    fileName = "reference_" + fileName
    # check for file existence
    if os.path.isfile( "./" + postFolder + "/" + fileName ):
        return

    d["dt"] = dt
    p = mhs.Problem()
    p.initialize( d )
    advanceUntilTfinal( p, Tfinal )
    pf = writePost( p, fileName, postFolder )
    print("Wrote reference solution to " + pf + ".")
    return pf


def computeDataSets(postFolder,
        checkForExistence=True,
        computeRefSolution=False):
    fileName = "l2.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    Tfinal = d["Tfinal"]
    #DTs = [100, 80, 60, 50, 40.0, 25, 20.0, 10.0, 5.0, 4.0, 2.5, 2, 1.0, 0.5, 0.1]
    DTs = [2, 1, 0.8, 0.5, 0.25, 0.2, 0.1]
    p = mhs.Problem()

    if computeRefSolution:
        refSoluton = writeReferenceSolution(postFolder, referenceDT=0.00001)

    for dt in DTs:
        d["dt"] = dt
        # prepare fine problem
        ##determine fine tstep size
        approxFine_dt = pow(dt, 2)
        # treat case dt > 1
        approxFine_dt = min( approxFine_dt, dt / 16.0 )
        fineStepsPerStep = int( np.ceil( dt / approxFine_dt ) )
        fine_dt = d["dt"] / float( fineStepsPerStep )
        # initialize fine iterator
        dFine = dict( d )
        dFine["dt"] = fine_dt
        pFine = mhs.Problem()
        pFine.initialize( dFine )
        ##do N blocks of fine time steps storing
        prevSols = [np.array(pFine.solution)]
        for istep in range( maxNstepsRequired ):
            for fineIstep in range(fineStepsPerStep):
                pFine.iterate()
            prevSols.insert( 0, np.array(pFine.solution) )
        arrPrevSols = np.transpose(np.stack( prevSols ))

        for timeIntegration in [1, 2, 3, 4]:
            d["timeIntegration"] = timeIntegration
            fileName = "{}_DT{}_T{}.csv".format(
                    labelsMapping[timeIntegration],
                    formatFloat(dt),
                    formatFloat( Tfinal ),
                    )
            if checkForExistence:
                if os.path.isfile( "./" + postFolder + "/" + fileName ):
                    print("Filename {} already exists".format(
                        fileName))
                    continue
            ##initialize time integrator for actual problem
            p = mhs.Problem()
            p.initialize( d )
            p.time = pFine.time
            p.currentPosition = pFine.mhs.currentPosition
            p.initializeIntegrator( arrPrevSols )
            advanceUntilTfinal( p, Tfinal )
            pf = writePost( p, fileName, postFolder )
            print("Wrote solution to " + pf + ".")

def retrieveDataSets(dataFolder):
    ##get it
    dataFiles = os.listdir( dataFolder )
    ##organize it
    output = {"reference" : {},
              "BE" : {},
              "BDF2" : {},
              "BDF3" : {},
              "BDF4" : {}}

    for f in dataFiles:
        fileName = dataFolder + "/" + f
        f = f.replace("_", ".")
        dt = re.search(r"DT(\d+(\.\d+)?)", f )
        if dt:
            dt = float(dt.group(1))
        else:
            print("Wrong file formatting for file {}".format(
                f))
            exit()
        p = loadProblem( fileName )
        p.dt = dt
        for key in output.keys():
            if f.startswith(key):
                print(f)
                output[key][dt] = p
    return output

def postProcess( simDict ):
    # retrieve reference solution
    referenceProblem = simDict.pop("reference", None)
    if referenceProblem:
        referenceProblem = referenceProblem[list(referenceProblem.keys())[0]]
    #process it
    L2Errors = dict( simDict )
    for TIlabel, TIdic in simDict.items():
        for dt, p in TIdic.items():
            L2Errors[TIlabel][dt] = computeL2Error( p, referenceProblem )
#
    ### plot
    colors = ["red", "blue", "magenta", "green"]
    plt.figure(dpi=200)
    for idx, (TIlabel, TIdic) in enumerate(L2Errors.items()):
        DTs = np.array( list(TIdic.keys()) )
        L2Error = np.array( list(TIdic.values()) )
        #sort
        p = np.argsort( DTs )
        DTs = DTs[p]
        L2Error = L2Error[p]
        #compute order
        order = np.log( L2Error[3] / L2Error[4] ) / np.log( DTs[3] / DTs[4] )
        cte = L2Error[4] / np.power(DTs[4], order)


        plt.loglog( DTs, L2Error,
                linestyle='None',
                color=colors[idx],
                marker="o",
                label=TIlabel+", order {}".format( round(order, 2) ) )
        plt.loglog( DTs, cte*np.power( DTs, order ),
                color=colors[idx],
                linestyle='--',
                linewidth=0.5,
                )
        plt.grid()
        plt.ylabel("L2")
        plt.xlabel("dt")
        plt.legend()
    plt.show()

if __name__=="__main__":
    #compute data
    dataFolder="data/fedR20_0T20_0"
    os.makedirs( dataFolder, exist_ok=True )
    computeDataSets(dataFolder,
            checkForExistence=True,
            computeRefSolution=False)
    ## get it
    output = retrieveDataSets(dataFolder)
    # post it
    postProcess( output )
