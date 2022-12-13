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

def formatFloat( fl ):
    return str(round(fl, 4)).replace(".", "_")

def computeAnalyticalSolution():
    fileName = "ode.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    T_0 = d["environmentTemperature"]
    cte = d["power"] / d["rho"] / d["specific_heat"]
    return (lambda t : T_0 + cte*(1 - np.exp(-t)))

def computeErrorAtPoint( p, f ):
    return abs( p.solution[0] - f(p.time) )
computeAnalyticalSolution
def advanceODE_untilTfinal( p, Tfinal):
    tol = 1e-4
    time = [p.time]
    temperature = [p.solution[0]]
    while( Tfinal - p.time > tol ):
        p.iterate()
        time.append( p.time )
        temperature.append( p.solution[0] )
    return time, temperature

def set_dt( d, nSteps ):
    Tf = d["Tfinal"]
    d["dt"] = Tf / float(nSteps)
    return d

def computeDataSets(postFolder, computeRefSolution=False):
    fileName = "ode.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    Tfinal = d["Tfinal"]
    time = 0.0
    nsteps = [20]
    p = mhs.Problem()

    plt.figure(dpi=200)
    times = []
    f = computeAnalyticalSolution()
    for timeIntegration in [1, 2, 3, 4]:
        d["timeIntegration"] = timeIntegration
        for N in nsteps:
            d = set_dt( d, N )
            dt = d["dt"]
            nnodes = int(d["nels"]) + 1
            maxStepBacks = 4
            prevSols = np.ones((nnodes, maxStepBacks))
            for stepBack in range(maxStepBacks):
                prevSols[:, stepBack] = f(time-dt*stepBack)*np.ones(nnodes)

            fileName = "{}_nsteps{}.csv".format(
                    labelsMapping[timeIntegration],
                    N,
                    )
            #if os.path.isfile(fullPath):
                #print("Filename {} already exists".format(
                    #fileName))
                #continue
            p = mhs.Problem()
            p.initialize( d )
            p.initializeIntegrator( prevSols )
            p.time = time
            times, temperature = advanceODE_untilTfinal( p, Tfinal)
            plt.scatter(times, temperature, label=labelsMapping[timeIntegration])
            #pd.DataFrame( {"t":time, "T":temperature} ).to_csv( fullPath,
                    #index=False)
            #print("Wrote reference solution to " + fileName + ".")


    plt.xlabel("Time")
    plt.ylabel("T")
    plt.plot( times, f(np.array(times)), linewidth=0.5, linestyle='--', label='Exact')
    plt.legend()
    plt.show()

def retrieveDataSets(dataFolder):
    ##get it
    dataFiles = os.listdir( dataFolder )
    ##organize it
    output = {"BE" : {},
              "BDF2" : {},
              "BDF3" : {},
              "BDF4" : {}}

    for f in dataFiles:
        fileName = dataFolder + "/" + f
        f = f.replace("_", ".")
        nsteps = re.search(r"nsteps(\d+)", f )
        if nsteps:
            nsteps = float(nsteps.group(1))
        else:
            print("Wrong file formatting for file {}".format(
                f))
            exit()
        df = pd.read_csv( fileName )
        for key in output.keys():
            if f.startswith(key):
                print(f)
                output[key][nsteps] = df
    return output

def convergenceOrder():
    time = 0.0
    fileName = "ode.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    nnodes = int(d["nels"]) + 1
    maxStepBacks = 4
    # get solution
    f = computeAnalyticalSolution()
    DTs = [0.25, 0.2, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.002, 0.001, 0.0005]
    output = {"BE" : {},
              "BDF2" : {},
              "BDF3" : {},
              "BDF4" : {}}
    for dt in DTs:
        prevSols = np.ones((nnodes, maxStepBacks))
        for stepBack in range(maxStepBacks):
            prevSols[:, stepBack] = f(time-dt*stepBack)*np.ones(nnodes)


        exactSol = f(d["Tfinal"])
        for TI in [1, 2, 3, 4]:
            d["timeIntegration"] = TI
            d["dt"] = dt
            p = mhs.Problem()
            p.initialize( d )
            p.initializeIntegrator( prevSols )
            p.time = time
            _, Ts = advanceODE_untilTfinal( p, d["Tfinal"])
            output[labelsMapping[TI]][dt] = abs( exactSol - Ts[-1] )
    colors = ["red", "blue", "magenta", "green"]
    plt.figure(dpi=200)
    for idx, (TIlabel, TIdic) in enumerate(output.items()):
        DTs = np.array( list(TIdic.keys()) )
        Errors = np.array( list(TIdic.values()) )
        #sort
        p = np.argsort( DTs )
        DTs = DTs[p]
        Errors = Errors[p]
        order = np.log( Errors[-1] / Errors[-2] ) / np.log( DTs[-1] / DTs[-2] )
        cte = Errors[-1] / np.power(DTs[-1], order)
        plt.loglog( DTs, Errors,
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
    plt.ylabel("Error")
    plt.xlabel("dt")
    plt.legend()
    plt.show()

def postProcess(simDict):
    # retrieve reference solution
    #L2Errors = dict( simDict )
    #analyticalSol = computeAnalyticalSolution()
    #for TIlabel, TI in simDict.items():
        #for CFL, p in TI.items():
            #L2Errors[TIlabel][CFL] = computeErrorAtPoint( p, analyticalSol )
    f = computeAnalyticalSolution()
    ### plot
    # plot one
    nsteps = 20
    plt.figure(dpi=200)
    for TI in simDict.keys():
        try:
            df = simDict[TI][nsteps]
        except KeyError:
            continue
        plt.scatter( df["t"], df["T"],
                label=TI)
    plt.plot(simDict["BDF2"][20]["t"], f(simDict["BDF2"][20]["t"]),
            linestyle='--',
            color='red',
            label="Analytical solution")
    plt.legend()
    plt.show()

if __name__=="__main__":
    #compute data
    dataFolder="data/debugBDF3"
    os.makedirs( dataFolder, exist_ok=True )
    computeDataSets(dataFolder, computeRefSolution=True)
    # get it
    #output = retrieveDataSets(dataFolder)
#
    #postProcess( output )
    #convergenceOrder()
