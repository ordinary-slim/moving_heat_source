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

def computeL2Error( p1, p2 ):
    # safety check
    if (p1.time != p2.time):
        print("Something wrong!")
        exit()
    return np.sqrt( np.trapz( np.square( p1.solution - p2.solution ), p1.mesh.pos ) )

def advanceUntilTfinal( p, Tfinal ):
    tol = 1e-4
    while( Tfinal - p.time > tol ):
        print( "Current time:", p.time )
        p.iterate()

def getCFL( d ):
    L = d["Right"] - d["Left"]
    nels = d["nels"]
    h = L/float(nels)
    dt = d["dt"]
    speed = d["speedX"]
    return (speed * dt) / h

def set_dt( d, CFL ):
    L = d["Right"] - d["Left"]
    nels = d["nels"]
    h = L/float(nels)
    speed = d["speedX"]
    d["dt"] = (CFL * h) / float(speed)
    print( "CFL = {}, dt = {}".format( CFL, d["dt"] ) )
    return d

def writeReferenceSolution(referenceCFL=1e-2):
    CFL = referenceCFL
    fileName = "input.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    Tfinal = d["Tfinal"]
    fileName = "referenceSolution_CFL{}_T{}.csv".format(
            formatFloat( CFL ),
            formatFloat( Tfinal ))
    postFolder="data"
    # check for file existence
    if os.path.isfile( "./" + postFolder + "/" + fileName ):
        return

    d = set_dt( d, CFL )
    p = mhs.Problem()
    p.initialize( d )
    advanceUntilTfinal( p, Tfinal )
    pf = writePost( p, fileName, postFolder )
    print("Wrote reference solution to " + pf + ".")
    return pf


def computeDataSets(computeRefSolution=False):
    fileName = "input.txt"
    postFolder = "data"
    d = formatInputFile( fileName )
    d = parseInput( d )
    Tfinal = d["Tfinal"]
    CFLs = [100, 80, 60, 50, 40.0, 25, 20.0, 10.0, 5.0, 4.0, 2.5, 2, 1.0, 0.5, 0.1, 0.01]
    #CFLs = [40.0, 20, 10, 5, 1]
    p = mhs.Problem()

    if computeRefSolution:
        refSoluton = writeReferenceSolution(referenceCFL=0.001)

    for timeIntegration in [1, 2, 3, 4]:
        d["timeIntegration"] = timeIntegration
        for CFL in CFLs:
            d = set_dt( d, CFL )
            fileName = "{}_CFL{}_T{}.csv".format(
                    labelsMapping[timeIntegration],
                    formatFloat(CFL),
                    formatFloat( Tfinal ),
                    )
            if os.path.isfile( "./" + postFolder + "/" + fileName ):
                print("Filename {} already exists".format(
                    fileName))
                continue
            p = mhs.Problem()
            p.initialize( d )
            advanceUntilTfinal( p, Tfinal )
            pf = writePost( p, fileName, postFolder )
            print("Wrote reference solution to " + pf + ".")

def retrieveDataSets():
    ##get it
    dataFolder = "data"
    dataFiles = os.listdir( "data" )
    ##organize it
    output = {"reference" : {},
              "BE" : {},
              "BDF2" : {},
              "BDF3" : {},
              "BDF4" : {}}

    for f in dataFiles:
        fileName = dataFolder + "/" + f
        f = f.replace("_", ".")
        CFL = re.search(r"CFL(\d+(\.\d+)?)", f )
        if CFL:
            CFL = float(CFL.group(1))
        else:
            print("Wrong file formatting for file {}".format(
                f))
            exit()
        p = loadProblem( fileName )
        p.CFL = CFL
        for key in output.keys():
            if f.startswith(key):
                print(f)
                output[key][CFL] = p
    return output

if __name__=="__main__":
    #compute data
    computeDataSets()
    # get it
    output = retrieveDataSets()
    # retrieve reference solution
    referenceProblem = output.pop("reference", None)
    if referenceProblem:
        referenceProblem = referenceProblem[list(referenceProblem.keys())[0]]
    #process it
    L2Errors = dict( output )
    for TIlabel, TI in output.items():
        for CFL, p in TI.items():
            L2Errors[TIlabel][CFL] = computeL2Error( p, referenceProblem )
#
    ### plot
    plt.figure(dpi=200)
    for TIlabel, TI in L2Errors.items():
        CFLs = np.array( list(TI.keys()) )
        L2Errors = np.array( list(TI.values()) )
        #sort
        p = np.argsort( CFLs )
        CFLs = CFLs[p]
        L2Errors = L2Errors[p]


        plt.loglog( CFLs, L2Errors,
                linestyle='None',
                marker="o",
                label=TIlabel)
        plt.grid()
        plt.ylabel("L2")
        plt.xlabel("CFL")
        plt.legend()
    powers = [4, 3, 2, 1]
    for power in powers:
        alpha = L2Errors[-1] / pow( CFLs[-1], power )
        plt.loglog( CFLs, alpha*np.power(CFLs, power),
                linestyle='--',
                label="Order {}".format( power ))
    plt.show()
