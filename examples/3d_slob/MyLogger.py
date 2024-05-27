import time
import argparse
import pickle

class IterationData:
    def __init__(self):
        self.isCoupled   = False
        self.isPrinting  = False
        self.ndofs  = -1
        self.startTime = 0.0
        self.endTime = +1e9
        self.duration = +1e9

class MyLogger:
    def __init__(self):
        self.iterations = []
    def iterate(self, driver):
        iterData = IterationData()
        iterData.startTime = time.time()
        driver.iterate()
        iterData.endTime = time.time()
        iterData.duration = iterData.endTime - iterData.startTime
        iterData.isCoupled = driver.isCoupled
        iterData.isPrinting = driver.getIsPrinting()
        iterData.ndofs = driver.getNdofs()
        self.iterations.append( iterData )

def processLog( logFile ):
    with open( logFile, 'rb' ) as lf:
        log = pickle.load( lf )
    numIters = len( log.iterations )
    numPrintIters = 0
    totTimePrinting = 0
    avDurationPrintingIter = 0
    numCoupledIters = 0
    avDurationCoupledIter = 0
    avDurationUncoupledIter = 0
    avNDofs = 0
    avNDofsUncoupled = 0
    avNDofsCoupled = 0
    for iteration in log.iterations:
        if iteration.isPrinting:
            numPrintIters += 1
            avDurationPrintingIter += iteration.duration
            avNDofs += iteration.ndofs
            if iteration.isCoupled:
                numCoupledIters += 1
                avDurationCoupledIter += iteration.duration
                avNDofsCoupled += iteration.ndofs
            else:
                avDurationUncoupledIter += iteration.duration
                avNDofsUncoupled += iteration.ndofs
    totTimePrinting = avDurationPrintingIter
    avDurationPrintingIter /= numPrintIters
    avNDofs /= numPrintIters
    try:
        avDurationCoupledIter /= numCoupledIters
        avNDofsCoupled /= numCoupledIters
    except ZeroDivisionError:
        avDurationCoupledIter = None
        avNDofsCoupled = None
    try:
        avDurationUncoupledIter /= (numPrintIters - numCoupledIters)
        avNDofsUncoupled /= (numPrintIters - numCoupledIters)
    except ZeroDivisionError:
        avDurationUncoupledIter = None
        avNDofsUncoupled = None
    print( "{} printings iters out of {}".format( numPrintIters, numIters ) )
    print( "{} coupled   iters out of {}".format( numCoupledIters, numIters ) )
    print( "Total time printing iterations = {}".format( totTimePrinting ) )
    print( "Average ndofs printing iter           = {}".format( avNDofs ) )
    print( "Average ndofs coupled  iter           = {}".format( avNDofsCoupled ) )
    print( "Average ndofs uncoupled printing iter = {}".format( avNDofsUncoupled ) )
    print( "Average duration printing iter           = {}".format( avDurationPrintingIter ) )
    print( "Average duration coupled  iter           = {}".format( avDurationCoupledIter ) )
    print( "Average duration uncoupled printing iter = {}".format( avDurationUncoupledIter ) )


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("logfile")
    args = parser.parse_args()
    processLog( args.logfile )
