import time
import argparse
import pickle

class IterationData:
    def __init__(self):
        self.isCoupled   = False
        self.isPrinting  = False
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
    for iteration in log.iterations:
        if iteration.isPrinting:
            numPrintIters += 1
            avDurationPrintingIter += iteration.duration
        if iteration.isCoupled:
            numCoupledIters += 1
            avDurationCoupledIter += iteration.duration
        elif iteration.isPrinting:
            avDurationUncoupledIter += iteration.duration
    totTimePrinting = avDurationPrintingIter
    avDurationPrintingIter /= numPrintIters
    try:
        avDurationCoupledIter /= numCoupledIters
    except ZeroDivisionError:
        avDurationCoupledIter = None
    avDurationUncoupledIter /= (numPrintIters - numCoupledIters)
    print( "{} printings iters out of {}".format( numPrintIters, numIters ) )
    print( "{} coupled   iters out of {}".format( numCoupledIters, numIters ) )
    print( "Total time printing iterations = {}".format( totTimePrinting ) )
    print( "Average duration printing iter           = {}".format( avDurationPrintingIter ) )
    print( "Average duration coupled  iter           = {}".format( avDurationCoupledIter ) )
    print( "Average duration uncoupled printing iter = {}".format( avDurationUncoupledIter ) )



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("logfile")
    args = parser.parse_args()
    processLog( args.logfile )
