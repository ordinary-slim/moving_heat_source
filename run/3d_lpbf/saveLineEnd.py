from paraview.simple import *
import argparse

def main(dataSet, isCoupled):
    dataSetName = dataSet.split(".")[0]
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'PVD Reader'
    coupled_safepvd = PVDReader(registrationName='coupled_safe.pvd', FileName=dataSet)
    coupled_safepvd.CellArrays = ['ActiveElements', 'material', 'physicalDomain']
    coupled_safepvd.PointArrays = ['T', 'Pulse', 'ActiveNodes', 'gammaNodes', 'forcedDofs']

    UpdatePipeline(time=20.100225, proxy=coupled_safepvd)

    # create a new 'Threshold'
    thresholdBy = 'physicalDomain'
    if not(isCoupled):
        thresholdBy = 'ActiveElements'
    threshold1 = Threshold(registrationName='Threshold1', Input=coupled_safepvd)
    threshold1.Scalars = ['CELLS',thresholdBy]
    threshold1.LowerThreshold = 0.1
    threshold1.UpperThreshold = 1.0

    UpdatePipeline(time=20.100225, proxy=threshold1)

    # get animation scene
    animationScene1 = GetAnimationScene()
    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # create a new 'Plot Over Line'
    plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=threshold1)
    plotOverLine1.Point1 = [-6.0, 0.0, 0.24999]
    plotOverLine1.Point2 = [6.0, 0.0, 0.24999]

    UpdatePipeline(time=20.100225, proxy=plotOverLine1)

    # save data
    PointDataArrays = ['T', 'material']
    if isCoupled:
        PointDataArrays.append( 'gammaNodes' )
    SaveData('/home/mslimani/Documents/dd_methodology_paper/figures/3d_lpbf/plots/endTenthLayer_{}.csv'.format(dataSetName), proxy=plotOverLine1, ChooseArraysToWrite=1,
        PointDataArrays=PointDataArrays, AddTime=1)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dataSet')
    parser.add_argument('-c', '--is-coupled', action='store_true')
    args = parser.parse_args()
    main( args.dataSet, args.is_coupled )
