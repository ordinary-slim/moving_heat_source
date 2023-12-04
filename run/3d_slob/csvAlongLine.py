# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
coupled_elsPerRad8_10Rinterfacepvd = PVDReader(registrationName='coupled_elsPerRad8_10Rinterface.pvd', FileName='/home/mslimani/acuario/moving_heat_source/run/3d_slob/coupled_elsPerRad8_10Rinterface.pvd')
coupled_elsPerRad8_10Rinterfacepvd.CellArrays = ['ActiveElements', 'physicalDomain']
coupled_elsPerRad8_10Rinterfacepvd.PointArrays = ['T', 'Pulse', 'ActiveNodes', 'gammaNodes', 'forcedDofs']

UpdatePipeline(time=0.0022, proxy=coupled_elsPerRad8_10Rinterfacepvd)

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=coupled_elsPerRad8_10Rinterfacepvd)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-0.9999999999999993, 0.5, -0.59375]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [-0.9999999999999993, 0.5, -0.59375]

# Properties modified on slice1
slice1.Triangulatetheslice = 0

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

UpdatePipeline(time=0.0022, proxy=slice1)

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=slice1)
plotOverLine1.Point1 = [-2.7999999999999994, 0.0, 0.0]
plotOverLine1.Point2 = [0.8000000000000007, 1.0, 0.0]
plotOverLine1.Resolution = 200

# Properties modified on plotOverLine1
plotOverLine1.Point1 = [-1.2, 0.0, 0.0]
plotOverLine1.Point2 = [0.5, 0.0, 0.0]

UpdatePipeline(time=0.0022, proxy=plotOverLine1)

# save data
SaveData('/home/mslimani/Documents/dd_methodology_paper/figures/3d_weldingLpbf/plots/farInterface.csv', proxy=plotOverLine1, ChooseArraysToWrite=1,
    PointDataArrays=['T', 'gammaNodes'])
