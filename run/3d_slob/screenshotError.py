# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11
#### import the simple module from the paraview
from paraview.simple import *
import argparse
import numpy as np

maxError = 100

def takeScreenshot( fileName ):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'PVD Reader'
    analyticalpvd = PVDReader(registrationName='analytical.pvd', FileName='/home/mslimani/acuario/moving_heat_source/run/3d_slob/analyticalSeparateMesh/analytical.pvd')
    analyticalpvd.PointArrays = ['T']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    analyticalpvdDisplay = Show(analyticalpvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    analyticalpvdDisplay.Representation = 'Surface'
    analyticalpvdDisplay.ColorArrayName = [None, '']
    analyticalpvdDisplay.SelectTCoordArray = 'None'
    analyticalpvdDisplay.SelectNormalArray = 'None'
    analyticalpvdDisplay.SelectTangentArray = 'None'
    analyticalpvdDisplay.OSPRayScaleArray = 'T'
    analyticalpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    analyticalpvdDisplay.SelectOrientationVectors = 'None'
    analyticalpvdDisplay.ScaleFactor = 0.05
    analyticalpvdDisplay.SelectScaleArray = 'None'
    analyticalpvdDisplay.GlyphType = 'Arrow'
    analyticalpvdDisplay.GlyphTableIndexArray = 'None'
    analyticalpvdDisplay.GaussianRadius = 0.0025
    analyticalpvdDisplay.SetScaleArray = ['POINTS', 'T']
    analyticalpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    analyticalpvdDisplay.OpacityArray = ['POINTS', 'T']
    analyticalpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    analyticalpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    analyticalpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    analyticalpvdDisplay.ScalarOpacityUnitDistance = 0.013226968584849428
    analyticalpvdDisplay.OpacityArrayName = ['POINTS', 'T']
    analyticalpvdDisplay.SelectInputVectors = [None, '']
    analyticalpvdDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    analyticalpvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    analyticalpvdDisplay.ScaleTransferFunction.Points = [25.000000019037056, 0.0, 0.5, 0.0, 2382.7323929821678, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    analyticalpvdDisplay.OpacityTransferFunction.Points = [25.000000019037056, 0.0, 0.5, 0.0, 2382.7323929821678, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=analyticalpvd)
    calculator1.Function = ''

    # Properties modified on calculator1
    calculator1.ResultArrayName = 'Tanal'
    calculator1.Function = 'T'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

    # get 2D transfer function for 'Tanal'
    tanalTF2D = GetTransferFunction2D('Tanal')

    # get color transfer function/color map for 'Tanal'
    tanalLUT = GetColorTransferFunction('Tanal')
    tanalLUT.TransferFunction2D = tanalTF2D
    tanalLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 405.06450680696855, 0.0, 0.0, 1.0, 810.1290136139371, 0.0, 1.0, 1.0, 1191.3661964910839, 0.0, 1.0, 0.0, 1596.4307032980525, 1.0, 1.0, 0.0, 2001.4952101050205, 1.0, 0.0, 0.0, 2382.7323929821678, 0.878431372549, 0.0, 1.0]
    tanalLUT.ColorSpace = 'RGB'
    tanalLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'Tanal'
    tanalPWF = GetOpacityTransferFunction('Tanal')
    tanalPWF.Points = [0.0, 0.0, 0.5, 0.0, 2382.7323929821678, 1.0, 0.5, 0.0]
    tanalPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Tanal']
    calculator1Display.LookupTable = tanalLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'Tanal'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'None'
    calculator1Display.ScaleFactor = 0.05
    calculator1Display.SelectScaleArray = 'Tanal'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Tanal'
    calculator1Display.GaussianRadius = 0.0025
    calculator1Display.SetScaleArray = ['POINTS', 'Tanal']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Tanal']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'
    calculator1Display.ScalarOpacityFunction = tanalPWF
    calculator1Display.ScalarOpacityUnitDistance = 0.013226968584849428
    calculator1Display.OpacityArrayName = ['POINTS', 'Tanal']
    calculator1Display.SelectInputVectors = [None, '']
    calculator1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    calculator1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [25.000000019037056, 0.0, 0.5, 0.0, 2382.7323929821678, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [25.000000019037056, 0.0, 0.5, 0.0, 2382.7323929821678, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(analyticalpvd, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'PVD Reader'
    dataSet = PVDReader(registrationName="dataSet", FileName=fileName)
    dataSet.CellArrays = ['ActiveElements', 'physicalDomain']
    dataSet.PointArrays = ['T', 'Pulse', 'ActiveNodes', 'gammaNodes', 'forcedDofs']

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # show data in view
    dataSetDisplay = Show(dataSet, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    dataSetDisplay.Representation = 'Surface'
    dataSetDisplay.ColorArrayName = [None, '']
    dataSetDisplay.SelectTCoordArray = 'None'
    dataSetDisplay.SelectNormalArray = 'None'
    dataSetDisplay.SelectTangentArray = 'None'
    dataSetDisplay.OSPRayScaleArray = 'ActiveNodes'
    dataSetDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    dataSetDisplay.SelectOrientationVectors = 'None'
    dataSetDisplay.ScaleFactor = 0.36
    dataSetDisplay.SelectScaleArray = 'ActiveNodes'
    dataSetDisplay.GlyphType = 'Arrow'
    dataSetDisplay.GlyphTableIndexArray = 'ActiveNodes'
    dataSetDisplay.GaussianRadius = 0.018
    dataSetDisplay.SetScaleArray = ['POINTS', 'ActiveNodes']
    dataSetDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    dataSetDisplay.OpacityArray = ['POINTS', 'ActiveNodes']
    dataSetDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    dataSetDisplay.DataAxesGrid = 'GridAxesRepresentation'
    dataSetDisplay.PolarAxes = 'PolarAxesRepresentation'
    dataSetDisplay.ScalarOpacityUnitDistance = 0.11913968780301959
    dataSetDisplay.OpacityArrayName = ['POINTS', 'ActiveNodes']
    dataSetDisplay.SelectInputVectors = [None, '']
    dataSetDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    dataSetDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Calculator'
    calculator2 = Calculator(registrationName='Calculator2', Input=dataSet)
    calculator2.Function = ''

    # Properties modified on renderView1
    renderView1.CenterOfRotation = [-0.14999999478459358, 0.0494999997317791, 0.0]
    renderView1.CameraPosition = [-0.11521622491233216, 0.041634347699672225, 0.21346981293248826]
    renderView1.CameraFocalPoint = [-0.11521622491233216, 0.041634347699672225, 0.0]
    renderView1.CameraParallelScale = 0.2538725016795355

    # Properties modified on calculator2
    calculator2.ResultArrayName = 'Th'
    calculator2.Function = 'T'

    # show data in view
    calculator2Display = Show(calculator2, renderView1, 'UnstructuredGridRepresentation')

    # get 2D transfer function for 'Th'
    thTF2D = GetTransferFunction2D('Th')

    # get color transfer function/color map for 'Th'
    thLUT = GetColorTransferFunction('Th')
    thLUT.TransferFunction2D = thTF2D
    thLUT.RGBPoints = [24.345052080314115, 1.0, 1.0, 1.0, 424.1884556287414, 0.0, 0.0, 1.0, 824.0318591771686, 0.0, 1.0, 1.0, 1200.3550625168648, 0.0, 1.0, 0.0, 1600.198466065292, 1.0, 1.0, 0.0, 2000.041869613719, 1.0, 0.0, 0.0, 2376.365072953415, 0.878431372549, 0.0, 1.0]
    thLUT.ColorSpace = 'RGB'
    thLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'Th'
    thPWF = GetOpacityTransferFunction('Th')
    thPWF.Points = [24.345052080314115, 0.0, 0.5, 0.0, 2376.365072953415, 1.0, 0.5, 0.0]
    thPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    calculator2Display.Representation = 'Surface'
    calculator2Display.ColorArrayName = ['POINTS', 'Th']
    calculator2Display.LookupTable = thLUT
    calculator2Display.SelectTCoordArray = 'None'
    calculator2Display.SelectNormalArray = 'None'
    calculator2Display.SelectTangentArray = 'None'
    calculator2Display.OSPRayScaleArray = 'Th'
    calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator2Display.SelectOrientationVectors = 'None'
    calculator2Display.ScaleFactor = 0.36
    calculator2Display.SelectScaleArray = 'Th'
    calculator2Display.GlyphType = 'Arrow'
    calculator2Display.GlyphTableIndexArray = 'Th'
    calculator2Display.GaussianRadius = 0.018
    calculator2Display.SetScaleArray = ['POINTS', 'Th']
    calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator2Display.OpacityArray = ['POINTS', 'Th']
    calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator2Display.PolarAxes = 'PolarAxesRepresentation'
    calculator2Display.ScalarOpacityFunction = thPWF
    calculator2Display.ScalarOpacityUnitDistance = 0.11913968780301959
    calculator2Display.OpacityArrayName = ['POINTS', 'Th']
    calculator2Display.SelectInputVectors = [None, '']
    calculator2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    calculator2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator2Display.ScaleTransferFunction.Points = [24.44252235453569, 0.0, 0.5, 0.0, 1784.9303468883584, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator2Display.OpacityTransferFunction.Points = [24.44252235453569, 0.0, 0.5, 0.0, 1784.9303468883584, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(dataSet, renderView1)

    # show color bar/color legend
    calculator2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(calculator1)

    # create a new 'Resample With Dataset'
    resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=calculator1,
        DestinationMesh=calculator2)
    resampleWithDataset1.CellLocator = 'Static Cell Locator'

    # Properties modified on resampleWithDataset1
    resampleWithDataset1.PassPointArrays = 1

    # show data in view
    resampleWithDataset1Display = Show(resampleWithDataset1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    resampleWithDataset1Display.Representation = 'Surface'
    resampleWithDataset1Display.ColorArrayName = ['POINTS', 'Tanal']
    resampleWithDataset1Display.LookupTable = tanalLUT
    resampleWithDataset1Display.SelectTCoordArray = 'None'
    resampleWithDataset1Display.SelectNormalArray = 'None'
    resampleWithDataset1Display.SelectTangentArray = 'None'
    resampleWithDataset1Display.OSPRayScaleArray = 'Tanal'
    resampleWithDataset1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    resampleWithDataset1Display.SelectOrientationVectors = 'None'
    resampleWithDataset1Display.ScaleFactor = 0.36
    resampleWithDataset1Display.SelectScaleArray = 'Tanal'
    resampleWithDataset1Display.GlyphType = 'Arrow'
    resampleWithDataset1Display.GlyphTableIndexArray = 'Tanal'
    resampleWithDataset1Display.GaussianRadius = 0.018
    resampleWithDataset1Display.SetScaleArray = ['POINTS', 'Tanal']
    resampleWithDataset1Display.ScaleTransferFunction = 'PiecewiseFunction'
    resampleWithDataset1Display.OpacityArray = ['POINTS', 'Tanal']
    resampleWithDataset1Display.OpacityTransferFunction = 'PiecewiseFunction'
    resampleWithDataset1Display.DataAxesGrid = 'GridAxesRepresentation'
    resampleWithDataset1Display.PolarAxes = 'PolarAxesRepresentation'
    resampleWithDataset1Display.ScalarOpacityFunction = tanalPWF
    resampleWithDataset1Display.ScalarOpacityUnitDistance = 0.11913968780301959
    resampleWithDataset1Display.OpacityArrayName = ['POINTS', 'Tanal']
    resampleWithDataset1Display.SelectInputVectors = [None, '']
    resampleWithDataset1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    resampleWithDataset1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleWithDataset1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2373.3980777660922, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleWithDataset1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2373.3980777660922, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(calculator2, renderView1)

    # hide data in view
    Hide(calculator1, renderView1)

    # show color bar/color legend
    resampleWithDataset1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Plane'
    plane1 = Plane(registrationName='Plane1')

    # Properties modified on plane1
    plane1.Origin = [-0.399, 0.0, 0.0]
    plane1.Point1 = [0.099, 0.0, 0.0]
    plane1.Point2 = [-0.399, 0.0999, 0.0]
    plane1.XResolution = 251
    plane1.YResolution = 101

    # show data in view
    plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    plane1Display.Representation = 'Surface'
    plane1Display.ColorArrayName = [None, '']
    plane1Display.SelectTCoordArray = 'TextureCoordinates'
    plane1Display.SelectNormalArray = 'Normals'
    plane1Display.SelectTangentArray = 'None'
    plane1Display.OSPRayScaleArray = 'Normals'
    plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plane1Display.SelectOrientationVectors = 'None'
    plane1Display.ScaleFactor = 0.04979999884963036
    plane1Display.SelectScaleArray = 'None'
    plane1Display.GlyphType = 'Arrow'
    plane1Display.GlyphTableIndexArray = 'None'
    plane1Display.GaussianRadius = 0.002489999942481518
    plane1Display.SetScaleArray = ['POINTS', 'Normals']
    plane1Display.ScaleTransferFunction = 'PiecewiseFunction'
    plane1Display.OpacityArray = ['POINTS', 'Normals']
    plane1Display.OpacityTransferFunction = 'PiecewiseFunction'
    plane1Display.DataAxesGrid = 'GridAxesRepresentation'
    plane1Display.PolarAxes = 'PolarAxesRepresentation'
    plane1Display.SelectInputVectors = ['POINTS', 'Normals']
    plane1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    plane1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    plane1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    plane1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(resampleWithDataset1)

    # create a new 'Resample With Dataset'
    resampleWithDataset2 = ResampleWithDataset(registrationName='ResampleWithDataset2', SourceDataArrays=resampleWithDataset1,
        DestinationMesh=plane1)
    resampleWithDataset2.CellLocator = 'Static Cell Locator'

    # Properties modified on resampleWithDataset2
    resampleWithDataset2.PassPointArrays = 1

    # show data in view
    resampleWithDataset2Display = Show(resampleWithDataset2, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    resampleWithDataset2Display.Representation = 'Surface'
    resampleWithDataset2Display.ColorArrayName = ['POINTS', 'Tanal']
    resampleWithDataset2Display.LookupTable = tanalLUT
    resampleWithDataset2Display.SelectTCoordArray = 'TextureCoordinates'
    resampleWithDataset2Display.SelectNormalArray = 'Normals'
    resampleWithDataset2Display.SelectTangentArray = 'None'
    resampleWithDataset2Display.OSPRayScaleArray = 'Tanal'
    resampleWithDataset2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    resampleWithDataset2Display.SelectOrientationVectors = 'None'
    resampleWithDataset2Display.ScaleFactor = 0.04979999884963036
    resampleWithDataset2Display.SelectScaleArray = 'Tanal'
    resampleWithDataset2Display.GlyphType = 'Arrow'
    resampleWithDataset2Display.GlyphTableIndexArray = 'Tanal'
    resampleWithDataset2Display.GaussianRadius = 0.002489999942481518
    resampleWithDataset2Display.SetScaleArray = ['POINTS', 'Tanal']
    resampleWithDataset2Display.ScaleTransferFunction = 'PiecewiseFunction'
    resampleWithDataset2Display.OpacityArray = ['POINTS', 'Tanal']
    resampleWithDataset2Display.OpacityTransferFunction = 'PiecewiseFunction'
    resampleWithDataset2Display.DataAxesGrid = 'GridAxesRepresentation'
    resampleWithDataset2Display.PolarAxes = 'PolarAxesRepresentation'
    resampleWithDataset2Display.SelectInputVectors = ['POINTS', 'Normals']
    resampleWithDataset2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    resampleWithDataset2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleWithDataset2Display.ScaleTransferFunction.Points = [1.3343963787527604, 0.0, 0.5, 0.0, 2371.60974577615, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleWithDataset2Display.OpacityTransferFunction.Points = [1.3343963787527604, 0.0, 0.5, 0.0, 2371.60974577615, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(plane1, renderView1)

    # hide data in view
    Hide(resampleWithDataset1, renderView1)

    # show color bar/color legend
    resampleWithDataset2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on renderView1
    renderView1.CameraPosition = [-0.11521622491233216, 0.041634347699672225, 0.21346981293248826]

    # Properties modified on renderView1
    renderView1.OrientationAxesVisibility = 0

    # create a new 'Calculator'
    calculator3 = Calculator(registrationName='Calculator3', Input=resampleWithDataset2)
    calculator3.Function = ''

    # Properties modified on calculator3
    calculator3.ResultArrayName = 'Error'
    calculator3.Function = '(Th - Tanal)^2'

    # show data in view
    calculator3Display = Show(calculator3, renderView1, 'GeometryRepresentation')

    # get 2D transfer function for 'Error'
    errorTF2D = GetTransferFunction2D('Error')

    # get color transfer function/color map for 'Error'
    errorLUT = GetColorTransferFunction('Error')
    errorLUT.TransferFunction2D = errorTF2D
    errorLUT.RGBPoints = [1.931156797141398e-05, 1.0, 1.0, 1.0, 185.64793887269494, 0.0, 0.0, 1.0, 371.29585843382193, 0.0, 1.0, 1.0, 546.0233121384119, 0.0, 1.0, 0.0, 731.6712316995389, 1.0, 1.0, 0.0, 917.3191512606656, 1.0, 0.0, 0.0, 1092.046604965256, 0.878431372549, 0.0, 1.0]
    errorLUT.ColorSpace = 'RGB'
    errorLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    calculator3Display.Representation = 'Surface'
    calculator3Display.ColorArrayName = ['POINTS', 'Error']
    calculator3Display.LookupTable = errorLUT
    calculator3Display.SelectTCoordArray = 'TextureCoordinates'
    calculator3Display.SelectNormalArray = 'Normals'
    calculator3Display.SelectTangentArray = 'None'
    calculator3Display.OSPRayScaleArray = 'Error'
    calculator3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator3Display.SelectOrientationVectors = 'None'
    calculator3Display.ScaleFactor = 0.04979999884963036
    calculator3Display.SelectScaleArray = 'Error'
    calculator3Display.GlyphType = 'Arrow'
    calculator3Display.GlyphTableIndexArray = 'Error'
    calculator3Display.GaussianRadius = 0.002489999942481518
    calculator3Display.SetScaleArray = ['POINTS', 'Error']
    calculator3Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator3Display.OpacityArray = ['POINTS', 'Error']
    calculator3Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator3Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator3Display.PolarAxes = 'PolarAxesRepresentation'
    calculator3Display.SelectInputVectors = ['POINTS', 'Normals']
    calculator3Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    calculator3Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator3Display.ScaleTransferFunction.Points = [1.931156797141398e-05, 0.0, 0.5, 0.0, 1092.046604965256, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator3Display.OpacityTransferFunction.Points = [1.931156797141398e-05, 0.0, 0.5, 0.0, 1092.046604965256, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(resampleWithDataset2, renderView1)

    # show color bar/color legend
    calculator3Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get opacity transfer function/opacity map for 'Error'
    errorPWF = GetOpacityTransferFunction('Error')
    errorPWF.Points = [1.931156797141398e-05, 0.0, 0.5, 0.0, 1092.046604965256, 1.0, 0.5, 0.0]
    errorPWF.ScalarRangeInitialized = 1

    # Properties modified on animationScene1
    animationScene1.AnimationTime = 0.0022

    # Rescale transfer function
    errorLUT.RescaleTransferFunction(0.0, 100.0)

    # Rescale transfer function
    errorPWF.RescaleTransferFunction(0.0, 100.0)

    # Rescale 2D transfer function
    errorTF2D.RescaleTransferFunction(0.0, 100.0, 0.0, 1.0)

    # get color legend/bar for errorLUT in view renderView1
    errorLUTColorBar = GetScalarBar(errorLUT, renderView1)
    errorLUTColorBar.Orientation = 'Vertical'
    errorLUTColorBar.WindowLocation = 'Any Location'
    errorLUTColorBar.Title = 'Error'
    errorLUTColorBar.ComponentTitle = ''
    errorLUTColorBar.TitleColor = [0.0, 0.0, 0.16]
    errorLUTColorBar.TitleFontFamily = 'Times'
    errorLUTColorBar.TitleFontSize = 30
    errorLUTColorBar.LabelColor = [0.0, 0.0, 0.16000610360875867]
    errorLUTColorBar.LabelFontFamily = 'Times'
    errorLUTColorBar.LabelFontSize = 30
    errorLUTColorBar.ScalarBarLength = 0.33000000000000096
    errorLUTColorBar.RangeLabelFormat = '%#.2f'

    # change scalar bar placement
    errorLUTColorBar.Position = [0.8868474148802018, 0.445]
    errorLUTColorBar.ScalarBarLength = 0.330000000000001

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=calculator3)
    contour1.ContourBy = ['POINTS', 'Error']
    contour1.Isosurfaces = [190.85209656718354]
    contour1.PointMergeMethod = 'Uniform Binning'

    # Properties modified on contour1
    contour1.ContourBy = ['POINTS', 'Th']
    contour1.Isosurfaces = [700.0, 1100.0, 1500.0, 1900.0, 2300.0, 2100.0]

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = ['POINTS', 'Th']
    contour1Display.LookupTable = thLUT
    contour1Display.SelectTCoordArray = 'TextureCoordinates'
    contour1Display.SelectNormalArray = 'Normals'
    contour1Display.SelectTangentArray = 'None'
    contour1Display.OSPRayScaleArray = 'Th'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 0.04318874478340149
    contour1Display.SelectScaleArray = 'Th'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'Th'
    contour1Display.GaussianRadius = 0.0021594372391700745
    contour1Display.SetScaleArray = ['POINTS', 'Th']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'Th']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'
    contour1Display.SelectInputVectors = ['POINTS', 'Normals']
    contour1Display.WriteLog = ''
    contour1Display.Opacity = 0.4
    contour1Display.LineWidth = 3

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    contour1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(calculator3, renderView1)

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # turn off scalar coloring
    ColorBy(contour1Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(thLUT, renderView1)

    # change solid color
    contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

    # set active source
    SetActiveSource(calculator3)

    # show data in view
    calculator3Display = Show(calculator3, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    calculator3Display.SetScalarBarVisibility(renderView1, True)

    # Properties modified on errorLUTColorBar
    errorLUTColorBar.Title = 'Error'

    # set active source
    SetActiveSource(plane1)

    # set active source
    SetActiveSource(calculator3)

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    errorLUT.ApplyPreset('Rainbow Uniform', True)

    # Rescale transfer function
    errorLUT.RescaleTransferFunction(0.0, maxError)

    # Rescale transfer function
    errorPWF.RescaleTransferFunction(0.0, maxError)

    # Rescale 2D transfer function
    errorTF2D.RescaleTransferFunction(0.0, maxError, 0.0, 1.0)

    # change scalar bar placement
    errorLUTColorBar.Position = [0.8868474148802018, 0.2325]
    errorLUTColorBar.ScalarBarLength = 0.542500000000001

    # Properties modified on errorLUTColorBar
    errorLUTColorBar.RangeLabelFormat = '%#.0f'

    # Show mesh
    # create a new 'Slice'
    sliceMesh = Slice(registrationName='sliceMesh', Input=dataSet)
    sliceMesh.SliceType = 'Plane'
    sliceMesh.HyperTreeGridSlicer = 'Plane'
    sliceMesh.SliceOffsetValues = [0.0]
    sliceMesh.SliceType.Origin = [0.0, 0.0, 0.0]
    sliceMesh.SliceType.Normal = [0.0, 0.0, 1.0]
    sliceMesh.Triangulatetheslice = 0
    Hide(sliceMesh, renderView1)
    clipMesh = Clip(registrationName='clipMesh', Input=sliceMesh)
    clipMesh.ClipType = 'Plane'
    clipMesh.HyperTreeGridClipper = 'Plane'
    clipMesh.Scalars = ['POINTS', 'ActiveNodes']
    clipMesh.Value = 0.5
    clipMesh.Invert = 0
    clipMesh.ClipType.Origin = [0.0999999, 0.0, 0.0]
    clipMesh.ClipType.Normal = [-1.0, 0.0, 0.0]
    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clipMesh.HyperTreeGridClipper.Origin = [-0.9999999999999993, 0.5, 0.0]
    Hide(clipMesh, renderView1)

    extractEdges1 = ExtractEdges(registrationName='ExtractEdges1', Input=clipMesh)
    extractEdges1Display = Show(extractEdges1, renderView1, 'GeometryRepresentation')
    extractEdges1Display.Representation = 'Surface'
    extractEdges1Display.ColorArrayName = [None, '']
    extractEdges1Display.Opacity = 0.2
    extractEdges1Display.DiffuseColor = [0.0, 0.0, 0.0]

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1586, 320)

    # current camera placement for renderView1
    renderView1.CameraPosition = [-0.11521622491233216, 0.041634347699672225, 0.21346981293248826]
    renderView1.CameraFocalPoint = [-0.11521622491233216, 0.041634347699672225, 0.0]
    renderView1.CameraParallelScale = 0.2538725016795355

    # save screenshot
    SaveScreenshot('./figures/err_{}.png'.format(
        fileName.split(".")[0] ),
                   renderView1,
                   ImageResolution=[1586, 320],
        TransparentBackground=1)

    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================

    #--------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(1586, 320)

    #-----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [-0.11521622491233216, 0.041634347699672225, 0.21346981293248826]
    renderView1.CameraFocalPoint = [-0.11521622491233216, 0.041634347699672225, 0.0]
    renderView1.CameraParallelScale = 0.2538725016795355

    #--------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).

    # Compute L2 error
    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator3)
    l2Err    = np.sqrt( integrateVariables1.PointData["Error"].GetRange()[0] )
    normAnal = np.sqrt( integrateVariables1.PointData["Tanal"].GetRange()[0] )
    with open("l2Errors.txt", "a") as l2File:
        l2File.write("{}: {}\n".format( fileName, l2Err ) )
    with open("relL2Errors.txt", "a") as l2File:
        l2File.write("{}: {}\n".format( fileName, 100*l2Err/normAnal ) )

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dataSet')
    args = parser.parse_args()
    takeScreenshot( args.dataSet )
