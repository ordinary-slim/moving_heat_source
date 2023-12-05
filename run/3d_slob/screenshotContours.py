import argparse
from paraview.simple import *

contourValues = [700.0, 1150.0, 1400.0, 1600.0, 1900.0, 2300.0, 2100.0]

def takeScreenshot( fileName, closeInterface=False, otherMeshFile="" ):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    dataSet = None
    dataSet = PVDReader(registrationName='dataSet', FileName=(fileName) )

    dataSet.CellArrays = ['ActiveElements', 'physicalDomain']
    dataSet.PointArrays = ['T', 'Pulse', 'ActiveNodes', 'gammaNodes', 'forcedDofs']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # hide data in view
    Hide(dataSet, renderView1)

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on animationScene1
    animationScene1.AnimationTime = 0.0022

    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=dataSet)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    slice1.SliceType.Origin = [0.0, 0.0, 0.0]
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]
    slice1.Triangulatetheslice = 0

    # create a new 'Clip'
    clip1 = Clip(registrationName='Clip1', Input=slice1)

    if not(otherMeshFile):
        extractEdges1 = ExtractEdges(registrationName='ExtractEdges1', Input=clip1)
        extractEdges1Display = Show(extractEdges1, renderView1, 'GeometryRepresentation')
        extractEdges1Display.Representation = 'Surface'
        extractEdges1Display.ColorArrayName = [None, '']
        extractEdges1Display.Opacity = 0.2
        extractEdges1Display.DiffuseColor = [0.0, 0.0, 0.0]
    else:
        dataSetMesh = PVDReader(registrationName='dataSetMesh', FileName=(otherMeshFile) )
        Hide(dataSetMesh, renderView1)
        # create a new 'Slice'
        sliceMesh = Slice(registrationName='sliceMesh', Input=dataSetMesh)
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


    clip1.ClipType = 'Plane'
    clip1.HyperTreeGridClipper = 'Plane'
    clip1.Scalars = ['POINTS', 'ActiveNodes']
    clip1.Value = 0.5
    clip1.Invert = 0
    clip1.ClipType.Origin = [0.0999999, 0.0, 0.0]
    clip1.ClipType.Normal = [-1.0, 0.0, 0.0]
    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip1.HyperTreeGridClipper.Origin = [-0.9999999999999993, 0.5, 0.0]

    # show data in view
    clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    clip1Display.Representation = 'Surface'
    clip1Display.ColorArrayName = [None, '']
    clip1Display.SelectTCoordArray = 'None'
    clip1Display.SelectNormalArray = 'None'
    clip1Display.SelectTangentArray = 'None'
    clip1Display.OSPRayScaleArray = 'ActiveNodes'
    clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip1Display.SelectOrientationVectors = 'None'
    clip1Display.ScaleFactor = 0.29
    clip1Display.SelectScaleArray = 'ActiveNodes'
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'ActiveNodes'
    clip1Display.GaussianRadius = 0.014499999999999997
    clip1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
    clip1Display.OpacityArray = ['POINTS', 'ActiveNodes']
    clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
    clip1Display.DataAxesGrid = 'GridAxesRepresentation'
    clip1Display.PolarAxes = 'PolarAxesRepresentation'
    clip1Display.ScalarOpacityUnitDistance = 0.2571277303864475
    clip1Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
    clip1Display.SelectInputVectors = [None, '']
    clip1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    clip1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip1.ClipType)

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=clip1)
    contour1.ContourBy = ['POINTS', 'ActiveNodes']
    contour1.Isosurfaces = [0.5]
    contour1.PointMergeMethod = 'Uniform Binning'

    # Properties modified on contour1
    contour1.ContourBy = ['POINTS', 'T']
    contour1.Isosurfaces = contourValues

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # get 2D transfer function for 'T'
    tTF2D = GetTransferFunction2D('T')
    tTF2D.ScalarRangeInitialized = 1
    tTF2D.Range = [25.0, 2400.0, 0.0, 1.0]

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')
    tLUT.AutomaticRescaleRangeMode = 'Never'
    tLUT.TransferFunction2D = tTF2D
    tLUT.RGBPoints = [25.0, 0.02, 0.3813, 0.9981, 81.54761904761904, 0.02000006, 0.424267768, 0.96906969, 138.09523809523807, 0.02, 0.467233763, 0.940033043, 194.64285714285714, 0.02, 0.5102, 0.911, 251.19047619047618, 0.02000006, 0.546401494, 0.872669438, 307.73809523809524, 0.02, 0.582600362, 0.83433295, 364.2857142857143, 0.02, 0.6188, 0.796, 420.8333333333333, 0.02000006, 0.652535156, 0.749802434, 477.38095238095235, 0.02, 0.686267004, 0.703599538, 533.9285714285713, 0.02, 0.72, 0.6574, 590.4761904761905, 0.02000006, 0.757035456, 0.603735359, 647.0238095238095, 0.02, 0.794067037, 0.55006613, 703.5714285714286, 0.02, 0.8311, 0.4964, 760.1190476190476, 0.021354336738172372, 0.8645368555261631, 0.4285579460761159, 816.6666666666666, 0.023312914349117714, 0.897999359924484, 0.36073871343115577, 873.2142857142858, 0.015976108242848862, 0.9310479513349017, 0.2925631815088092, 929.7619047619047, 0.27421074700988196, 0.952562960995083, 0.15356836602739213, 986.3095238095239, 0.4933546281681699, 0.9619038625309482, 0.11119493614749336, 1042.8571428571427, 0.6439, 0.9773, 0.0469, 1099.404761904762, 0.762401813, 0.984669591, 0.034600153, 1155.952380952381, 0.880901185, 0.992033407, 0.022299877, 1212.5, 0.9995285432627147, 0.9995193706781492, 0.0134884641450013, 1269.047619047619, 0.999402998, 0.955036376, 0.079066628, 1325.5952380952383, 0.9994, 0.910666223, 0.148134024, 1382.142857142857, 0.9994, 0.8663, 0.2172, 1438.6904761904761, 0.999269665, 0.818035981, 0.217200652, 1495.2380952380952, 0.999133332, 0.769766184, 0.2172, 1551.7857142857144, 0.999, 0.7215, 0.2172, 1608.3333333333333, 0.99913633, 0.673435546, 0.217200652, 1664.8809523809523, 0.999266668, 0.625366186, 0.2172, 1721.4285714285716, 0.9994, 0.5773, 0.2172, 1777.9761904761906, 0.999402998, 0.521068455, 0.217200652, 1834.5238095238094, 0.9994, 0.464832771, 0.2172, 1891.0714285714284, 0.9994, 0.4086, 0.2172, 1947.6190476190477, 0.9947599917687346, 0.33177297300202935, 0.2112309638520206, 2004.1666666666667, 0.9867129505479589, 0.2595183410914934, 0.19012239549291934, 2060.7142857142853, 0.9912458875646419, 0.14799417507952672, 0.21078892136920357, 2117.2619047619046, 0.949903037, 0.116867171, 0.252900603, 2173.809523809524, 0.903199533, 0.078432949, 0.291800389, 2230.357142857143, 0.8565, 0.04, 0.3307, 2286.904761904762, 0.798902627, 0.04333345, 0.358434298, 2343.4523809523807, 0.741299424, 0.0466667, 0.386166944, 2400.0, 0.6837, 0.05, 0.4139]
    tLUT.ColorSpace = 'RGB'
    tLUT.NanColor = [1.0, 0.0, 0.0]
    tLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = ['POINTS', 'T']
    contour1Display.LookupTable = tLUT
    contour1Display.SelectTCoordArray = 'None'
    contour1Display.SelectNormalArray = 'None'
    contour1Display.SelectTangentArray = 'None'
    contour1Display.OSPRayScaleArray = 'T'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 0.10342125543781429
    contour1Display.SelectScaleArray = 'T'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'T'
    contour1Display.GaussianRadius = 0.005171062771890714
    contour1Display.SetScaleArray = ['POINTS', 'T']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'T']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'
    contour1Display.SelectInputVectors = [None, '']
    contour1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    contour1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')
    tPWF.Points = [25.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]
    tPWF.ScalarRangeInitialized = 1

    # set active source
    SetActiveSource(clip1)

    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=clip1)
    calculator1.Function = ''

    # Properties modified on calculator1
    calculator1.ResultArrayName = 'roundedT'
    calculator1.Function = 'floor(T*10)/10'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

    # get 2D transfer function for 'roundedT'
    roundedTTF2D = GetTransferFunction2D('roundedT')

    # get color transfer function/color map for 'roundedT'
    roundedTLUT = GetColorTransferFunction('roundedT')
    roundedTLUT.TransferFunction2D = roundedTTF2D
    roundedTLUT.RGBPoints = [25.0, 1.0, 1.0, 1.0, 424.8060083010807, 0.0, 0.0, 1.0, 824.6120166021614, 0.0, 1.0, 1.0, 1200.9000244140625, 0.0, 1.0, 0.0, 1600.706032715143, 1.0, 1.0, 0.0, 2000.512041016224, 1.0, 0.0, 0.0, 2376.800048828125, 0.878431372549, 0.0, 1.0]
    roundedTLUT.ColorSpace = 'RGB'
    roundedTLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'roundedT'
    roundedTPWF = GetOpacityTransferFunction('roundedT')
    roundedTPWF.Points = [25.0, 0.0, 0.5, 0.0, 2376.800048828125, 1.0, 0.5, 0.0]
    roundedTPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'roundedT']
    calculator1Display.LookupTable = roundedTLUT
    calculator1Display.SelectTCoordArray = 'None'
    calculator1Display.SelectNormalArray = 'None'
    calculator1Display.SelectTangentArray = 'None'
    calculator1Display.OSPRayScaleArray = 'roundedT'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'None'
    calculator1Display.ScaleFactor = 0.29
    calculator1Display.SelectScaleArray = 'roundedT'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'roundedT'
    calculator1Display.GaussianRadius = 0.014499999999999997
    calculator1Display.SetScaleArray = ['POINTS', 'roundedT']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'roundedT']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'
    calculator1Display.ScalarOpacityFunction = roundedTPWF
    calculator1Display.ScalarOpacityUnitDistance = 0.2571277303864475
    calculator1Display.OpacityArrayName = ['POINTS', 'roundedT']
    calculator1Display.SelectInputVectors = [None, '']
    calculator1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    calculator1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [25.0, 0.0, 0.5, 0.0, 2376.3, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [25.0, 0.0, 0.5, 0.0, 2376.3, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(clip1, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on renderView1
    renderView1.CenterOfRotation = [-0.14999999478459358, 0.0494999997317791, 0.0]
    renderView1.CameraPosition = [-0.11521622491233216, 0.041634347699672225, 0.21346981293248826]
    renderView1.CameraFocalPoint = [-0.11521622491233216, 0.041634347699672225, 0.0]
    renderView1.CameraParallelScale = 0.2538725016795355

    # set active source
    SetActiveSource(contour1)

    # create a new 'Clip'
    clip2 = Clip(registrationName='Clip2', Input=contour1)
    clip2.ClipType = 'Plane'
    clip2.HyperTreeGridClipper = 'Plane'
    clip2.Scalars = ['POINTS', 'T']
    clip2.Value = 1500.0

    # init the 'Plane' selected for 'ClipType'
    clip2.ClipType.Origin = [-0.48421881822110746, 0.037545388657709525, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip2.HyperTreeGridClipper.Origin = [-0.48421881822110746, 0.037545388657709525, 0.0]

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip2.ClipType)

    # Properties modified on clip2.ClipType
    clip2.ClipType.Origin = [-0.05, 0.0, 0.0]
    clip2.ClipType.Normal = [-0.5, 1.0, 0.0]

    # show data in view
    clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    clip2Display.Representation = 'Surface'
    clip2Display.ColorArrayName = ['POINTS', 'T']
    clip2Display.LookupTable = tLUT
    clip2Display.SelectTCoordArray = 'None'
    clip2Display.SelectNormalArray = 'None'
    clip2Display.SelectTangentArray = 'None'
    clip2Display.OSPRayScaleArray = 'T'
    clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip2Display.SelectOrientationVectors = 'None'
    clip2Display.ScaleFactor = 0.008044822643816821
    clip2Display.SelectScaleArray = 'T'
    clip2Display.GlyphType = 'Arrow'
    clip2Display.GlyphTableIndexArray = 'T'
    clip2Display.GaussianRadius = 0.000402241132190841
    clip2Display.SetScaleArray = ['POINTS', 'T']
    clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
    clip2Display.OpacityArray = ['POINTS', 'T']
    clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
    clip2Display.DataAxesGrid = 'GridAxesRepresentation'
    clip2Display.PolarAxes = 'PolarAxesRepresentation'
    clip2Display.ScalarOpacityFunction = tPWF
    clip2Display.ScalarOpacityUnitDistance = 0.040650473364930315
    clip2Display.OpacityArrayName = ['POINTS', 'T']
    clip2Display.SelectInputVectors = [None, '']
    clip2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    clip2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    clip2Display.ScaleTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    clip2Display.OpacityTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(contour1, renderView1)

    # show color bar/color legend
    clip2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Slice'
    slice2 = Slice(registrationName='Slice2', Input=clip2)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice2.SliceType.Origin = [-0.007336654251120126, 0.01730698788025198, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice2.HyperTreeGridSlicer.Origin = [-0.007336654251120126, 0.01730698788025198, 0.0]

    # Properties modified on slice2.SliceType
    slice2.SliceType.Origin = [-0.05, 0.0, 0.0]
    slice2.SliceType.Normal = [-1.0, -0.5, 0.0]

    # show data in view
    slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice2Display.Representation = 'Surface'
    slice2Display.ColorArrayName = [None, '']
    slice2Display.SelectTCoordArray = 'None'
    slice2Display.SelectNormalArray = 'None'
    slice2Display.SelectTangentArray = 'None'
    slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice2Display.SelectOrientationVectors = 'None'
    slice2Display.ScaleFactor = -0.2
    slice2Display.SelectScaleArray = 'None'
    slice2Display.GlyphType = 'Arrow'
    slice2Display.GlyphTableIndexArray = 'None'
    slice2Display.GaussianRadius = -0.01
    slice2Display.SetScaleArray = [None, '']
    slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice2Display.OpacityArray = [None, '']
    slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice2Display.DataAxesGrid = 'GridAxesRepresentation'
    slice2Display.PolarAxes = 'PolarAxesRepresentation'
    slice2Display.SelectInputVectors = [None, '']
    slice2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    slice2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(clip2, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(clip2)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice2.SliceType)

    # Properties modified on clip2
    clip2.Invert = 0

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(slice2)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=slice2.SliceType)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice2.SliceType)

    # set active source
    SetActiveSource(contour1)

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    # turn off scalar coloring
    ColorBy(contour1Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(tLUT, renderView1)

    # change solid color
    contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

    # set active source
    SetActiveSource(calculator1)

    # create a query selection
    QuerySelect(QueryString='(roundedT == max(roundedT))', FieldType='POINT', InsideOut=0)

    # create a new 'Extract Selection'
    extractSelection1 = ExtractSelection(registrationName='ExtractSelection1', Input=calculator1)

    # show data in view
    extractSelection1Display = Show(extractSelection1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    extractSelection1Display.Representation = 'Surface'
    extractSelection1Display.ColorArrayName = ['POINTS', 'roundedT']
    extractSelection1Display.LookupTable = roundedTLUT
    extractSelection1Display.SelectTCoordArray = 'None'
    extractSelection1Display.SelectNormalArray = 'None'
    extractSelection1Display.SelectTangentArray = 'None'
    extractSelection1Display.OSPRayScaleArray = 'roundedT'
    extractSelection1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSelection1Display.SelectOrientationVectors = 'None'
    extractSelection1Display.ScaleFactor = 0.1
    extractSelection1Display.SelectScaleArray = 'roundedT'
    extractSelection1Display.GlyphType = 'Arrow'
    extractSelection1Display.GlyphTableIndexArray = 'roundedT'
    extractSelection1Display.GaussianRadius = 0.005
    extractSelection1Display.SetScaleArray = ['POINTS', 'roundedT']
    extractSelection1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSelection1Display.OpacityArray = ['POINTS', 'roundedT']
    extractSelection1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSelection1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSelection1Display.PolarAxes = 'PolarAxesRepresentation'
    extractSelection1Display.ScalarOpacityFunction = roundedTPWF
    extractSelection1Display.ScalarOpacityUnitDistance = 0.0
    extractSelection1Display.OpacityArrayName = ['POINTS', 'roundedT']
    extractSelection1Display.SelectInputVectors = [None, '']
    extractSelection1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    extractSelection1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSelection1Display.ScaleTransferFunction.Points = [2376.3, 0.0, 0.5, 0.0, 2376.800048828125, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSelection1Display.OpacityTransferFunction.Points = [2376.3, 0.0, 0.5, 0.0, 2376.800048828125, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(calculator1, renderView1)

    # show color bar/color legend
    extractSelection1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Annotate Attribute Data'
    annotateAttributeData1 = AnnotateAttributeData(registrationName='AnnotateAttributeData1', Input=extractSelection1)
    annotateAttributeData1.SelectInputArray = ['POINTS', 'roundedT']

    # Properties modified on annotateAttributeData1
    annotateAttributeData1.Prefix = ''

    # show data in view
    annotateAttributeData1Display = Show(annotateAttributeData1, renderView1, 'TextSourceRepresentation')

    # trace defaults for the display properties.
    annotateAttributeData1Display.WindowLocation = 'Any Location'
    annotateAttributeData1Display.FontSize = 34

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.Position = [0.55, 0.009375000000000001]

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.Position = [0.6, 0.009375000000000001]

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.Position = [0.575, 0.009375000000000001]

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.Position = [0.575, 0.001]

    # set active source
    SetActiveSource(clip2)

    # show data in view
    clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

    # show color bar/color legend
    clip2Display.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    Hide(clip2, renderView1)

    # set active source
    SetActiveSource(clip1)

    # show data in view
    clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')
    # set scalar coloring
    ColorBy(clip1Display, ('POINTS', 'T'))
    # rescale color and/or opacity maps used to include current data range
    clip1Display.RescaleTransferFunctionToDataRange(True, False)
    # show color bar/color legend
    clip1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(extractSelection1)

    # set active source
    SetActiveSource(annotateAttributeData1)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=annotateAttributeData1Display)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=annotateAttributeData1Display)

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.Position = [0.55, 0.001]

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.FontFamily = 'Times'
    annotateAttributeData1Display.Color = [0.0, 0.0, 0.0]

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.FontSize = 3

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.FontSize = 30

    # set active source
    SetActiveSource(extractSelection1)

    # turn off scalar coloring
    ColorBy(extractSelection1Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(roundedTLUT, renderView1)

    # change solid color
    extractSelection1Display.AmbientColor = [1.0, 0.0, 0.0]
    extractSelection1Display.DiffuseColor = [1.0, 0.0, 0.0]

    # Properties modified on extractSelection1Display
    extractSelection1Display.PointSize = 8.0

    # set active source
    SetActiveSource(slice2)

    # create a query selection
    QuerySelect(QueryString='(T <= 1e9)', FieldType='POINT', InsideOut=0)

    # Properties modified on slice2Display
    slice2Display.SelectionPointFieldDataArrayName = 'T'
    slice2Display.SelectionPointLabelVisibility = 1

    # Properties modified on slice2Display
    slice2Display.SelectionPointLabelColor = [0.0, 0.0, 0.16000610360875867]
    slice2Display.SelectionPointLabelFontFamily = 'Times'
    slice2Display.SelectionPointLabelFontSize = 22
    slice2Display.SelectionPointLabelFormat = '%.0f'

    if (closeInterface):
        #Hide(extractEdges1, renderView1)

        threshold1 = Threshold(registrationName='Threshold1', Input=extractEdges1)
        threshold1.Scalars = ['POINTS', 'gammaNodes']
        threshold1.LowerThreshold = 0.1
        threshold1.UpperThreshold = 1.0
        threshold1.UpdatePipeline()
        threshold1Display = Show(threshold1, renderView1)
        text1 = Text(registrationName='text1', Text = '$\\Gamma$',)
        text1Display = Show(text1, renderView1, 'TextSourceRepresentation')
        threshold1Display.DiffuseColor = [0.0, 0.0, 0.0]
        threshold1Display.LineWidth = 2.0
        text1Display.FontSize = 30
        text1Display.FontFamily = 'Times'
        text1Display.Position = [0.145, 0.018]

    # get layout
    layout1 = GetLayout()

    if layout1 == None:
        layout1 = CreateLayout()

    # layout/tab size in pixels
    layout1.SetSize(1586, 320)

    # color bar
    # get 2D transfer function for 'T'
    SetActiveSource(clip1)

    # Rescale transfer function
    # get 2D transfer function for 'T'
    tTF2D = GetTransferFunction2D('T')
    tTF2D.ScalarRangeInitialized = 1
    tTF2D.Range = [25.0, 2400.0, 0.0, 1.0]

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')
    tLUT.AutomaticRescaleRangeMode = 'Never'
    tLUT.TransferFunction2D = tTF2D
    tLUT.RGBPoints = [25.0, 0.02, 0.3813, 0.9981, 81.54761904761904, 0.02000006, 0.424267768, 0.96906969, 138.09523809523807, 0.02, 0.467233763, 0.940033043, 194.64285714285714, 0.02, 0.5102, 0.911, 251.19047619047618, 0.02000006, 0.546401494, 0.872669438, 307.73809523809524, 0.02, 0.582600362, 0.83433295, 364.2857142857143, 0.02, 0.6188, 0.796, 420.8333333333333, 0.02000006, 0.652535156, 0.749802434, 477.38095238095235, 0.02, 0.686267004, 0.703599538, 533.9285714285713, 0.02, 0.72, 0.6574, 590.4761904761905, 0.02000006, 0.757035456, 0.603735359, 647.0238095238095, 0.02, 0.794067037, 0.55006613, 703.5714285714286, 0.02, 0.8311, 0.4964, 760.1190476190476, 0.021354336738172372, 0.8645368555261631, 0.4285579460761159, 816.6666666666666, 0.023312914349117714, 0.897999359924484, 0.36073871343115577, 873.2142857142858, 0.015976108242848862, 0.9310479513349017, 0.2925631815088092, 929.7619047619047, 0.27421074700988196, 0.952562960995083, 0.15356836602739213, 986.3095238095239, 0.4933546281681699, 0.9619038625309482, 0.11119493614749336, 1042.8571428571427, 0.6439, 0.9773, 0.0469, 1099.404761904762, 0.762401813, 0.984669591, 0.034600153, 1155.952380952381, 0.880901185, 0.992033407, 0.022299877, 1212.5, 0.9995285432627147, 0.9995193706781492, 0.0134884641450013, 1269.047619047619, 0.999402998, 0.955036376, 0.079066628, 1325.5952380952383, 0.9994, 0.910666223, 0.148134024, 1382.142857142857, 0.9994, 0.8663, 0.2172, 1438.6904761904761, 0.999269665, 0.818035981, 0.217200652, 1495.2380952380952, 0.999133332, 0.769766184, 0.2172, 1551.7857142857144, 0.999, 0.7215, 0.2172, 1608.3333333333333, 0.99913633, 0.673435546, 0.217200652, 1664.8809523809523, 0.999266668, 0.625366186, 0.2172, 1721.4285714285716, 0.9994, 0.5773, 0.2172, 1777.9761904761906, 0.999402998, 0.521068455, 0.217200652, 1834.5238095238094, 0.9994, 0.464832771, 0.2172, 1891.0714285714284, 0.9994, 0.4086, 0.2172, 1947.6190476190477, 0.9947599917687346, 0.33177297300202935, 0.2112309638520206, 2004.1666666666667, 0.9867129505479589, 0.2595183410914934, 0.19012239549291934, 2060.7142857142853, 0.9912458875646419, 0.14799417507952672, 0.21078892136920357, 2117.2619047619046, 0.949903037, 0.116867171, 0.252900603, 2173.809523809524, 0.903199533, 0.078432949, 0.291800389, 2230.357142857143, 0.8565, 0.04, 0.3307, 2286.904761904762, 0.798902627, 0.04333345, 0.358434298, 2343.4523809523807, 0.741299424, 0.0466667, 0.386166944, 2400.0, 0.6837, 0.05, 0.4139]
    tLUT.ColorSpace = 'RGB'
    tLUT.NanColor = [1.0, 0.0, 0.0]
    tLUT.ScalarRangeInitialized = 1.0

    # get color legend/bar for tLUT in view renderView1
    tLUTColorBar = GetScalarBar(tLUT, renderView1)
    tLUTColorBar.Orientation = 'Vertical'
    tLUTColorBar.WindowLocation = 'Any Location'
    tLUTColorBar.Title = 'T'
    tLUTColorBar.ComponentTitle = ''
    tLUTColorBar.TitleColor = [0.0, 0.0, 0.16]
    tLUTColorBar.TitleFontFamily = 'Times'
    tLUTColorBar.TitleFontSize = 30
    tLUTColorBar.LabelColor = [0.0, 0.0, 0.16000610360875867]
    tLUTColorBar.LabelFontFamily = 'Times'
    tLUTColorBar.LabelFontSize = 30
    tLUTColorBar.ScalarBarLength = 0.6081250000000009
    tLUTColorBar.RangeLabelFormat = '%#.0f'
    # change scalar bar placement
    tLUTColorBar.Position = [0.9013493064312735, 0.1906250000000001]
    renderView1.OrientationAxesVisibility = 0

    # save screenshot
    SaveScreenshot('/home/mslimani/acuario/moving_heat_source/run/3d_slob/figures/' + fileName.split(".")[0] + ".png", renderView1, ImageResolution=[1586, 320],
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

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dataSet')
    parser.add_argument('-i', '--interface', action='store_true')
    parser.add_argument('--other-mesh-file', default="")
    args = parser.parse_args()
    takeScreenshot( args.dataSet, args.interface, args.other_mesh_file )
