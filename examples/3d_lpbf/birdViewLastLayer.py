from paraview.simple import *
import argparse

contours = [1600, 1800, 2000, 2200, 2400]

def main(dataSet, isCoupled):
    dataSetName = dataSet.split(".")[0]
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'PVD Reader'
    coupled_safepvd = PVDReader(registrationName='coupled_safe.pvd', FileName=dataSet)
    coupled_safepvd.CellArrays = ['ActiveElements', 'material', 'physicalDomain']
    coupled_safepvd.PointArrays = ['T', 'Pulse', 'ActiveNodes', 'gammaNodes', 'forcedDofs']

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    coupled_safepvdDisplay = Show(coupled_safepvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    coupled_safepvdDisplay.Representation = 'Surface'
    coupled_safepvdDisplay.ColorArrayName = [None, '']
    coupled_safepvdDisplay.SelectTCoordArray = 'None'
    coupled_safepvdDisplay.SelectNormalArray = 'None'
    coupled_safepvdDisplay.SelectTangentArray = 'None'
    coupled_safepvdDisplay.OSPRayScaleArray = 'ActiveNodes'
    coupled_safepvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    coupled_safepvdDisplay.SelectOrientationVectors = 'None'
    coupled_safepvdDisplay.ScaleFactor = 1.2000000000000002
    coupled_safepvdDisplay.SelectScaleArray = 'ActiveNodes'
    coupled_safepvdDisplay.GlyphType = 'Arrow'
    coupled_safepvdDisplay.GlyphTableIndexArray = 'ActiveNodes'
    coupled_safepvdDisplay.GaussianRadius = 0.06
    coupled_safepvdDisplay.SetScaleArray = ['POINTS', 'ActiveNodes']
    coupled_safepvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    coupled_safepvdDisplay.OpacityArray = ['POINTS', 'ActiveNodes']
    coupled_safepvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    coupled_safepvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    coupled_safepvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    coupled_safepvdDisplay.ScalarOpacityUnitDistance = 0.1839324947367901
    coupled_safepvdDisplay.OpacityArrayName = ['POINTS', 'ActiveNodes']
    coupled_safepvdDisplay.SelectInputVectors = [None, '']
    coupled_safepvdDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    coupled_safepvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    animationScene1.GoToLast()

    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=coupled_safepvd)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.0, 0.0, -0.375]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, -0.375]

    # Properties modified on slice1
    slice1.Triangulatetheslice = 0

    # Properties modified on slice1.SliceType
    slice1.SliceType.Origin = [0.0, 0.0, 0.25]
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = [None, '']
    slice1Display.SelectTCoordArray = 'None'
    slice1Display.SelectNormalArray = 'None'
    slice1Display.SelectTangentArray = 'None'
    slice1Display.OSPRayScaleArray = 'ActiveNodes'
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'None'
    slice1Display.ScaleFactor = 1.2000000000000002
    slice1Display.SelectScaleArray = 'ActiveNodes'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'ActiveNodes'
    slice1Display.GaussianRadius = 0.06
    slice1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = ['POINTS', 'ActiveNodes']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.PolarAxes = 'PolarAxesRepresentation'
    slice1Display.SelectInputVectors = [None, '']
    slice1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    slice1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(coupled_safepvd, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice1.SliceType)

    # Properties modified on renderView1
    renderView1.CenterOfRotation = [0.0, 0.0, 0.25]
    renderView1.CameraPosition = [-4.537542623855386, -0.013677557920216403, 0.6791527055335171]
    renderView1.CameraFocalPoint = [-4.537542623855386, -0.013677557920216403, 0.25]
    renderView1.CameraParallelScale = 6.082762530298219

    # set scalar coloring
    ColorBy(slice1Display, ('POINTS', 'T'))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # get 2D transfer function for 'T'
    tTF2D = GetTransferFunction2D('T')

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')
    tLUT.TransferFunction2D = tTF2D
    tLUT.RGBPoints = [55.277515103665884, 1.0, 1.0, 1.0, 468.1055581981267, 0.0, 0.0, 1.0, 880.9336012925875, 0.0, 1.0, 1.0, 1269.4776418520798, 0.0, 1.0, 0.0, 1682.3056849465406, 1.0, 1.0, 0.0, 2095.133728041001, 1.0, 0.0, 0.0, 2483.677768600494, 0.878431372549, 0.0, 1.0]
    tLUT.ColorSpace = 'RGB'
    tLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')
    tPWF.Points = [55.277515103665884, 0.0, 0.5, 0.0, 2483.677768600494, 1.0, 0.5, 0.0]
    tPWF.ScalarRangeInitialized = 1

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    tLUT.ApplyPreset('Rainbow Uniform', True)

    # create a query selection
    QuerySelect(QueryString='(T == max(T))', FieldType='POINT', InsideOut=0)

    # create a new 'Extract Selection'
    extractSelection1 = ExtractSelection(registrationName='ExtractSelection1', Input=slice1)

    # show data in view
    extractSelection1Display = Show(extractSelection1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    extractSelection1Display.Representation = 'Surface'
    extractSelection1Display.ColorArrayName = ['POINTS', 'T']
    extractSelection1Display.LookupTable = tLUT
    extractSelection1Display.SelectTCoordArray = 'None'
    extractSelection1Display.SelectNormalArray = 'None'
    extractSelection1Display.SelectTangentArray = 'None'
    extractSelection1Display.OSPRayScaleArray = 'ActiveNodes'
    extractSelection1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSelection1Display.SelectOrientationVectors = 'None'
    extractSelection1Display.ScaleFactor = 0.1
    extractSelection1Display.SelectScaleArray = 'ActiveNodes'
    extractSelection1Display.GlyphType = 'Arrow'
    extractSelection1Display.GlyphTableIndexArray = 'ActiveNodes'
    extractSelection1Display.GaussianRadius = 0.005
    extractSelection1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    extractSelection1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSelection1Display.OpacityArray = ['POINTS', 'ActiveNodes']
    extractSelection1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSelection1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSelection1Display.PolarAxes = 'PolarAxesRepresentation'
    extractSelection1Display.ScalarOpacityFunction = tPWF
    extractSelection1Display.ScalarOpacityUnitDistance = 0.0
    extractSelection1Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
    extractSelection1Display.SelectInputVectors = [None, '']
    extractSelection1Display.WriteLog = ''

    extractSelection1Display.AmbientColor = [0.0, 0.0, 0.0]
    extractSelection1Display.DiffuseColor = [0.0, 0.0, 0.0]

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    extractSelection1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSelection1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSelection1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(slice1, renderView1)

    # show color bar/color legend
    extractSelection1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Annotate Attribute Data'
    annotateAttributeData1 = AnnotateAttributeData(registrationName='AnnotateAttributeData1', Input=extractSelection1)
    annotateAttributeData1.SelectInputArray = ['POINTS', 'ActiveNodes']

    # set active source
    SetActiveSource(extractSelection1)

    # destroy annotateAttributeData1
    Delete(annotateAttributeData1)
    del annotateAttributeData1

    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=extractSelection1)
    calculator1.Function = ''

    # Properties modified on calculator1
    calculator1.ResultArrayName = 'roundedT'
    calculator1.Function = 'round(T*10)/10'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

    # get 2D transfer function for 'roundedT'
    roundedTTF2D = GetTransferFunction2D('roundedT')

    # get color transfer function/color map for 'roundedT'
    roundedTLUT = GetColorTransferFunction('roundedT')
    roundedTLUT.TransferFunction2D = roundedTTF2D
    roundedTLUT.RGBPoints = [2483.7, 1.0, 1.0, 1.0, 2483.7849916992186, 0.0, 0.0, 1.0, 2483.8699833984374, 0.0, 1.0, 1.0, 2483.9499755859374, 0.0, 1.0, 0.0, 2484.034967285156, 1.0, 1.0, 0.0, 2484.119958984375, 1.0, 0.0, 0.0, 2484.199951171875, 0.878431372549, 0.0, 1.0]
    roundedTLUT.ColorSpace = 'RGB'
    roundedTLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'roundedT'
    roundedTPWF = GetOpacityTransferFunction('roundedT')
    roundedTPWF.Points = [2483.7, 0.0, 0.5, 0.0, 2484.199951171875, 1.0, 0.5, 0.0]
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
    calculator1Display.ScaleFactor = 0.1
    calculator1Display.SelectScaleArray = 'roundedT'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'roundedT'
    calculator1Display.GaussianRadius = 0.005
    calculator1Display.SetScaleArray = ['POINTS', 'roundedT']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'roundedT']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'
    calculator1Display.ScalarOpacityFunction = roundedTPWF
    calculator1Display.ScalarOpacityUnitDistance = 0.0
    calculator1Display.OpacityArrayName = ['POINTS', 'roundedT']
    calculator1Display.SelectInputVectors = [None, '']
    calculator1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    calculator1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    calculator1Display.ScaleTransferFunction.Points = [2483.7, 0.0, 0.5, 0.0, 2484.199951171875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [2483.7, 0.0, 0.5, 0.0, 2484.199951171875, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(extractSelection1, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Annotate Attribute Data'
    annotateAttributeData1 = AnnotateAttributeData(registrationName='AnnotateAttributeData1', Input=calculator1)
    annotateAttributeData1.SelectInputArray = ['POINTS', 'roundedT']

    # Properties modified on annotateAttributeData1
    annotateAttributeData1.Prefix = ''

    # show data in view
    annotateAttributeData1Display = Show(annotateAttributeData1, renderView1, 'TextSourceRepresentation')

    # trace defaults for the display properties.
    annotateAttributeData1Display.WindowLocation = 'Any Location'

    # update the view to ensure updated data information
    renderView1.Update()

    # hide data in view
    Hide(calculator1, renderView1)

    # set active source
    SetActiveSource(extractSelection1)

    # show data in view
    extractSelection1Display = Show(extractSelection1, renderView1, 'UnstructuredGridRepresentation')

    # show color bar/color legend
    extractSelection1Display.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    Hide(extractSelection1, renderView1)

    # set active source
    SetActiveSource(slice1)

    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(extractSelection1)

    # show data in view
    extractSelection1Display = Show(extractSelection1, renderView1, 'UnstructuredGridRepresentation')

    # show color bar/color legend
    extractSelection1Display.SetScalarBarVisibility(renderView1, True)

    # turn off scalar coloring
    ColorBy(extractSelection1Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(tLUT, renderView1)

    # Properties modified on extractSelection1Display
    extractSelection1Display.PointSize = 8.0

    # set active source
    SetActiveSource(annotateAttributeData1)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=annotateAttributeData1Display)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=annotateAttributeData1Display)

    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.FontFamily = 'Times'

    # Properties modified on annotateAttributeData1Display
    # Properties modified on annotateAttributeData1Display
    annotateAttributeData1Display.FontSize = 40
    annotateAttributeData1Display.Position = [0.07, 0.57]
    annotateAttributeData1Display.Color = [1.0, 1.0, 1.0]
    #annotateAttributeData1Display.Bold = 1.0

    # set active source
    SetActiveSource(slice1)

    # Properties modified on renderView1
    renderView1.OrientationAxesVisibility = 0

    # get color legend/bar for tLUT in view renderView1
    tLUTColorBar = GetScalarBar(tLUT, renderView1)
    tLUTColorBar.WindowLocation = 'Any Location'
    tLUTColorBar.Title = 'T'
    tLUTColorBar.ComponentTitle = ''
    tLUTColorBar.TitleColor = [1.0, 1.0, 1.0]
    tLUTColorBar.TitleFontFamily = 'Times'
    tLUTColorBar.TitleFontSize = 38
    tLUTColorBar.LabelColor = [1.0, 1.0, 1.0]
    tLUTColorBar.LabelFontFamily = 'Times'
    tLUTColorBar.LabelFontSize = 38
    tLUTColorBar.ScalarBarLength = 0.33000000000000096
    tLUTColorBar.RangeLabelFormat = '%#.2f'

    # change scalar bar placement
    tLUTColorBar.Orientation = 'Horizontal'
    tLUTColorBar.Position = [0.620794044665012, 0.0]
    tLUTColorBar.ScalarBarLength = 0.33000000000000085

    # Rescale transfer function
    tLUT.RescaleTransferFunction(25.0, 2600.0)

    # Rescale transfer function
    tPWF.RescaleTransferFunction(25.0, 2600.0)

    # Rescale 2D transfer function
    tTF2D.RescaleTransferFunction(25.0, 2600.0, 0.0, 1.0)

    # Properties modified on tLUTColorBar
    tLUTColorBar.RangeLabelFormat = '%#.0f'

    # set active source
    SetActiveSource(extractSelection1)

    # set active source
    SetActiveSource(annotateAttributeData1)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=annotateAttributeData1Display)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=annotateAttributeData1Display)

    # set active source
    SetActiveSource(slice1)

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=slice1)
    contour1.ContourBy = ['POINTS', 'ActiveNodes']
    contour1.Isosurfaces = [0.5]
    contour1.PointMergeMethod = 'Uniform Binning'

    # Properties modified on contour1
    contour1.ContourBy = ['POINTS', 'T']
    contour1.Isosurfaces = contours

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

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
    contour1Display.ScaleFactor = 0.04223311054095902
    contour1Display.SelectScaleArray = 'T'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'T'
    contour1Display.GaussianRadius = 0.002111655527047951
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
    contour1Display.ScaleTransferFunction.Points = [1900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [1900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]

    contour1Display.Opacity = 0.5
    contour1Display.LineWidth = 3

    # hide data in view
    Hide(slice1, renderView1)

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # turn off scalar coloring
    ColorBy(contour1Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(tLUT, renderView1)

    # set active source
    SetActiveSource(slice1)

    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(slice1)

    # set active source
    SetActiveSource(contour1)

    # create a new 'Clip'
    clip1 = Clip(registrationName='Clip1', Input=contour1)
    clip1.ClipType = 'Plane'
    clip1.HyperTreeGridClipper = 'Plane'
    clip1.Scalars = ['POINTS', 'T']
    clip1.Value = 2150.0

    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [-4.770246272595661, -9.158908908610264e-06, 0.25]

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip1.HyperTreeGridClipper.Origin = [-4.770246272595661, -9.158908908610264e-06, 0.25]

    # Properties modified on clip1
    clip1.Invert = 0

    # Properties modified on clip1.ClipType
    clip1.ClipType.Origin = [-5.0, 0.0, 0.0]

    # show data in view
    clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    clip1Display.Representation = 'Surface'
    clip1Display.ColorArrayName = ['POINTS', 'T']
    clip1Display.LookupTable = tLUT
    clip1Display.SelectTCoordArray = 'None'
    clip1Display.SelectNormalArray = 'None'
    clip1Display.SelectTangentArray = 'None'
    clip1Display.OSPRayScaleArray = 'T'
    clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip1Display.SelectOrientationVectors = 'None'
    clip1Display.ScaleFactor = 0.04223311054095902
    clip1Display.SelectScaleArray = 'T'
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'T'
    clip1Display.GaussianRadius = 0.002111655527047951
    clip1Display.SetScaleArray = ['POINTS', 'T']
    clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
    clip1Display.OpacityArray = ['POINTS', 'T']
    clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
    clip1Display.DataAxesGrid = 'GridAxesRepresentation'
    clip1Display.PolarAxes = 'PolarAxesRepresentation'
    clip1Display.ScalarOpacityFunction = tPWF
    clip1Display.ScalarOpacityUnitDistance = 0.10434854781474205
    clip1Display.OpacityArrayName = ['POINTS', 'T']
    clip1Display.SelectInputVectors = [None, '']
    clip1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    clip1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    clip1Display.ScaleTransferFunction.Points = [1900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    clip1Display.OpacityTransferFunction.Points = [1900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(contour1, renderView1)

    # show color bar/color legend
    clip1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on clip1.ClipType
    clip1.ClipType.Origin = [-4.9, 0.0, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip1.ClipType)

    # create a new 'Slice'
    slice2 = Slice(registrationName='Slice2', Input=clip1)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice2.SliceType.Origin = [-4.7295403599454335, -9.158908908610264e-06, 0.25]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice2.HyperTreeGridSlicer.Origin = [-4.7295403599454335, -9.158908908610264e-06, 0.25]

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice2.SliceType)

    # Properties modified on slice2.SliceType
    slice2.SliceType.Origin = [0.0, 0.0, 0.0]
    slice2.SliceType.Normal = [0.0, 1.0, 0.0]

    # show data in view
    slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice2Display.Representation = 'Surface'
    slice2Display.ColorArrayName = ['POINTS', 'T']
    slice2Display.LookupTable = tLUT
    slice2Display.SelectTCoordArray = 'None'
    slice2Display.SelectNormalArray = 'None'
    slice2Display.SelectTangentArray = 'None'
    slice2Display.OSPRayScaleArray = 'T'
    slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice2Display.SelectOrientationVectors = 'None'
    slice2Display.ScaleFactor = 0.030884002364524845
    slice2Display.SelectScaleArray = 'T'
    slice2Display.GlyphType = 'Arrow'
    slice2Display.GlyphTableIndexArray = 'T'
    slice2Display.GaussianRadius = 0.001544200118226242
    slice2Display.SetScaleArray = ['POINTS', 'T']
    slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice2Display.OpacityArray = ['POINTS', 'T']
    slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice2Display.DataAxesGrid = 'GridAxesRepresentation'
    slice2Display.PolarAxes = 'PolarAxesRepresentation'
    slice2Display.SelectInputVectors = [None, '']
    slice2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    slice2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    slice2Display.ScaleTransferFunction.Points = [1900.0, 0.0, 0.5, 0.0, 2400.0000000000005, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice2Display.OpacityTransferFunction.Points = [1900.0, 0.0, 0.5, 0.0, 2400.0000000000005, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(clip1, renderView1)

    # show color bar/color legend
    slice2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a query selection
    QuerySelect(QueryString='(id >= 0)', FieldType='POINT', InsideOut=0)

    # create a query selection
    QuerySelect(QueryString='(id >= 0)', FieldType='POINT', InsideOut=0)

    # Properties modified on slice2Display
    slice2Display.SelectionPointLabelVisibility = 1
    slice2Display.SelectionPointLabelFontSize = 29
    slice2Display.SelectionPointLabelBold = 1

    # Properties modified on slice2Display
    slice2Display.SelectionPointLabelColor = [0.0, 0.0, 0.16000610360875867]
    slice2Display.SelectionPointLabelFontFamily = 'Times'
    slice2Display.SelectionPointLabelFormat = '%.0f'

    # set active source
    SetActiveSource(contour1)

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # change solid color
    contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

    # set active source
    SetActiveSource(contour1)

    # set active source
    SetActiveSource(slice1)

    # create a new 'Extract Edges'
    extractEdges1 = ExtractEdges(registrationName='ExtractEdges1', Input=slice1)
    extractEdges1Display = Show(extractEdges1, renderView1, 'GeometryRepresentation')
    extractEdges1Display.Representation = 'Surface'
    extractEdges1Display.ColorArrayName = [None, '']
    extractEdges1Display.Opacity = 0.2
    extractEdges1Display.DiffuseColor = [0.0, 0.0, 0.0]

    if isCoupled:
        # update the view to ensure updated data information
        renderView1.Update()

        # create a new 'Threshold'
        threshold1 = Threshold(registrationName='Threshold1', Input=extractEdges1)
        threshold1.Scalars = ['POINTS', 'gammaNodes']
        threshold1.LowerThreshold = 0.1
        threshold1.UpperThreshold = 1.0


        # show data in view
        threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

        # trace defaults for the display properties.
        threshold1Display.Representation = 'Surface'
        threshold1Display.ColorArrayName = ['POINTS', 'T']
        threshold1Display.LookupTable = tLUT
        threshold1Display.SelectTCoordArray = 'None'
        threshold1Display.SelectNormalArray = 'None'
        threshold1Display.SelectTangentArray = 'None'
        threshold1Display.OSPRayScaleArray = 'ActiveNodes'
        threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        threshold1Display.SelectOrientationVectors = 'None'
        threshold1Display.ScaleFactor = 0.05916317701339722
        threshold1Display.SelectScaleArray = 'ActiveNodes'
        threshold1Display.GlyphType = 'Arrow'
        threshold1Display.GlyphTableIndexArray = 'ActiveNodes'
        threshold1Display.GaussianRadius = 0.0029581588506698607
        threshold1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
        threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
        threshold1Display.OpacityArray = ['POINTS', 'ActiveNodes']
        threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
        threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
        threshold1Display.PolarAxes = 'PolarAxesRepresentation'
        threshold1Display.ScalarOpacityFunction = tPWF
        threshold1Display.ScalarOpacityUnitDistance = 0.19850482381536444
        threshold1Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
        threshold1Display.SelectInputVectors = [None, '']
        threshold1Display.WriteLog = ''
        threshold1Display.LineWidth = 2

        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
        threshold1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        threshold1Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        threshold1Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

        # show color bar/color legend
        threshold1Display.SetScalarBarVisibility(renderView1, True)

        # update the view to ensure updated data information
        renderView1.Update()

        # turn off scalar coloring
        ColorBy(threshold1Display, None)

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(tLUT, renderView1)

        # change solid color
        threshold1Display.DiffuseColor = [0.0, 0.0, 0.0]

        # create a new 'Text'
        text1 = Text(registrationName='Text1')

        # Properties modified on text1
        text1.Text = '$\\Gamma$'

        # show data in view
        text1Display = Show(text1, renderView1, 'TextSourceRepresentation')

        # trace defaults for the display properties.
        text1Display.WindowLocation = 'Any Location'
        text1Display.Color = [1.0, 1.0, 1.0]
        text1Display.Position = [0.34, 0.1]

        # update the view to ensure updated data information
        renderView1.Update()

        # Properties modified on text1Display
        text1Display.FontFamily = 'Times'

        # Properties modified on text1Display
        text1Display.FontSize = 40
        text1Display.Bold = 1

    # set active source
    SetActiveSource(slice1)

    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(slice1)


    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1612, 368)

    # current camera placement for renderView1
    renderView1.CameraPosition = [-4.537542623855386, -0.013677557920216403, 0.6791527055335171]
    renderView1.CameraFocalPoint = [-4.537542623855386, -0.013677557920216403, 0.25]
    renderView1.CameraParallelScale = 6.082762530298219

    # save screenshot
    SaveScreenshot('./figures/birdViewLastTstep_{}.png'.format( dataSetName ), renderView1, ImageResolution=[1612, 368],
        TransparentBackground=1, 
        # PNG options
        CompressionLevel='3')

    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================

    #--------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(1612, 368)

    #-----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [-4.537542623855386, -0.013677557920216403, 0.6791527055335171]
    renderView1.CameraFocalPoint = [-4.537542623855386, -0.013677557920216403, 0.25]
    renderView1.CameraParallelScale = 6.082762530298219

    #--------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dataSet')
    parser.add_argument('-c', '--is-coupled', action='store_true')
    args = parser.parse_args()
    main( args.dataSet, args.is_coupled )
