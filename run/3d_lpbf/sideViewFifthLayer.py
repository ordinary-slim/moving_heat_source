from paraview.simple import *
import argparse

#contourValues = [900.0, 1150.0, 1400.0, 1700.0, 1900.0, 2200.0, 2400.0]
contourValues = [1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0]

def main(dataSet, isCoupled):
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
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
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
    # create a new 'Threshold'
    threshold1 = Threshold(registrationName='Threshold1', Input=coupled_safepvd)
    threshold1.Scalars = ['POINTS', 'ActiveNodes']
    threshold1.UpperThreshold = 1.0
    # Properties modified on threshold1
    thresholdBy = "ActiveElements"
    if isCoupled:
        thresholdBy = "physicalDomain"
    threshold1.Scalars = ['CELLS',thresholdBy]
    threshold1.LowerThreshold = 0.1
    # show data in view
    threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')
    # trace defaults for the display properties.
    threshold1Display.Representation = 'Surface'
    threshold1Display.ColorArrayName = [None, '']
    threshold1Display.SelectTCoordArray = 'None'
    threshold1Display.SelectNormalArray = 'None'
    threshold1Display.SelectTangentArray = 'None'
    threshold1Display.OSPRayScaleArray = 'ActiveNodes'
    threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold1Display.SelectOrientationVectors = 'None'
    threshold1Display.ScaleFactor = 1.2000000000000002
    threshold1Display.SelectScaleArray = 'ActiveNodes'
    threshold1Display.GlyphType = 'Arrow'
    threshold1Display.GlyphTableIndexArray = 'ActiveNodes'
    threshold1Display.GaussianRadius = 0.06
    threshold1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
    threshold1Display.OpacityArray = ['POINTS', 'ActiveNodes']
    threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
    threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold1Display.PolarAxes = 'PolarAxesRepresentation'
    threshold1Display.ScalarOpacityUnitDistance = 0.22204196676543572
    threshold1Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
    threshold1Display.SelectInputVectors = [None, '']
    threshold1Display.WriteLog = ''
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    threshold1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    threshold1Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    threshold1Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
    # hide data in view
    Hide(coupled_safepvd, renderView1)
    # update the view to ensure updated data information
    renderView1.Update()
    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=threshold1)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.0, 0.0, -0.5]
    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, -0.5]
    # Properties modified on slice1.SliceType
    slice1.SliceType.Origin = [0.0, 0.0, 0.0]
    slice1.SliceType.Normal = [0.0, 1.0, 0.0]
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
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    slice1Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice1Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
    # hide data in view
    Hide(threshold1, renderView1)
    # update the view to ensure updated data information
    renderView1.Update()
    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice1.SliceType)
    # Properties modified on animationScene1
    animationScene1.AnimationTime = 10.0484
    # Properties modified on renderView1
    renderView1.CenterOfRotation = [0.0, 0.0, 0.25]
    # Properties modified on renderView1
    renderView1.CenterOfRotation = [0.0, 0.0, -0.4375]
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
    tLUT.RGBPoints = [50.92100100541105, 1.0, 1.0, 1.0, 473.93415622284846, 0.0, 0.0, 1.0, 896.947311440286, 0.0, 1.0, 1.0, 1295.0773398802269, 0.0, 1.0, 0.0, 1718.0904950976644, 1.0, 1.0, 0.0, 2141.1036503151017, 1.0, 0.0, 0.0, 2539.233678755043, 0.878431372549, 0.0, 1.0]
    tLUT.ColorSpace = 'RGB'
    tLUT.ScalarRangeInitialized = 1.0
    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')
    tPWF.Points = [50.92100100541105, 0.0, 0.5, 0.0, 2539.233678755043, 1.0, 0.5, 0.0]
    tPWF.ScalarRangeInitialized = 1
    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=slice1)
    contour1.ContourBy = ['POINTS', 'ActiveNodes']
    contour1.Isosurfaces = [0.5]
    contour1.PointMergeMethod = 'Uniform Binning'
    # Properties modified on contour1
    contour1.ContourBy = ['POINTS', 'T']
    contour1.Isosurfaces = contourValues
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
    contour1Display.ScaleFactor = 0.31429170664548556
    contour1Display.SelectScaleArray = 'T'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'T'
    contour1Display.GaussianRadius = 0.015714585332274277
    contour1Display.SetScaleArray = ['POINTS', 'T']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'T']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'
    contour1Display.SelectInputVectors = [None, '']
    contour1Display.WriteLog = ''
    contour1Display.DiffuseColor = [0.0, 0.0, 0.0]
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    contour1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]
    contour1Display.LineWidth = 3
    contour1Display.Opacity = 0.5
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
    # create a new 'Clip'
    clip1 = Clip(registrationName='Clip1', Input=contour1)
    clip1.ClipType = 'Plane'
    clip1.HyperTreeGridClipper = 'Plane'
    clip1.Scalars = ['POINTS', 'ActiveNodes']
    clip1.Value = 0.5
    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [0.0, 0.0, -0.4375]
    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip1.HyperTreeGridClipper.Origin = [0.0, 0.0, -0.4375]
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    # Properties modified on clip1.ClipType
    clip1.ClipType.Origin = [3.225, 0.0, 0.0]
    clip1.ClipType.Normal = [-1.0, 0.0, 0.0]
    # show data in view
    clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')
    # trace defaults for the display properties.
    clip1Display.Representation = 'Surface'
    clip1Display.ColorArrayName = ['POINTS', 'T']
    clip1Display.LookupTable = tLUT
    clip1Display.SelectTCoordArray = 'None'
    clip1Display.SelectNormalArray = 'None'
    clip1Display.SelectTangentArray = 'None'
    clip1Display.OSPRayScaleArray = 'ActiveNodes'
    clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip1Display.SelectOrientationVectors = 'None'
    clip1Display.ScaleFactor = 0.2775
    clip1Display.SelectScaleArray = 'ActiveNodes'
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'ActiveNodes'
    clip1Display.GaussianRadius = 0.013875
    clip1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
    clip1Display.OpacityArray = ['POINTS', 'ActiveNodes']
    clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
    clip1Display.DataAxesGrid = 'GridAxesRepresentation'
    clip1Display.PolarAxes = 'PolarAxesRepresentation'
    clip1Display.ScalarOpacityFunction = tPWF
    clip1Display.ScalarOpacityUnitDistance = 0.20706728679881936
    clip1Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
    clip1Display.SelectInputVectors = [None, '']
    clip1Display.WriteLog = ''
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    clip1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]
    # hide data in view
    Hide(slice1, renderView1)
    # show color bar/color legend
    clip1Display.SetScalarBarVisibility(renderView1, True)
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on clip1
    clip1.Invert = 0
    # update the view to ensure updated data information
    renderView1.Update()
    # set active source
    SetActiveSource(slice1)
    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)
    # set active source
    SetActiveSource(clip1)
    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=clip1.ClipType)
    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip1.ClipType)
    # create a new 'Slice'
    slice2 = Slice(registrationName='Slice2', Input=clip1)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    # init the 'Plane' selected for 'SliceType'
    slice2.SliceType.Origin = [1.7034577076803987, 0.0, 0.07820883881408089]
    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice2.HyperTreeGridSlicer.Origin = [1.7034577076803987, 0.0, 0.07820883881408089]
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    renderView1.InteractionMode = 'Selection'
    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice2.SliceType)
    # Properties modified on slice2.SliceType
    slice2.SliceType.Origin = [3.225, 0.0, 0.125]
    slice2.SliceType.Normal = [1.0, 0.0, -2.0]
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
    slice2Display.ScaleFactor = 0.014019490651772727
    slice2Display.SelectScaleArray = 'T'
    slice2Display.GlyphType = 'Arrow'
    slice2Display.GlyphTableIndexArray = 'T'
    slice2Display.GaussianRadius = 0.0007009745325886363
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
    slice2Display.ScaleTransferFunction.Points = [900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice2Display.OpacityTransferFunction.Points = [900.0, 0.0, 0.5, 0.0, 2400.0, 1.0, 0.5, 0.0]
    # hide data in view
    Hide(clip1, renderView1)
    # show color bar/color legend
    slice2Display.SetScalarBarVisibility(renderView1, True)
    # update the view to ensure updated data information
    renderView1.Update()
    # create a query selection
    SetActiveSource( slice1 )
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
    # set active source
    SetActiveSource(slice1)
    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)
    # set active source
    SetActiveSource(slice1)
    # set active source
    SetActiveSource(slice2)
    # create a query selection
    QuerySelect(QueryString='(id >= 0)', FieldType='POINT', InsideOut=0)
    # create a query selection
    QuerySelect(QueryString='(id >= 0)', FieldType='POINT', InsideOut=0)
    # Properties modified on slice2.SliceType
    slice2.SliceType.Normal = [1.0, 0.0, -3.0]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2Display
    slice2Display.SelectionPointLabelVisibility = 1
    # Properties modified on slice2Display
    slice2Display.SelectionPointLabelBold = 1
    slice2Display.SelectionPointLabelFontFamily = 'Times'
    slice2Display.SelectionPointLabelFormat = '%.0f'
    slice2Display.SelectionPointLabelFontSize = 29
    slice2Display.SelectionPointLabelColor = [0.0, 0.0, 0.16000610360875867]
    # set active source
    SetActiveSource(slice2)
    # Properties modified on slice2.SliceType
    slice2.SliceType.Normal = [1.0, 0.0, -4.0]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Origin = [3.225, 0.0, 0.124]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Origin = [3.225, 0.0, 0.122]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Origin = [3.225, 0.0, 0.1]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Normal = [1.0, 0.0, -8.0]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Normal = [1.0, 0.0, -10.0]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Normal = [1.0, 0.0, -12.0]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Normal = [1.0, 0.0, -14.0]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2.SliceType
    slice2.SliceType.Normal = [1.0, 0.0, -16.0]
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on slice2Display
    slice2Display.SelectionPointLabelJustification = 'Right'
    # set active source
    SetActiveSource(extractSelection1)
    # turn off scalar coloring
    ColorBy(extractSelection1Display, None)
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(tLUT, renderView1)
    # Properties modified on extractSelection1Display
    extractSelection1Display.PointSize = 8.0
    # change solid color
    extractSelection1Display.DiffuseColor = [0.0, 0.0, 0.0]
    renderView1.Update()
    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=extractSelection1)
    calculator1.Function = ''
    # Properties modified on calculator1
    calculator1.ResultArrayName = 'roundedT'
    calculator1.Function = 'round(10*T)/10'
    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')
    SetActiveSource( calculator1 )
    # get 2D transfer function for 'roundedT'
    roundedTTF2D = GetTransferFunction2D('roundedT')
    # get color transfer function/color map for 'roundedT'
    roundedTLUT = GetColorTransferFunction('roundedT')
    roundedTLUT.TransferFunction2D = roundedTTF2D
    roundedTLUT.RGBPoints = [25.0, 1.0, 1.0, 1.0, 452.4989916992645, 0.0, 0.0, 1.0, 879.997983398529, 0.0, 1.0, 1.0, 1282.3499755859375, 0.0, 1.0, 0.0, 1709.8489672852018, 1.0, 1.0, 0.0, 2137.3479589844665, 1.0, 0.0, 0.0, 2539.699951171875, 0.878431372549, 0.0, 1.0]
    roundedTLUT.ColorSpace = 'RGB'
    roundedTLUT.ScalarRangeInitialized = 1.0
    # get opacity transfer function/opacity map for 'roundedT'
    roundedTPWF = GetOpacityTransferFunction('roundedT')
    roundedTPWF.Points = [25.0, 0.0, 0.5, 0.0, 2539.699951171875, 1.0, 0.5, 0.0]
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
    calculator1Display.ScaleTransferFunction.Points = [25.0, 0.0, 0.5, 0.0, 25.00390625, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [25.0, 0.0, 0.5, 0.0, 25.00390625, 1.0, 0.5, 0.0]
    # hide data in view
    Hide( calculator1, renderView1 )
    # update the view to ensure updated data information
    renderView1.Update()
    # create a new 'Annotate Attribute Data'
    annotateAttributeData1 = AnnotateAttributeData(registrationName='AnnotateAttributeData1', Input=calculator1)
    annotateAttributeData1.Prefix = ''
    annotateAttributeData1.SelectInputArray = ['POINTS', 'roundedT']
    # show data in view
    annotateAttributeData1Display = Show(annotateAttributeData1, renderView1, 'TextSourceRepresentation')
    # trace defaults for the display properties.
    annotateAttributeData1Display.WindowLocation = 'Any Location'
    annotateAttributeData1Display.FontSize = 40
    annotateAttributeData1Display.FontFamily = 'Times'
    annotateAttributeData1Display.Position = [0.87, 0.75]
    annotateAttributeData1Display.Color = [1.0, 1.0, 1.0]
    renderView1.Update()
    # set active source
    SetActiveSource(slice1)
    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    tLUT.ApplyPreset('Rainbow Uniform', True)
    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    tLUT.ApplyPreset('Rainbow Uniform', True)
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
    tLUTColorBar.Position = [0.3861329431438122, 0.25]
    tLUTColorBar.ScalarBarLength = 0.3300000000000011
    # Properties modified on tLUTColorBar
    tLUTColorBar.RangeLabelFormat = '%#.0f'
    # Rescale transfer function
    tLUT.RescaleTransferFunction(25.0, 2600.0)
    # Rescale transfer function
    tPWF.RescaleTransferFunction(25.0, 2600.0)
    # Rescale 2D transfer function
    tTF2D.RescaleTransferFunction(25.0, 2600.0, 0.0, 1.0)

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
        threshold2 = Threshold(registrationName='Threshold2', Input=extractEdges1)
        threshold2.Scalars = ['POINTS', 'ActiveNodes']
        threshold2.UpperThreshold = 1.0
        # Properties modified on threshold2
        threshold2.Scalars = ['POINTS', 'gammaNodes']
        threshold2.LowerThreshold = 0.1
        # show data in view
        threshold2Display = Show(threshold2, renderView1, 'UnstructuredGridRepresentation')
        # trace defaults for the display properties.
        threshold2Display.Representation = 'Surface'
        threshold2Display.ColorArrayName = ['POINTS', 'T']
        threshold2Display.LookupTable = tLUT
        threshold2Display.SelectTCoordArray = 'None'
        threshold2Display.SelectNormalArray = 'None'
        threshold2Display.SelectTangentArray = 'None'
        threshold2Display.OSPRayScaleArray = 'ActiveNodes'
        threshold2Display.OSPRayScaleFunction = 'PiecewiseFunction'
        threshold2Display.SelectOrientationVectors = 'None'
        threshold2Display.ScaleFactor = 0.5299999952316284
        threshold2Display.SelectScaleArray = 'ActiveNodes'
        threshold2Display.GlyphType = 'Arrow'
        threshold2Display.GlyphTableIndexArray = 'ActiveNodes'
        threshold2Display.GaussianRadius = 0.026499999761581423
        threshold2Display.SetScaleArray = ['POINTS', 'ActiveNodes']
        threshold2Display.ScaleTransferFunction = 'PiecewiseFunction'
        threshold2Display.OpacityArray = ['POINTS', 'ActiveNodes']
        threshold2Display.OpacityTransferFunction = 'PiecewiseFunction'
        threshold2Display.DataAxesGrid = 'GridAxesRepresentation'
        threshold2Display.PolarAxes = 'PolarAxesRepresentation'
        threshold2Display.ScalarOpacityFunction = tPWF
        threshold2Display.ScalarOpacityUnitDistance = 0.8615535710745331
        threshold2Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
        threshold2Display.SelectInputVectors = [None, '']
        threshold2Display.WriteLog = ''
        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
        threshold2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]
        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        threshold2Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        threshold2Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
        # show color bar/color legend
        threshold2Display.SetScalarBarVisibility(renderView1, True)
        # change solid color
        threshold2Display.DiffuseColor = [0.0, 0.0, 0.0]
        # Properties modified on threshold2Display
        threshold2Display.LineWidth = 2.0
        # turn off scalar coloring
        ColorBy(threshold2Display, None)
        # create a new 'Text'
        text1 = Text(registrationName='Text1')
        # Properties modified on text1
        text1.Text = '$\\Gamma$'
        # show data in view
        text1Display = Show(text1, renderView1, 'TextSourceRepresentation')
        # trace defaults for the display properties.
        text1Display.WindowLocation = 'Any Location'
        # update the view to ensure updated data information
        renderView1.Update()
        # set active source
        SetActiveSource(text1)
        # toggle interactive widget visibility (only when running from the GUI)
        ShowInteractiveWidgets(proxy=text1Display)
        # toggle interactive widget visibility (only when running from the GUI)
        HideInteractiveWidgets(proxy=text1Display)
        # Properties modified on text1Display
        text1Display.FontFamily = 'Times'
        # Properties modified on text1Display
        text1Display.FontSize = 40
        text1Display.Bold = 1
        # Properties modified on text1Display
        text1Display.Color = [1.0, 1.0, 1.0]
        text1Display.Position = [0.12, 0.07]

    # update the view to ensure updated data information
    renderView1.Update()
    # set active source
    SetActiveSource(slice1)
    # Properties modified on slice1
    slice1.Triangulatetheslice = 0
    # update the view to ensure updated data information
    renderView1.Update()
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(tLUT, renderView1)
    # set active source
    SetActiveSource(slice1)
    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)
    # change scalar bar placement
    # get layout
    layout1 = GetLayout()
    # layout/tab size in pixels
    layout1.SetSize(1612, 368)
    # saving layout sizes for layouts
    '''
    Old
    renderView1.CameraPosition = [2.754331707978817, -0.6224911282488886, -0.04289387845433411]
    renderView1.CameraFocalPoint = [2.754331707978817, 0.0, -0.04289387845433411]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 6.026309504995574
    '''
    renderView1.CameraPosition = [2.8533576368104205, -0.413066684355839, 0.014219890113180509]
    renderView1.CameraFocalPoint = [2.8533576368104205, 0.0, 0.014219890113180509]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 6.026309504995574
    # Properties modified on renderView1
    renderView1.OrientationAxesVisibility = 0
    renderView1.Update()

    # save screenshot
    SaveScreenshot('figures/sideView_fifthlayer_{}.png'.format(dataSet.split(".")[0]), renderView1, ImageResolution=[1612, 368],
        TransparentBackground=1)
    # set active source
    SetActiveSource(slice1)
    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================
    #--------------------------------
    # layout/tab size in pixels
    layout1.SetSize(1612, 368)
    #-----------------------------------
    # saving camera placements for views
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
