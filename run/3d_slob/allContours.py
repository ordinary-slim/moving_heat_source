# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

import argparse
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def setTime():
    animationScene1 = GetAnimationScene()
    animationScene1.AnimationTime = 0.0022
    renderView1 = GetActiveViewOrCreate('RenderView')
    renderView1.ViewTime = 0.0022
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    renderView1.Update()

def main(coupledDataSet, refDataSet):
    # create a new 'PVD Reader'
    analyticalpvd = PVDReader(registrationName='analytical.pvd', FileName='/home/mslimani/acuario/moving_heat_source/run/3d_slob/analyticalSeparateMesh/analytical.pvd')
    analyticalpvd.PointArrays = ['T']


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
    analyticalpvdDisplay.ScaleFactor = 0.06999999999999999
    analyticalpvdDisplay.SelectScaleArray = 'None'
    analyticalpvdDisplay.GlyphType = 'Arrow'
    analyticalpvdDisplay.GlyphTableIndexArray = 'None'
    analyticalpvdDisplay.GaussianRadius = 0.0034999999999999996
    analyticalpvdDisplay.SetScaleArray = ['POINTS', 'T']
    analyticalpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    analyticalpvdDisplay.OpacityArray = ['POINTS', 'T']
    analyticalpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    analyticalpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    analyticalpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    analyticalpvdDisplay.ScalarOpacityUnitDistance = 0.012895422730486595
    analyticalpvdDisplay.OpacityArrayName = ['POINTS', 'T']
    analyticalpvdDisplay.SelectInputVectors = [None, '']
    analyticalpvdDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    analyticalpvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    analyticalpvdDisplay.ScaleTransferFunction.Points = [25.0, 0.0, 0.5, 0.0, 2383.5577393095336, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    analyticalpvdDisplay.OpacityTransferFunction.Points = [25.0, 0.0, 0.5, 0.0, 2383.5577393095336, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Threshold'
    threshold1 = Threshold(registrationName='Threshold1', Input=analyticalpvd)
    threshold1.Scalars = ['POINTS', 'T']
    threshold1.LowerThreshold = 25.0
    threshold1.UpperThreshold = 2383.5577393095336

    # set active source
    SetActiveSource(analyticalpvd)

    # destroy threshold1
    Delete(threshold1)
    del threshold1

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=analyticalpvd)
    contour1.ContourBy = ['POINTS', 'T']
    contour1.Isosurfaces = [700.0, 1100.0, 1500.0, 1900.0, 2300.0, 2100.0]
    contour1.PointMergeMethod = 'Uniform Binning'

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # get 2D transfer function for 'T'
    tTF2D = GetTransferFunction2D('T')

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')
    tLUT.TransferFunction2D = tTF2D
    tLUT.RGBPoints = [700.0, 1.0, 1.0, 1.0, 972.0, 0.0, 0.0, 1.0, 1244.0, 0.0, 1.0, 1.0, 1500.0, 0.0, 1.0, 0.0, 1772.0, 1.0, 1.0, 0.0, 2043.9999999999998, 1.0, 0.0, 0.0, 2300.0, 0.878431372549, 0.0, 1.0]
    tLUT.ColorSpace = 'RGB'
    tLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = ['POINTS', 'T']
    contour1Display.LookupTable = tLUT
    contour1Display.SelectTCoordArray = 'None'
    contour1Display.SelectNormalArray = 'Normals'
    contour1Display.SelectTangentArray = 'None'
    contour1Display.OSPRayScaleArray = 'T'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 0.05330797110465644
    contour1Display.SelectScaleArray = 'T'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'T'
    contour1Display.GaussianRadius = 0.0026653985552328218
    contour1Display.SetScaleArray = ['POINTS', 'T']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'T']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'
    contour1Display.SelectInputVectors = ['POINTS', 'Normals']
    contour1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    contour1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(analyticalpvd, renderView1)

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')
    tPWF.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]
    tPWF.ScalarRangeInitialized = 1

    # set active source
    SetActiveSource(analyticalpvd)

    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=analyticalpvd)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.0, 0.0, 0.0]
    slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.0]
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = [None, '']
    slice1Display.SelectTCoordArray = 'None'
    slice1Display.SelectNormalArray = 'None'
    slice1Display.SelectTangentArray = 'None'
    slice1Display.OSPRayScaleArray = 'T'
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'None'
    slice1Display.ScaleFactor = 0.06999999999999999
    slice1Display.SelectScaleArray = 'None'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'None'
    slice1Display.GaussianRadius = 0.0034999999999999996
    slice1Display.SetScaleArray = ['POINTS', 'T']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = ['POINTS', 'T']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.PolarAxes = 'PolarAxesRepresentation'
    slice1Display.SelectInputVectors = [None, '']
    slice1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    slice1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    slice1Display.ScaleTransferFunction.Points = [25.000000000000064, 0.0, 0.5, 0.0, 2383.5577393095336, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice1Display.OpacityTransferFunction.Points = [25.000000000000064, 0.0, 0.5, 0.0, 2383.5577393095336, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(analyticalpvd, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on slice1
    slice1.Triangulatetheslice = 0

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(contour1)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice1.SliceType)

    # Properties modified on contour1
    contour1.Input = slice1

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice1.SliceType)

    # hide data in view
    Hide(slice1, renderView1)

    # create a new 'PVD Reader'
    reference_elsPerRad4_tstepsPerRad4pvd = PVDReader(registrationName='reference_elsPerRad4_tstepsPerRad4.pvd', FileName=refDataSet)
    reference_elsPerRad4_tstepsPerRad4pvd.CellArrays = ['ActiveElements']
    reference_elsPerRad4_tstepsPerRad4pvd.PointArrays = ['T', 'Pulse', 'ActiveNodes']

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # show data in view
    reference_elsPerRad4_tstepsPerRad4pvdDisplay = Show(reference_elsPerRad4_tstepsPerRad4pvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.Representation = 'Surface'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.ColorArrayName = [None, '']
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.SelectTCoordArray = 'None'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.SelectNormalArray = 'None'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.SelectTangentArray = 'None'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.OSPRayScaleArray = 'ActiveNodes'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.SelectOrientationVectors = 'None'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.ScaleFactor = 0.3600000000000001
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.SelectScaleArray = 'ActiveNodes'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.GlyphType = 'Arrow'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.GlyphTableIndexArray = 'ActiveNodes'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.GaussianRadius = 0.018000000000000002
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.SetScaleArray = ['POINTS', 'ActiveNodes']
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.OpacityArray = ['POINTS', 'ActiveNodes']
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.ScalarOpacityUnitDistance = 0.11913968780301962
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.OpacityArrayName = ['POINTS', 'ActiveNodes']
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.SelectInputVectors = [None, '']
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    reference_elsPerRad4_tstepsPerRad4pvdDisplay.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Slice'
    slice2 = Slice(registrationName='Slice2', Input=reference_elsPerRad4_tstepsPerRad4pvd)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice2.SliceType.Origin = [0.0, 0.0, 0.0]
    slice2.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.0]
    slice2.SliceType.Normal = [0.0, 0.0, 1.0]


    # set active source
    SetActiveSource(slice1)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice2.SliceType)

    # set active source
    SetActiveSource(slice2)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=slice2.SliceType)

    # show data in view
    slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice2Display.Representation = 'Surface'
    slice2Display.ColorArrayName = [None, '']
    slice2Display.SelectTCoordArray = 'None'
    slice2Display.SelectNormalArray = 'None'
    slice2Display.SelectTangentArray = 'None'
    slice2Display.OSPRayScaleArray = 'ActiveNodes'
    slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice2Display.SelectOrientationVectors = 'None'
    slice2Display.ScaleFactor = 0.3600000000000001
    slice2Display.SelectScaleArray = 'ActiveNodes'
    slice2Display.GlyphType = 'Arrow'
    slice2Display.GlyphTableIndexArray = 'ActiveNodes'
    slice2Display.GaussianRadius = 0.018000000000000002
    slice2Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice2Display.OpacityArray = ['POINTS', 'ActiveNodes']
    slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice2Display.DataAxesGrid = 'GridAxesRepresentation'
    slice2Display.PolarAxes = 'PolarAxesRepresentation'
    slice2Display.SelectInputVectors = [None, '']
    slice2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    slice2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    slice2Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice2Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(reference_elsPerRad4_tstepsPerRad4pvd, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice2.SliceType)

    # set active source
    SetActiveSource(contour1)

    # set active source
    SetActiveSource(slice2)

    # create a new 'Contour'
    contour2 = Contour(registrationName='Contour2', Input=slice2)
    contour2.ContourBy = ['POINTS', 'T']
    contour2.Isosurfaces = [700.0, 1100.0, 1500.0, 1900.0, 2300.0, 2100.0]
    contour2.PointMergeMethod = 'Uniform Binning'

    # show data in view
    contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    contour2Display.Representation = 'Surface'
    contour2Display.ColorArrayName = ['POINTS', 'T']
    contour2Display.LookupTable = tLUT
    contour2Display.SelectTCoordArray = 'None'
    contour2Display.SelectNormalArray = 'None'
    contour2Display.SelectTangentArray = 'None'
    contour2Display.OSPRayScaleArray = 'T'
    contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour2Display.SelectOrientationVectors = 'None'
    contour2Display.ScaleFactor = 0.09869538770761602
    contour2Display.SelectScaleArray = 'T'
    contour2Display.GlyphType = 'Arrow'
    contour2Display.GlyphTableIndexArray = 'T'
    contour2Display.GaussianRadius = 0.0049347693853808005
    contour2Display.SetScaleArray = ['POINTS', 'T']
    contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour2Display.OpacityArray = ['POINTS', 'T']
    contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour2Display.DataAxesGrid = 'GridAxesRepresentation'
    contour2Display.PolarAxes = 'PolarAxesRepresentation'
    contour2Display.SelectInputVectors = [None, '']
    contour2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    contour2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour2Display.ScaleTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour2Display.OpacityTransferFunction.Points = [700.0, 0.0, 0.5, 0.0, 2300.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(slice2, renderView1)

    # show color bar/color legend
    contour2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'PVD Reader'
    coupled_elsPerRad8pvd = PVDReader(registrationName='coupled_elsPerRad8.pvd', FileName=coupledDataSet)
    coupled_elsPerRad8pvd.CellArrays = ['ActiveElements', 'physicalDomain']
    coupled_elsPerRad8pvd.PointArrays = ['T', 'Pulse', 'ActiveNodes', 'gammaNodes', 'forcedDofs']


    # show data in view
    coupled_elsPerRad8pvdDisplay = Show(coupled_elsPerRad8pvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    coupled_elsPerRad8pvdDisplay.Representation = 'Surface'
    coupled_elsPerRad8pvdDisplay.ColorArrayName = [None, '']
    coupled_elsPerRad8pvdDisplay.SelectTCoordArray = 'None'
    coupled_elsPerRad8pvdDisplay.SelectNormalArray = 'None'
    coupled_elsPerRad8pvdDisplay.SelectTangentArray = 'None'
    coupled_elsPerRad8pvdDisplay.OSPRayScaleArray = 'ActiveNodes'
    coupled_elsPerRad8pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    coupled_elsPerRad8pvdDisplay.SelectOrientationVectors = 'None'
    coupled_elsPerRad8pvdDisplay.ScaleFactor = 0.36000000000000004
    coupled_elsPerRad8pvdDisplay.SelectScaleArray = 'ActiveNodes'
    coupled_elsPerRad8pvdDisplay.GlyphType = 'Arrow'
    coupled_elsPerRad8pvdDisplay.GlyphTableIndexArray = 'ActiveNodes'
    coupled_elsPerRad8pvdDisplay.GaussianRadius = 0.018000000000000002
    coupled_elsPerRad8pvdDisplay.SetScaleArray = ['POINTS', 'ActiveNodes']
    coupled_elsPerRad8pvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    coupled_elsPerRad8pvdDisplay.OpacityArray = ['POINTS', 'ActiveNodes']
    coupled_elsPerRad8pvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    coupled_elsPerRad8pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    coupled_elsPerRad8pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    coupled_elsPerRad8pvdDisplay.ScalarOpacityUnitDistance = 0.060946451134645946
    coupled_elsPerRad8pvdDisplay.OpacityArrayName = ['POINTS', 'ActiveNodes']
    coupled_elsPerRad8pvdDisplay.SelectInputVectors = [None, '']
    coupled_elsPerRad8pvdDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    coupled_elsPerRad8pvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Slice'
    slice3 = Slice(registrationName='Slice3', Input=coupled_elsPerRad8pvd)
    slice3.SliceType = 'Plane'
    slice3.HyperTreeGridSlicer = 'Plane'
    slice3.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice3.SliceType.Origin = [0.0, 0.0, 0.0]
    slice3.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.0]
    slice3.SliceType.Normal = [0.0, 0.0, 1.0]


    # Properties modified on slice3
    slice3.Triangulatetheslice = 0

    # Properties modified on slice3.SliceType
    slice3.SliceType.Origin = [0.0, 0.0, 0.0]
    slice3.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice3Display = Show(slice3, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice3Display.Representation = 'Surface'
    slice3Display.ColorArrayName = [None, '']
    slice3Display.SelectTCoordArray = 'None'
    slice3Display.SelectNormalArray = 'None'
    slice3Display.SelectTangentArray = 'None'
    slice3Display.OSPRayScaleArray = 'ActiveNodes'
    slice3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice3Display.SelectOrientationVectors = 'None'
    slice3Display.ScaleFactor = 0.36000000000000004
    slice3Display.SelectScaleArray = 'ActiveNodes'
    slice3Display.GlyphType = 'Arrow'
    slice3Display.GlyphTableIndexArray = 'ActiveNodes'
    slice3Display.GaussianRadius = 0.018000000000000002
    slice3Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    slice3Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice3Display.OpacityArray = ['POINTS', 'ActiveNodes']
    slice3Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice3Display.DataAxesGrid = 'GridAxesRepresentation'
    slice3Display.PolarAxes = 'PolarAxesRepresentation'
    slice3Display.SelectInputVectors = [None, '']
    slice3Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    slice3Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(coupled_elsPerRad8pvd, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice3.SliceType)

    # create a new 'Clip'
    clip1 = Clip(registrationName='Clip1', Input=slice3)
    clip1.ClipType = 'Plane'
    clip1.HyperTreeGridClipper = 'Plane'
    clip1.Scalars = ['POINTS', 'ActiveNodes']
    clip1.Value = 0.5

    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [-0.9999999999999993, 0.5, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip1.HyperTreeGridClipper.Origin = [-0.9999999999999993, 0.5, 0.0]

    # Properties modified on clip1
    clip1.Invert = 0

    # Properties modified on clip1.ClipType
    clip1.ClipType.Origin = [0.09999, 0.0, 0.0]
    clip1.ClipType.Normal = [-1.0, 0.0, 0.0]

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
    clip1Display.ScaleFactor = 0.28999899999999995
    clip1Display.SelectScaleArray = 'ActiveNodes'
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'ActiveNodes'
    clip1Display.GaussianRadius = 0.014499949999999998
    clip1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
    clip1Display.OpacityArray = ['POINTS', 'ActiveNodes']
    clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
    clip1Display.DataAxesGrid = 'GridAxesRepresentation'
    clip1Display.PolarAxes = 'PolarAxesRepresentation'
    clip1Display.ScalarOpacityUnitDistance = 0.1643614956606805
    clip1Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
    clip1Display.SelectInputVectors = [None, '']
    clip1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    clip1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(slice3, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip1.ClipType)

    # Properties modified on renderView1
    renderView1.CenterOfRotation = [-0.14999999478459358, 0.0494999997317791, 0.0]
    renderView1.CameraPosition = [-0.11521622491233216, 0.041634347699672225, 0.21346981293248826]
    renderView1.CameraFocalPoint = [-0.11521622491233216, 0.041634347699672225, 0.0]
    renderView1.CameraParallelScale = 0.2538725016795355

    # set active source
    SetActiveSource(contour1)

    # turn off scalar coloring
    ColorBy(contour1Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(tLUT, renderView1)

    # change solid color
    contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

    # Properties modified on contour1Display
    contour1Display.LineWidth = 4.0

    # Properties modified on contour1Display
    contour1Display.Opacity = 0.9

    # Properties modified on contour1Display
    contour1Display.LineWidth = 3.0

    # set active source
    SetActiveSource(contour2)

    # turn off scalar coloring
    ColorBy(contour2Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(tLUT, renderView1)

    # change solid color
    contour2Display.AmbientColor = [0.0, 0.6666666666666666, 1.0]
    contour2Display.DiffuseColor = [0.0, 0.6666666666666666, 1.0]

    # Properties modified on contour2Display
    contour2Display.LineWidth = 3.0

    # Properties modified on contour2Display
    contour2Display.Opacity = 0.9

    # set active source
    SetActiveSource(clip1)

    # create a new 'Contour'
    contour3 = Contour(registrationName='Contour3', Input=clip1)
    contour3.ContourBy = ['POINTS', 'T']
    contour3.Isosurfaces = [700.0, 1100.0, 1500.0, 1900.0, 2300.0, 2100.0]
    contour3.PointMergeMethod = 'Uniform Binning'

    # show data in view
    contour3Display = Show(contour3, renderView1, 'GeometryRepresentation')

    # get 2D transfer function for 'ActiveNodes'
    activeNodesTF2D = GetTransferFunction2D('ActiveNodes')

    # get color transfer function/color map for 'ActiveNodes'
    activeNodesLUT = GetColorTransferFunction('ActiveNodes')
    activeNodesLUT.TransferFunction2D = activeNodesTF2D
    activeNodesLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.998828272471228e-39, 0.0, 0.0, 1.0, 3.997656544942456e-39, 0.0, 1.0, 1.0, 5.878906683738906e-39, 0.0, 1.0, 0.0, 7.877734956210135e-39, 1.0, 1.0, 0.0, 9.876563228681361e-39, 1.0, 0.0, 0.0, 1.1757813367477812e-38, 0.878431372549, 0.0, 1.0]
    activeNodesLUT.ColorSpace = 'RGB'
    activeNodesLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    contour3Display.Representation = 'Surface'
    contour3Display.ColorArrayName = ['POINTS', 'ActiveNodes']
    contour3Display.LookupTable = activeNodesLUT
    contour3Display.SelectTCoordArray = 'None'
    contour3Display.SelectNormalArray = 'None'
    contour3Display.SelectTangentArray = 'None'
    contour3Display.OSPRayScaleArray = 'ActiveNodes'
    contour3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour3Display.SelectOrientationVectors = 'None'
    contour3Display.ScaleFactor = 0.07062400000018257
    contour3Display.SelectScaleArray = 'ActiveNodes'
    contour3Display.GlyphType = 'Arrow'
    contour3Display.GlyphTableIndexArray = 'ActiveNodes'
    contour3Display.GaussianRadius = 0.003531200000009128
    contour3Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour3Display.OpacityArray = ['POINTS', 'ActiveNodes']
    contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour3Display.DataAxesGrid = 'GridAxesRepresentation'
    contour3Display.PolarAxes = 'PolarAxesRepresentation'
    contour3Display.SelectInputVectors = [None, '']
    contour3Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    contour3Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(clip1, renderView1)

    # show color bar/color legend
    contour3Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get opacity transfer function/opacity map for 'ActiveNodes'
    activeNodesPWF = GetOpacityTransferFunction('ActiveNodes')
    activeNodesPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
    activeNodesPWF.ScalarRangeInitialized = 1

    # turn off scalar coloring
    ColorBy(contour3Display, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(activeNodesLUT, renderView1)

    # change solid color
    contour3Display.AmbientColor = [0.3333333333333333, 0.6666666666666666, 0.0]
    contour3Display.DiffuseColor = [0.3333333333333333, 0.6666666666666666, 0.0]

    # Properties modified on contour3Display
    contour3Display.LineWidth = 3.0

    # Properties modified on contour3Display
    contour3Display.Opacity = 0.9

    # set active source
    SetActiveSource(clip1)

    # show data in view
    clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

    # set active source
    SetActiveSource(clip1)

    # create a new 'Extract Edges'
    extractEdges1 = ExtractEdges(registrationName='ExtractEdges1', Input=clip1)

    # show data in view
    extractEdges1Display = Show(extractEdges1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    extractEdges1Display.Representation = 'Surface'
    extractEdges1Display.ColorArrayName = [None, '']
    extractEdges1Display.SelectTCoordArray = 'None'
    extractEdges1Display.SelectNormalArray = 'None'
    extractEdges1Display.SelectTangentArray = 'None'
    extractEdges1Display.OSPRayScaleArray = 'ActiveNodes'
    extractEdges1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractEdges1Display.SelectOrientationVectors = 'None'
    extractEdges1Display.ScaleFactor = 0.28999899551272396
    extractEdges1Display.SelectScaleArray = 'ActiveNodes'
    extractEdges1Display.GlyphType = 'Arrow'
    extractEdges1Display.GlyphTableIndexArray = 'ActiveNodes'
    extractEdges1Display.GaussianRadius = 0.014499949775636196
    extractEdges1Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    extractEdges1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractEdges1Display.OpacityArray = ['POINTS', 'ActiveNodes']
    extractEdges1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractEdges1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractEdges1Display.PolarAxes = 'PolarAxesRepresentation'
    extractEdges1Display.SelectInputVectors = [None, '']
    extractEdges1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    extractEdges1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(clip1, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(clip1)

    # show data in view
    clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

    # set active source
    SetActiveSource(slice1)

    # Properties modified on extractEdges1Display
    extractEdges1Display.Opacity = 0.3

    # create a query selection
    QuerySelect(QueryString='(T == max(T))', FieldType='POINT', InsideOut=0)

    # create a new 'Extract Selection'
    extractSelection1 = ExtractSelection(registrationName='ExtractSelection1', Input=slice1)

    # show data in view
    extractSelection1Display = Show(extractSelection1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    extractSelection1Display.Representation = 'Surface'
    extractSelection1Display.ColorArrayName = [None, '']
    extractSelection1Display.SelectTCoordArray = 'None'
    extractSelection1Display.SelectNormalArray = 'None'
    extractSelection1Display.SelectTangentArray = 'None'
    extractSelection1Display.OSPRayScaleArray = 'T'
    extractSelection1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSelection1Display.SelectOrientationVectors = 'None'
    extractSelection1Display.ScaleFactor = 0.1
    extractSelection1Display.SelectScaleArray = 'None'
    extractSelection1Display.GlyphType = 'Arrow'
    extractSelection1Display.GlyphTableIndexArray = 'None'
    extractSelection1Display.GaussianRadius = 0.005
    extractSelection1Display.SetScaleArray = ['POINTS', 'T']
    extractSelection1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSelection1Display.OpacityArray = ['POINTS', 'T']
    extractSelection1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSelection1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSelection1Display.PolarAxes = 'PolarAxesRepresentation'
    extractSelection1Display.ScalarOpacityUnitDistance = 0.0
    extractSelection1Display.OpacityArrayName = ['POINTS', 'T']
    extractSelection1Display.SelectInputVectors = [None, '']
    extractSelection1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    extractSelection1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSelection1Display.ScaleTransferFunction.Points = [2383.5577393095336, 0.0, 0.5, 0.0, 2384.057861328125, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSelection1Display.OpacityTransferFunction.Points = [2383.5577393095336, 0.0, 0.5, 0.0, 2384.057861328125, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(slice1, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Extract Selection'
    extractSelection2 = ExtractSelection(registrationName='ExtractSelection2', Input=slice1)

    # show data in view
    extractSelection2Display = Show(extractSelection2, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    extractSelection2Display.Representation = 'Surface'
    extractSelection2Display.ColorArrayName = [None, '']
    extractSelection2Display.SelectTCoordArray = 'None'
    extractSelection2Display.SelectNormalArray = 'None'
    extractSelection2Display.SelectTangentArray = 'None'
    extractSelection2Display.OSPRayScaleArray = 'T'
    extractSelection2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSelection2Display.SelectOrientationVectors = 'None'
    extractSelection2Display.ScaleFactor = 0.1
    extractSelection2Display.SelectScaleArray = 'None'
    extractSelection2Display.GlyphType = 'Arrow'
    extractSelection2Display.GlyphTableIndexArray = 'None'
    extractSelection2Display.GaussianRadius = 0.005
    extractSelection2Display.SetScaleArray = ['POINTS', 'T']
    extractSelection2Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSelection2Display.OpacityArray = ['POINTS', 'T']
    extractSelection2Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSelection2Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSelection2Display.PolarAxes = 'PolarAxesRepresentation'
    extractSelection2Display.ScalarOpacityUnitDistance = 0.0
    extractSelection2Display.OpacityArrayName = ['POINTS', 'T']
    extractSelection2Display.SelectInputVectors = [None, '']
    extractSelection2Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    extractSelection2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSelection2Display.ScaleTransferFunction.Points = [2383.5577393095336, 0.0, 0.5, 0.0, 2384.057861328125, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSelection2Display.OpacityTransferFunction.Points = [2383.5577393095336, 0.0, 0.5, 0.0, 2384.057861328125, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(slice1, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on extractSelection2
    extractSelection2.Input = slice2

    # set active source
    SetActiveSource(extractSelection2)

    # set active source
    SetActiveSource(slice3)

    # set active source
    SetActiveSource(slice3)

    # show data in view
    slice3Display = Show(slice3, renderView1, 'GeometryRepresentation')

    # create a query selection
    QuerySelect(QueryString='(T == max(T))', FieldType='POINT', InsideOut=0)

    # create a new 'Extract Selection'
    extractSelection3 = ExtractSelection(registrationName='ExtractSelection3', Input=slice3)

    # show data in view
    extractSelection3Display = Show(extractSelection3, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    extractSelection3Display.Representation = 'Surface'
    extractSelection3Display.ColorArrayName = [None, '']
    extractSelection3Display.SelectTCoordArray = 'None'
    extractSelection3Display.SelectNormalArray = 'None'
    extractSelection3Display.SelectTangentArray = 'None'
    extractSelection3Display.OSPRayScaleArray = 'ActiveNodes'
    extractSelection3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSelection3Display.SelectOrientationVectors = 'None'
    extractSelection3Display.ScaleFactor = 0.1
    extractSelection3Display.SelectScaleArray = 'ActiveNodes'
    extractSelection3Display.GlyphType = 'Arrow'
    extractSelection3Display.GlyphTableIndexArray = 'ActiveNodes'
    extractSelection3Display.GaussianRadius = 0.005
    extractSelection3Display.SetScaleArray = ['POINTS', 'ActiveNodes']
    extractSelection3Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSelection3Display.OpacityArray = ['POINTS', 'ActiveNodes']
    extractSelection3Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSelection3Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSelection3Display.PolarAxes = 'PolarAxesRepresentation'
    extractSelection3Display.ScalarOpacityUnitDistance = 0.0
    extractSelection3Display.OpacityArrayName = ['POINTS', 'ActiveNodes']
    extractSelection3Display.SelectInputVectors = [None, '']
    extractSelection3Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    extractSelection3Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.9208114632984963e-05, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSelection3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSelection3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(slice3, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(extractSelection2)

    # set active source
    SetActiveSource(extractSelection3)

    # change representation type
    extractSelection3Display.SetRepresentationType('Points')

    # Properties modified on extractSelection3Display
    extractSelection3Display.PointSize = 8.0

    # change solid color
    extractSelection3Display.AmbientColor = [0.054901960784313725, 0.47058823529411764, 0.6745098039215687]
    extractSelection3Display.DiffuseColor = [0.054901960784313725, 0.47058823529411764, 0.6745098039215687]

    # set active source
    SetActiveSource(extractSelection2)

    # change solid color
    extractSelection2Display.AmbientColor = [0.054901960784313725, 0.47058823529411764, 0.6745098039215687]
    extractSelection2Display.DiffuseColor = [0.054901960784313725, 0.47058823529411764, 0.6745098039215687]

    # Properties modified on extractSelection2Display
    extractSelection2Display.PointSize = 8.0

    # Properties modified on extractSelection2Display
    extractSelection2Display.Opacity = 0.9

    # set active source
    SetActiveSource(extractSelection3)

    # change solid color
    extractSelection3Display.AmbientColor = [0.25882352941176473, 0.47058823529411764, 0.050980392156862744]
    extractSelection3Display.DiffuseColor = [0.25882352941176473, 0.47058823529411764, 0.050980392156862744]

    # change solid color
    extractSelection3Display.AmbientColor = [0.0, 0.6666666666666666, 0.0]
    extractSelection3Display.DiffuseColor = [0.0, 0.6666666666666666, 0.0]

    # set active source
    SetActiveSource(extractSelection2)

    # change solid color
    extractSelection2Display.AmbientColor = [0.0, 0.3333333333333333, 1.0]
    extractSelection2Display.DiffuseColor = [0.0, 0.3333333333333333, 1.0]

    # set active source
    SetActiveSource(extractSelection1)

    # change representation type
    extractSelection1Display.SetRepresentationType('Points')

    # change solid color
    extractSelection1Display.DiffuseColor = [0.0, 0.0, 0.0]

    # Properties modified on extractSelection1Display
    extractSelection1Display.PointSize = 8.0

    # Properties modified on extractSelection1Display
    extractSelection1Display.Opacity = 0.9

    renderView1.OrientationAxesVisibility = 0
    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================

    setTime()
    # get layout
    layout1 = GetLayout()

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

    # save screenshot
    SaveScreenshot('/home/mslimani/acuario/moving_heat_source/run/3d_slob/figures/' + "allContours.png", renderView1, ImageResolution=[1586, 320],
        TransparentBackground=1)


    #--------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('coupled')
    parser.add_argument('ref')
    args = parser.parse_args()
    main(args.coupled, args.ref)
