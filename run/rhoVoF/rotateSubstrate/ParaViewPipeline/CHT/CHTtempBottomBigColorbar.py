# state file generated using paraview version 5.13.2
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1750, 1140]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.07300000078976154, 0.0, 0.0]
renderView1.UseLight = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.07300000078976154, 0.0, -0.06556921995160327]
renderView1.CameraFocalPoint = [0.07300000078976154, 0.0, 0.0]
renderView1.CameraViewUp = [-1.0, 2.220446049250313e-16, 0.0]
renderView1.CameraViewAngle = 22.748538011695906
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.01697056289599111
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [640, 480]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'Grid Axes 3D Actor'
renderView2.CenterOfRotation = [0.0729999989271164, -5.00003807246685e-06, 0.006000000052154064]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [0.115270066467944, -5.00003807246685e-06, 0.006000000052154064]
renderView2.CameraFocalPoint = [0.0729999989271164, -5.00003807246685e-06, 0.006000000052154064]
renderView2.CameraViewUp = [0.0, 0.0, 1.0]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 0.010940298517336059
renderView2.LegendGrid = 'Legend Grid Actor'
renderView2.PolarGrid = 'Polar Grid Actor'
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1750, 1140)

# create new layout object 'Layout #2'
layout2 = CreateLayout(name='Layout #2')
layout2.AssignView(0, renderView2)
layout2.SetSize(640, 480)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Open FOAM Reader'
a12Dfoam = OpenFOAMReader(registrationName='12D.foam', FileName='/work/scratch/lm74veso/reducedCaseCHTtestRealVeloSR/common/CHT.foam')
a12Dfoam.CaseType = 'Decomposed Case'
a12Dfoam.MeshRegions = ['/fluid/internalMesh', '/solid/internalMesh']
a12Dfoam.CellArrays = ['K_', 'T', 'T.air', 'T.water', 'U', 'alpha.air', 'alpha.water', 'boilingHeatFlow', 'cellLevel', 'interfaceCentre.water', 'interfaceNormal.water', 'massSource_', 'p', 'p_rgh', 'psi0_', 'reconstructedDistanceFunction', 'rho']

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

tk = GetTimeKeeper()
tsteps = tk.TimestepValues
#print(tsteps)
i=0
print(tk)
print(tsteps)
for t in tsteps:
    tk.Time=tsteps[i]   #(luis) changes the timestep


    # create a new 'Plane'
    plane1 = Plane(registrationName='Plane1')
    plane1.Origin = [0.073, 0.0073959, 0.005607]
    plane1.Point1 = [0.073, 0.0073959, -0.002937]
    plane1.Point2 = [0.073, -0.0062745, 0.005607]

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=a12Dfoam)
    contour1.ContourBy = ['POINTS', 'alpha.water']
    contour1.Isosurfaces = [0.5]
    contour1.PointMergeMethod = 'Uniform Binning'

    # create a new 'Texture Map to Plane'
    textureMaptoPlane1 = TextureMaptoPlane(registrationName='TextureMaptoPlane1', Input=plane1)

    # create a new 'Slice'
    slice4 = Slice(registrationName='Slice4', Input=a12Dfoam)
    slice4.SliceType = 'Plane'
    slice4.HyperTreeGridSlicer = 'Plane'
    slice4.Triangulatetheslice = 0
    slice4.PointMergeMethod = 'Uniform Binning'

    # init the 'Plane' selected for 'SliceType'
    slice4.SliceType.Origin = [0.07300000078976154, 0.0, 0.0]
    slice4.SliceType.Normal = [0.0, 0.0, 1.0]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice4.HyperTreeGridSlicer.Origin = [0.07300000078976154, 0.0, 0.004200000083073974]

    # create a new 'Slice'
    slice3 = Slice(registrationName='Slice3', Input=contour1)
    slice3.SliceType = 'Plane'
    slice3.HyperTreeGridSlicer = 'Plane'
    slice3.SliceOffsetValues = [0.0]
    slice3.PointMergeMethod = 'Uniform Binning'

    # init the 'Plane' selected for 'SliceType'
    slice3.SliceType.Origin = [0.07300000078976154, -0.00010568194556981325, 0.0003461308078840375]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice3.HyperTreeGridSlicer.Origin = [0.07299929484724998, -0.00010568194556981325, 0.0003461308078840375]

    # create a new 'Transform'
    transform2 = Transform(registrationName='Transform2', Input=slice3)
    transform2.Transform = 'RotateAroundOriginTransform'

    # init the 'RotateAroundOriginTransform' selected for 'Transform'
    transform2.Transform.Originofrotation = [0.0729999989271164, 0.0, 0.0]
    transform2.Transform.Rotate = [0.0, 0.0, 180.0]

    # create a new 'Transform'
    transform3 = Transform(registrationName='Transform3', Input=transform2)
    transform3.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform3.Transform.Translate = [0.0, 5e-05, 5e-05]
    transform3.Transform.Rotate = [-1.0, 0.0, 0.0]

    # create a new 'Slice'
    slice2 = Slice(registrationName='Slice2', Input=a12Dfoam)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    slice2.PointMergeMethod = 'Uniform Binning'

    # init the 'Plane' selected for 'SliceType'
    slice2.SliceType.Origin = [0.07300000078976154, 0.0, 0.006000000052154064]
    slice2.SliceType.Normal = [0.0, 1.0, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice2.HyperTreeGridSlicer.Origin = [0.07300000078976154, 0.0, 0.006000000052154064]

    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=a12Dfoam)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.Triangulatetheslice = 0
    slice1.SliceOffsetValues = [0.0]
    slice1.PointMergeMethod = 'Uniform Binning'

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.07300000078976154, 0.0, 0.006000000052154064]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.07300000078976154, 0.0, 0.006000000052154064]

    # create a new 'Clip'
    clip1 = Clip(registrationName='Clip1', Input=slice1)
    clip1.ClipType = 'Plane'
    clip1.HyperTreeGridClipper = 'Plane'
    clip1.Scalars = ['POINTS', 'p']
    clip1.Value = 43.6705322265625

    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [0.0729999989271164, -0.0068352, 0.006000000052154064]
    clip1.ClipType.Normal = [0.0, -1.0, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip1.HyperTreeGridClipper.Origin = [0.0729999989271164, 0.0, 0.006000000052154064]

    # create a new 'Clip'
    clip2 = Clip(registrationName='Clip2', Input=clip1)
    clip2.ClipType = 'Plane'
    clip2.HyperTreeGridClipper = 'Plane'
    clip2.Scalars = ['POINTS', 'p']
    clip2.Value = 43.6705322265625

    # init the 'Plane' selected for 'ClipType'
    clip2.ClipType.Origin = [0.0729999989271164, 0.0068252, 0.006000000052154064]
    clip2.ClipType.Normal = [0.0, 1.0, 0.0]

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip2.HyperTreeGridClipper.Origin = [0.0729999989271164, 0.002582400105893612, 0.006000000052154064]

    # create a new 'Clip'
    clip3 = Clip(registrationName='Clip3', Input=clip2)
    clip3.ClipType = 'Plane'
    clip3.HyperTreeGridClipper = 'Plane'
    clip3.Scalars = ['POINTS', 'p']
    clip3.Value = 43.6705322265625

    # init the 'Plane' selected for 'ClipType'
    clip3.ClipType.Origin = [0.0729999989271164, -5.00003807246685e-06, 0.005607]
    clip3.ClipType.Normal = [0.0, 0.0, 1.0]

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip3.HyperTreeGridClipper.Origin = [0.0729999989271164, -5.00003807246685e-06, 0.006000000052154064]

    # create a new 'Transform'
    transform1 = Transform(registrationName='Transform1', Input=clip3)
    transform1.Transform = 'RotateAroundOriginTransform'

    # init the 'RotateAroundOriginTransform' selected for 'Transform'
    transform1.Transform.Originofrotation = [0.0729999989271164, 0.0, 0.0]
    transform1.Transform.Rotate = [0.0, 0.0, 180.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from slice4
    slice4Display = Show(slice4, renderView1, 'GeometryRepresentation')

    # get 2D transfer function for 'T'
    tTF2D = GetTransferFunction2D('T')
    tTF2D.ScalarRangeInitialized = 1
    tTF2D.Range = [296.7, 298.85, 0.0, 1.0]

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')
    tLUT.AutomaticRescaleRangeMode = 'Clamp and update every timestep'
    tLUT.RescaleOnVisibilityChange = 1
    tLUT.TransferFunction2D = tTF2D
    tLUT.RGBPoints = [298.1495361328125, 0.23137254902, 0.298039215686, 0.752941176471, 298.14976501464844, 0.865, 0.865, 0.865, 298.1499938964844, 0.705882352941, 0.0156862745098, 0.149019607843]
    tLUT.NumberOfTableValues = 265
    tLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    slice4Display.Representation = 'Surface With Edges'
    slice4Display.ColorArrayName = ['CELLS', 'T']
    slice4Display.LookupTable = tLUT
    slice4Display.SelectNormalArray = 'None'
    slice4Display.SelectTangentArray = 'None'
    slice4Display.SelectTCoordArray = 'None'
    slice4Display.TextureTransform = 'Transform2'
    slice4Display.OSPRayScaleArray = 'p'
    slice4Display.OSPRayScaleFunction = 'Piecewise Function'
    slice4Display.Assembly = 'Hierarchy'
    slice4Display.SelectedBlockSelectors = ['']
    slice4Display.SelectOrientationVectors = 'None'
    slice4Display.ScaleFactor = 0.002400000020861626
    slice4Display.SelectScaleArray = 'p'
    slice4Display.GlyphType = 'Arrow'
    slice4Display.GlyphTableIndexArray = 'p'
    slice4Display.GaussianRadius = 0.00012000000104308129
    slice4Display.SetScaleArray = ['POINTS', 'p']
    slice4Display.ScaleTransferFunction = 'Piecewise Function'
    slice4Display.OpacityArray = ['POINTS', 'p']
    slice4Display.OpacityTransferFunction = 'Piecewise Function'
    slice4Display.DataAxesGrid = 'Grid Axes Representation'
    slice4Display.PolarAxes = 'Polar Axes Representation'
    slice4Display.SelectInputVectors = [None, '']
    slice4Display.WriteLog = ''

    # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
    slice4Display.ScaleTransferFunction.Points = [100000.0, 0.0, 0.5, 0.0, 100016.0, 1.0, 0.5, 0.0]

    # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
    slice4Display.OpacityTransferFunction.Points = [100000.0, 0.0, 0.5, 0.0, 100016.0, 1.0, 0.5, 0.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for tLUT in view renderView1
    tLUTColorBar = GetScalarBar(tLUT, renderView1)
    tLUTColorBar.Title = 'T'
    tLUTColorBar.ComponentTitle = ''
    tLUTColorBar.TitleBold = 1
    tLUTColorBar.TitleFontSize = 40
    tLUTColorBar.LabelBold = 1
    tLUTColorBar.LabelFontSize = 40
    tLUTColorBar.ScalarBarLength = 0.5

    # set color bar visibility
    tLUTColorBar.Visibility = 1

    # show color legend
    slice4Display.SetScalarBarVisibility(renderView1, True)
    # rescale color and/or opacity maps used to include current data range
    slice4Display.RescaleTransferFunctionToDataRange(True)
    #slice4Display.RescaleTransferFunctionToVisibleRange(renderView1, True)
    # ----------------------------------------------------------------
    # setup color maps and opacity maps used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')
    tPWF.Points = [298.1495361328125, 0.0, 0.5, 0.0, 298.1499938964844, 1.0, 0.5, 0.0]
    tPWF.ScalarRangeInitialized = 1
    
    Render(renderView1)
    # ----------------------------------------------------------------
    # setup animation scene, tracks and keyframes
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # initialize the timekeeper

    # get time animation track
    timeAnimationCue1 = GetTimeTrack()

    # initialize the animation track

    # get animation scene
    animationScene1 = GetAnimationScene()

    # initialize the animation scene
    animationScene1.ViewModules = [renderView1, renderView2]
    animationScene1.Cues = timeAnimationCue1
    animationScene1.AnimationTime = 5.09021e-05
    animationScene1.StartTime = 5.09021e-05
    animationScene1.EndTime = 0.00196971
    animationScene1.PlayMode = 'Snap To TimeSteps'

    # initialize the animation scene

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(slice4)
    # ----------------------------------------------------------------
    print(i)
    SaveScreenshot("/home/lm74veso/Results/PICs_SYNCCEDtoLOCAL_/CHTlowCellROTbottom/CHTlowCellROTbottom"+str(i)+".png")
    i=i+1

##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
