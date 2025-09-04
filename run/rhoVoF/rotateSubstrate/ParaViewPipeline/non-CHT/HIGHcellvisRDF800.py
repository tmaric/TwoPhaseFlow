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
renderView1.CenterOfRotation = [0.0729999989271164, 0.000560699962079525, 0.00133500003721565]
renderView1.UseLight = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.10414296760236164, 0.000560699962079525, 0.00133500003721565]
renderView1.CameraFocalPoint = [0.0729999989271164, 0.000560699962079525, 0.00133500003721565]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraViewAngle = 17.771428571428572
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.00806039341418469
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

# create a new 'Plane'
plane1 = Plane(registrationName='Plane1')
plane1.Origin = [0.073, 0.0073959, 0.005607]
plane1.Point1 = [0.073, 0.0073959, -0.002937]
plane1.Point2 = [0.073, -0.0062745, 0.005607]

# create a new 'Texture Map to Plane'
textureMaptoPlane1 = TextureMaptoPlane(registrationName='TextureMaptoPlane1', Input=plane1)

# create a new 'Open FOAM Reader'
a12Dfoam = OpenFOAMReader(registrationName='12D.foam', FileName='/work/scratch/lm74veso/HIGHreducedCaseFullPLICaRDF800/common/H800.foam')
a12Dfoam.CaseType = 'Decomposed Case'
a12Dfoam.MeshRegions = ['internalMesh']
a12Dfoam.CellArrays = ['K_', 'U', 'alpha.water', 'cellLevel', 'interfaceCentre.water', 'interfaceNormal.water', 'p', 'p_rgh', 'rAU', 'reconstructedDistanceFunction', 'visRDF']

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

tk = GetTimeKeeper()
tsteps = tk.TimestepValues
#print(tsteps)
i=0
p=180  #(luis) start frame needs to be matched manually
print(tk)
print(tsteps)
for t in tsteps:
    tk.Time=tsteps[i]   #(luis) changes the timestep
    rpm800_Layer2_Set40180 = CreateTexture(filename='/home/lm74veso/ExperimentalData/ExperimentsLuis/BW800RPM/Rpm800_Layer2_Set40'+str(p)+'.png')

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=a12Dfoam)
    contour1.ContourBy = ['POINTS', 'visRDF']
    contour1.Isosurfaces = [0.0]
    contour1.PointMergeMethod = 'Uniform Binning'

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

    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=a12Dfoam)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
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

    # create a new 'Transform'
    transform3 = Transform(registrationName='Transform3', Input=transform2)
    transform3.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform3.Transform.Translate = [0.0, 5e-05, 5e-05]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from textureMaptoPlane1
    textureMaptoPlane1Display = Show(textureMaptoPlane1, renderView1, 'GeometryRepresentation')

    # a texture
    #rpm800_Layer2_Set40180 = FindTextureOrCreate(registrationName='Rpm800_Layer2_Set40180', filename='/home/luigi/OpenFOAM/ExperimentalData/BW800RPM/Rpm800_Layer2_Set40180.png')

    # trace defaults for the display properties.
    textureMaptoPlane1Display.Representation = 'Surface'
    textureMaptoPlane1Display.ColorArrayName = [None, '']
    textureMaptoPlane1Display.SelectNormalArray = 'Normals'
    textureMaptoPlane1Display.SelectTangentArray = 'None'
    textureMaptoPlane1Display.SelectTCoordArray = 'Texture Coordinates'
    textureMaptoPlane1Display.Texture = rpm800_Layer2_Set40180
    textureMaptoPlane1Display.TextureTransform = 'Transform2'
    textureMaptoPlane1Display.OSPRayScaleArray = 'Normals'
    textureMaptoPlane1Display.OSPRayScaleFunction = 'Piecewise Function'
    textureMaptoPlane1Display.Assembly = ''
    textureMaptoPlane1Display.SelectedBlockSelectors = ['']
    textureMaptoPlane1Display.SelectOrientationVectors = 'None'
    textureMaptoPlane1Display.ScaleFactor = 0.001367039978504181
    textureMaptoPlane1Display.SelectScaleArray = 'None'
    textureMaptoPlane1Display.GlyphType = 'Arrow'
    textureMaptoPlane1Display.GlyphTableIndexArray = 'None'
    textureMaptoPlane1Display.GaussianRadius = 6.835199892520905e-05
    textureMaptoPlane1Display.SetScaleArray = ['POINTS', 'Normals']
    textureMaptoPlane1Display.ScaleTransferFunction = 'Piecewise Function'
    textureMaptoPlane1Display.OpacityArray = ['POINTS', 'Normals']
    textureMaptoPlane1Display.OpacityTransferFunction = 'Piecewise Function'
    textureMaptoPlane1Display.DataAxesGrid = 'Grid Axes Representation'
    textureMaptoPlane1Display.PolarAxes = 'Polar Axes Representation'
    textureMaptoPlane1Display.SelectInputVectors = ['POINTS', 'Normals']
    textureMaptoPlane1Display.WriteLog = ''

    # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
    textureMaptoPlane1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

    # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
    textureMaptoPlane1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

    # show data from transform3
    transform3Display = Show(transform3, renderView1, 'GeometryRepresentation')

    # get 2D transfer function for 'visRDF'
    visRDFTF2D = GetTransferFunction2D('visRDF')

    # get color transfer function/color map for 'visRDF'
    visRDFLUT = GetColorTransferFunction('visRDF')
    visRDFLUT.TransferFunction2D = visRDFTF2D
    visRDFLUT.RGBPoints = [-100.0, 0.0, 0.0, 0.0, -20.0, 0.901960784314, 0.0, 0.0, 60.0, 0.901960784314, 0.901960784314, 0.0, 100.0, 1.0, 1.0, 1.0]
    visRDFLUT.ColorSpace = 'RGB'
    visRDFLUT.NanColor = [0.0, 0.498039215686, 1.0]
    visRDFLUT.NumberOfTableValues = 1
    visRDFLUT.ScalarRangeInitialized = 1.0

    # trace defaults for the display properties.
    transform3Display.Representation = 'Surface'
    transform3Display.ColorArrayName = ['POINTS', 'visRDF']
    transform3Display.LookupTable = visRDFLUT
    transform3Display.LineWidth = 3.0
    transform3Display.SelectNormalArray = 'Normals'
    transform3Display.SelectTangentArray = 'None'
    transform3Display.SelectTCoordArray = 'None'
    transform3Display.TextureTransform = 'Transform2'
    transform3Display.OSPRayScaleArray = 'visRDF'
    transform3Display.OSPRayScaleFunction = 'Piecewise Function'
    transform3Display.Assembly = 'Hierarchy'
    transform3Display.SelectedBlockSelectors = ['']
    transform3Display.SelectOrientationVectors = 'U'
    transform3Display.ScaleFactor = 0.00023667984642088415
    transform3Display.SelectScaleArray = 'visRDF'
    transform3Display.GlyphType = 'Arrow'
    transform3Display.GlyphTableIndexArray = 'visRDF'
    transform3Display.GaussianRadius = 1.1833992321044207e-05
    transform3Display.SetScaleArray = ['POINTS', 'visRDF']
    transform3Display.ScaleTransferFunction = 'Piecewise Function'
    transform3Display.OpacityArray = ['POINTS', 'visRDF']
    transform3Display.OpacityTransferFunction = 'Piecewise Function'
    transform3Display.DataAxesGrid = 'Grid Axes Representation'
    transform3Display.PolarAxes = 'Polar Axes Representation'
    transform3Display.SelectInputVectors = ['POINTS', 'U']
    transform3Display.WriteLog = ''

    # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
    transform3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
    transform3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

    # ----------------------------------------------------------------
    # setup color maps and opacity maps used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # get opacity transfer function/opacity map for 'visRDF'
    visRDFPWF = GetOpacityTransferFunction('visRDF')
    visRDFPWF.Points = [-100.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]
    visRDFPWF.ScalarRangeInitialized = 1

    # ----------------------------------------------------------------
    # setup animation scene, tracks and keyframes
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # get time animation track
    timeAnimationCue1 = GetTimeTrack()

    # initialize the animation scene

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # initialize the timekeeper

    # initialize the animation track

    # get animation scene
    animationScene1 = GetAnimationScene()

    # initialize the animation scene
    animationScene1.ViewModules = [renderView1, renderView2]
    animationScene1.Cues = timeAnimationCue1
    animationScene1.AnimationTime = 5.04805e-05
    animationScene1.StartTime = 5.04805e-05
    animationScene1.EndTime = 0.0019696
    animationScene1.PlayMode = 'Snap To TimeSteps'

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(plane1)
    # ----------------------------------------------------------------
    print(i)
    SaveScreenshot("/home/lm74veso/Results/PICs_SYNCCEDtoLOCAL_/HighCellcountvisRDF/HighCellvisRDFContour"+str(i)+".png")
    i=i+1
    p=p+5 

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
