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
renderView1.ViewSize = [853, 543]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.0729999989271164, 0.000560699962079525, 0.00133500003721565]
renderView1.UseLight = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.10414296760236164, 0.000560699962079525, 0.00133500003721565]
renderView1.CameraFocalPoint = [0.0729999989271164, 0.000560699962079525, 0.00133500003721565]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraViewAngle = 17.35052754982415
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
layout1.SetSize(853, 543)

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

# create a new 'XML PolyData Reader'
sequencedVTK = XMLPolyDataReader(registrationName='sequencedVTK', FileName=['/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK'])
sequencedVTK.TimeArray = 'None'

# create a new 'XML PolyData Reader'
freeSurf0001vtp = XMLPolyDataReader(registrationName='freeSurf.0001.vtp*', FileName=['/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0001.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0002.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0003.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0004.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0005.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0006.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0007.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0008.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0009.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0010.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0011.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0012.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0013.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0014.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0015.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0016.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0017.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0018.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0019.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0020.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0021.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0022.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0023.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0024.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0025.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0026.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0027.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0028.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0029.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0030.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0031.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0032.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0033.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0034.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0035.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0036.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0037.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0038.vtp', '/work/scratch/lm74veso/reducedCasePLICaRDF800/common/sequencedVTK/freeSurf.0039.vtp'])
freeSurf0001vtp.CellArrayStatus = ['p', 'U']

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

tk = GetTimeKeeper()
tsteps = tk.TimestepValues
#print(tsteps)
i=0
p=178  #(luis) start frame needs to be matched manually
print(tk)
print(tsteps)
for t in tsteps:
    tk.Time=tsteps[i]   #(luis) changes the timestep
    rpm800_Layer2_Setpng = CreateTexture(filename='/home/lm74veso/ExperimentalData/ExperimentsLuis/BW800RPM/Rpm800_Layer2_Set40'+str(p)+'.png')

    # create a new 'Slice'
    slice4 = Slice(registrationName='Slice4', Input=freeSurf0001vtp)
    slice4.SliceType = 'Plane'
    slice4.HyperTreeGridSlicer = 'Plane'
    slice4.SliceOffsetValues = [0.0]
    slice4.PointMergeMethod = 'Uniform Binning'

    # init the 'Plane' selected for 'SliceType'
    slice4.SliceType.Origin = [0.07300000078976154, -0.00010568194556981325, 0.0003461308078840375]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice4.HyperTreeGridSlicer.Origin = [0.07299929484724998, -0.00010568194556981325, 0.0003461308078840375]

    # create a new 'Transform'
    transform3 = Transform(registrationName='Transform3', Input=slice4)
    transform3.Transform = 'RotateAroundOriginTransform'

    # init the 'RotateAroundOriginTransform' selected for 'Transform'
    transform3.Transform.Originofrotation = [0.0729999989271164, 0.0, 0.0]
    transform3.Transform.Rotate = [0.0, 0.0, 180.0]

    # create a new 'Transform'
    transform4 = Transform(registrationName='Transform4', Input=transform3)
    transform4.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform4.Transform.Translate = [0.0, 5e-05, 5e-05]

    # create a new 'Plane'
    plane1 = Plane(registrationName='Plane1')
    plane1.Origin = [0.073, 0.0073959, 0.005607]
    plane1.Point1 = [0.073, 0.0073959, -0.002937]
    plane1.Point2 = [0.073, -0.0062745, 0.005607]

    # create a new 'Texture Map to Plane'
    textureMaptoPlane1 = TextureMaptoPlane(registrationName='TextureMaptoPlane1', Input=plane1)

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from textureMaptoPlane1
    textureMaptoPlane1Display = Show(textureMaptoPlane1, renderView1, 'GeometryRepresentation')

    # a texture
    #rpm800_Layer2_Setpng = FindTextureOrCreate(registrationName='Rpm800_Layer2_Set..png', filename='/home/local/CSI/lm74veso/DIRP/ExperimentsLuis/BW800RPM/Rpm800_Layer2_Set40178.png')

    # trace defaults for the display properties.
    textureMaptoPlane1Display.Representation = 'Surface'
    textureMaptoPlane1Display.ColorArrayName = [None, '']
    textureMaptoPlane1Display.SelectNormalArray = 'Normals'
    textureMaptoPlane1Display.SelectTangentArray = 'None'
    textureMaptoPlane1Display.SelectTCoordArray = 'Texture Coordinates'
    textureMaptoPlane1Display.Texture = rpm800_Layer2_Setpng
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

    # show data from transform4
    transform4Display = Show(transform4, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    transform4Display.Representation = 'Surface'
    transform4Display.AmbientColor = [1.0, 0.3333333333333333, 0.0]
    transform4Display.ColorArrayName = [None, '']
    transform4Display.DiffuseColor = [1.0, 0.3333333333333333, 0.0]
    transform4Display.LineWidth = 3.0
    transform4Display.SelectNormalArray = 'None'
    transform4Display.SelectTangentArray = 'None'
    transform4Display.SelectTCoordArray = 'None'
    transform4Display.TextureTransform = 'Transform2'
    transform4Display.OSPRayScaleFunction = 'Piecewise Function'
    transform4Display.Assembly = ''
    transform4Display.SelectedBlockSelectors = ['']
    transform4Display.SelectOrientationVectors = 'None'
    transform4Display.ScaleFactor = 0.0002413730719126761
    transform4Display.SelectScaleArray = 'None'
    transform4Display.GlyphType = 'Arrow'
    transform4Display.GlyphTableIndexArray = 'None'
    transform4Display.GaussianRadius = 1.2068653595633806e-05
    transform4Display.SetScaleArray = [None, '']
    transform4Display.ScaleTransferFunction = 'Piecewise Function'
    transform4Display.OpacityArray = [None, '']
    transform4Display.OpacityTransferFunction = 'Piecewise Function'
    transform4Display.DataAxesGrid = 'Grid Axes Representation'
    transform4Display.PolarAxes = 'Polar Axes Representation'
    transform4Display.SelectInputVectors = [None, '']
    transform4Display.WriteLog = ''

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
    animationScene1.AnimationTime = 5.0499998906161636e-05
    animationScene1.StartTime = 5.0499998906161636e-05
    animationScene1.EndTime = 0.0019694999791681767
    animationScene1.PlayMode = 'Snap To TimeSteps'

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(freeSurf0001vtp)
    # ----------------------------------------------------------------
    print(i)
    SaveScreenshot("/home/lm74veso/Results/PICs_SYNCCEDtoLOCAL_/LowCellPLIC/PLICContour"+str(i)+".png")
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
