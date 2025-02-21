# state file generated using paraview version 5.12.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Open FOAM Reader'
ratio100foam = OpenFOAMReader(registrationName='ratio100.foam', FileName='/work/groups/da_mma_b/jun/rhoVOF/TwoPhaseFlow/run/rhoVoF/risingBubbleAcc3D/N80_3D_noST_CylinderV9-timeStep-MatthiasSetup_00001_templateCase/cylinderv9.foam')
ratio100foam.MeshRegions = ['internalMesh']
ratio100foam.CellArrays = ['K_', 'alpha.water', 'p', 'p_rgh', 'reconstructedDistanceFunction', 'visRDF']

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=ratio100foam)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.Triangulatetheslice = 0
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [1e-06, 0.0, 0.0]

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=slice1)
calculator1.ResultArrayName = 'CellCenter'
calculator1.Function = 'coords'

# create a new 'Point Data to Cell Data'
pointDatatoCellData1 = PointDatatoCellData(registrationName='PointDatatoCellData1', Input=calculator1)
pointDatatoCellData1.ProcessAllArrays = 0
pointDatatoCellData1.PointDataArraytoprocess = ['CellCenter']
pointDatatoCellData1.PassPointData = 1

# create a new 'Slice'
slice4 = Slice(registrationName='Slice4', Input=pointDatatoCellData1)
slice4.SliceType = 'Plane'
slice4.HyperTreeGridSlicer = 'Plane'
slice4.SliceOffsetValues = [0.0]
slice4.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice4.SliceType.Origin = [0.0, 1e-07, 0.0]
slice4.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice4.HyperTreeGridSlicer.Origin = [9.999999974752427e-07, 0.0, 0.0]

# create a new 'Slice'
slice3 = Slice(registrationName='Slice3', Input=pointDatatoCellData1)
slice3.SliceType = 'Plane'
slice3.HyperTreeGridSlicer = 'Plane'
slice3.SliceOffsetValues = [0.0]
slice3.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice3.SliceType.Origin = [0.0, 0.00055501, 0.0]
slice3.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice3.HyperTreeGridSlicer.Origin = [9.999999974752427e-07, 0.0, 0.0]

# create a new 'Slice'
slice2 = Slice(registrationName='Slice2', Input=pointDatatoCellData1)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]
slice2.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [0.0, 0.0041137, 0.0]
slice2.SliceType.Normal = [0.0, 1.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [9.999999974752427e-07, 0.0, 0.0]

# ----------------------------------------------------------------
# restore active source
# SetActiveSource(slice4)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Save data
SaveData("farCylinder.csv", slice2, Precision=12, FieldAssociation='Cell Data')
SaveData("nearCylinder.csv", slice3, Precision=12, FieldAssociation='Cell Data')
SaveData("zAxis.csv", slice4, Precision=12, FieldAssociation='Cell Data')

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
## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
