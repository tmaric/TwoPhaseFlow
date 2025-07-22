# state file generated using paraview version 5.13.0
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

import argparse
import os

#### import the simple module from the paraview
from paraview.simple import *
import sys
#### disable automatic camera reset on 'Show'
#paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------
def find_or_create_foam_file(directory_path):
    """
    This function checks if a file with the `.foam` extension exists in the given directory.
    If found, it returns the file name. If not, it creates a file named `caseFoam.foam`.

    Args:
    directory_path (str): The relative or absolute path to the directory.

    Returns:
    str: The path of the found or created .foam file.
    """
    # Convert relative path to absolute path
    absolute_path = os.path.abspath(directory_path)

    # Check if the path exists and is a directory
    if not os.path.exists(absolute_path):
        raise FileNotFoundError(f"The directory at {absolute_path} does not exist.")
    if not os.path.isdir(absolute_path):
        raise NotADirectoryError(f"{absolute_path} is not a directory.")

    # Search for any .foam file in the directory
    for file_name in os.listdir(absolute_path):
        if file_name.endswith(".foam"):
            print(f"Found .foam file: {file_name}")
            return os.path.join(absolute_path, file_name)

    # If no .foam file is found, create `caseFoam.foam`
    new_foam_file = os.path.join(absolute_path, "caseFoam.foam")
    with open(new_foam_file, 'w') as file:
        file.write("This is a new .foam file.\n")

    print(f"No .foam file found. Created new file: caseFoam.foam")
    return new_foam_file

def extractPressure(case_path):
	foam_path = find_or_create_foam_file(case_path)
	# create a new 'Open FOAM Reader'
	Casefoam = OpenFOAMReader(registrationName='Case.foam', FileName=foam_path)
	Casefoam.MeshRegions = ['internalMesh']
	Casefoam.CellArrays = ['K_', 'U', 'alpha.water', 'interfaceCentre.water', 'interfaceNormal.water', 'p', 'p_rgh', 'reconstructedDistanceFunction', 'visRDF']

	# create a new 'Extract Cells By Region'
	extractCellsByRegion1 = ExtractCellsByRegion(registrationName='ExtractCellsByRegion1', Input=Casefoam)
	extractCellsByRegion1.IntersectWith = 'Box'

	# init the 'Box' selected for 'IntersectWith'
	extractCellsByRegion1.IntersectWith.Position = [-0.005, -0.005, -0.005]
	extractCellsByRegion1.IntersectWith.Length = [0.03, 0.01, 0.01]

	# create a new 'Calculator'
	calculator1 = Calculator(registrationName='Calculator1', Input=extractCellsByRegion1)
	calculator1.ResultArrayName = 'CellCenters'
	calculator1.Function = 'coords'

	# create a new 'Point Data to Cell Data'
	pointDatatoCellData1 = PointDatatoCellData(registrationName='PointDatatoCellData1', Input=calculator1)
	pointDatatoCellData1.ProcessAllArrays = 0
	pointDatatoCellData1.PointDataArraytoprocess = ['CellCenters']
	pointDatatoCellData1.PassPointData = 1

	# create a new 'Extract Cells Along Line'
	extractCellsAlongLine1 = ExtractCellsAlongLine(registrationName='ExtractCellsAlongLine1', Input=pointDatatoCellData1)
	extractCellsAlongLine1.Point1 = [0.02401, 1e-05, -0.005]
	extractCellsAlongLine1.Point2 = [0.02401, 1e-05, 0.005]

	# ----------------------------------------------------------------
	# restore active source
	#SetActiveSource(extractCellsAlongLine1)
	# ----------------------------------------------------------------
	absolute_path = os.path.abspath(case_path)
	# save data
	SaveData(absolute_path + '/Pressure_nearWall.csv', proxy=extractCellsAlongLine1, ChooseArraysToWrite=1,
	    CellDataArrays=['CellCenters', 'p', 'p_rgh'],
	    Precision=12,
	    FieldAssociation='Cell Data')
	#Disconnect()
	#sys.exit(0)
	os._exit(0)

if __name__ == "__main__":
    # Create an ArgumentParser object to handle the command-line argument
    parser = argparse.ArgumentParser(description="Find or create a .foam file in a given directory.")

    # Add an argument for the directory path
    parser.add_argument("directory_path", type=str, help="The relative or absolute path to the directory.")

    # Parse the command-line arguments
    args = parser.parse_args()

    extractPressure(args.directory_path)
    
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
