# testReportCore.py

import pandas as pd
import numpy as np
import fnmatch
import os
import sys

#------------------------------------------------------------------------------
#   functions related to finding files and directories
#------------------------------------------------------------------------------
def list_of_parameter_files():
    """Find all files in current directory that end with '.parameter'"""

    parameterFiles = []
    workingDirectory = os.getcwd()

    for file in os.listdir(workingDirectory):
        if file.endswith(".parameter"):
            parameterFiles.append(file)

    return parameterFiles


def pattern_cases(pattern):
    """Finds all directories that fit a pattern"""

    patternCases = []
    workingDirectory = os.getcwd()

    for root, dirNames, fileNames in os.walk(workingDirectory):
      for dirName in fnmatch.filter(dirNames, pattern):
        patternCases.append(os.path.join(root, dirName))

    return patternCases


def parameter_study_directories(parameterStudyName):
    """Return a list containing all directory names beloging to
       the given study name.
    """

    directories = pattern_cases(parameterStudyName+"_0*")

    if not directories:
        print("Error: no directories for study ",parameterStudyName,"found. Exiting")
        sys.exit()

    return directories

def result_file_name(dirName):
    """If present, return the name of a result file in the given directory."""

    for file in os.listdir(dirName):
        if file.endswith("Results.csv"):
            return file

    return False


def result_path(envVariable, resultType):
    """Return a path for storing test results. If 'envVariable' is not set,
       the current folder is used.
    """

    path = ""
    typeSuffix = ""

    try:
        path = os.environ[envVariable]
    except(KeyError):
        path = os.path.abspath(os.curdir)

    if resultType == "table":
        typeSuffix = "tables"
    elif resultType == "figure":
        typeSuffix == "figures"

    path = os.path.join(path, typeSuffix)

    # Ensure folder exists
    if not os.path.exists(path):
        os.makedirs(path)

    return path


#------------------------------------------------------------------------------
#   functions related to extracting the parameter vector
#------------------------------------------------------------------------------
def contains_parameter(line):
    """Check if a string contains a meaningful parameter : value pair"""

    if line[0] == '/':
        return False
    elif line.isspace():
        return False

    # PyFoam also lists the the functions used in 'derivedParameters.py'
    # in the parameters file.
    # Remove them since only the computed values are of interest.
    if "<function" in line and "0x" in line:
        return False

    return True


def filter_parameter_vector(parameterVector, entriesToRemove):
    """Remove entries from a dictionary representing a parameter vector."""

    for key in entriesToRemove:
        del parameterVector[key]

    return parameterVector


def read_parameter_and_value(line):
    """Read a parameter : value pair from a string"""

    key,value = line.split(" ", 1)
    value_string = value.replace(";\n", "")

    # If the parameter has numeric values, convert them accordingly
    try:
        value = float(value_string)
    except ValueError:
        value = value_string

    return key,value


def read_parameter_vector(casePath):
    """Read the parameter vector of a given variation"""

    parameterFile = open(os.path.join(casePath,"PyFoamPrepareCaseParameters"))
    parameterVector = {}

    for line in parameterFile:
        if contains_parameter(line):
            parameter,value = read_parameter_and_value(line)
            parameterVector[parameter] = value

    # TODO: make this selectable, e.g. as a function argument (TT)
    columnsToRemove = ["casePath", "foamVersion", "foamFork", "caseName"]
    parameterVector = filter_parameter_vector(parameterVector, columnsToRemove)

    return parameterVector


#------------------------------------------------------------------------------
#   functions related to assembling the pandas dataframe
#------------------------------------------------------------------------------
def dataframe_from_csv(fileName):
    """Read dataframe from csv, ignore columns with empty header"""

    # Remember: by default, pandas names a headerless column 'Unnamed N' 
    df = pd.read_csv(fileName, usecols=lambda x: not ('Unnamed' in x))

    return df


def assemble_dataframe_header(casePath, additionalColumns=[]):
    """Assemble the header for the pandas dataframe. The header consists of
       of the parameter names found in 'PyFoamPrepareCaseParameters'
       and those given by 'additionalColumns'.
    """

    header = additionalColumns
    
    parameterVector = read_parameter_vector(casePath)

    for key in parameterVector:
        header.append(key)

    return header


def extract_variation_number_from_path(path):
    """Extract the variation number from the case directory name"""

    directoryName = path.rsplit('/',1)[1]

    variationNumber = ""

    for subString in directoryName.split("_"):
        if subString.isdigit():
            variationNumber = subString

    return variationNumber


def columns_with_varying_data(dataFrame):
    """Find columns in data frame which have varying data, e.g. determine 
       those parameters which are actually varied.
       Return a list of names of those columns.
    """

    columnNames = []

    for column in dataFrame.columns:
        refValue = dataFrame.loc[0,column]
        valuesVary = False

        for value in dataFrame.loc[:,column] :
            if value != refValue:
                valuesVary = True
                break
        
        if valuesVary:
            columnNames.append(column)

    return columnNames


def only_keep_given_columns(dataFrame, columnsToKeep):
    """Return a reduced dataframe which only contains the specified columns."""

    return dataFrame.loc[:,columnsToKeep]


def completed_runs(directoryList):
    """Return a list of all directories containing a result file."""

    completedRuns = []

    for directory in directoryList:
        if determine_run_status(directory) == "ok":
            completedRuns.append(directory)

    return completedRuns


def parameter_space(studyName):
    """Return a dataframe containing the parameter names as columns and the
       corresponding values.
    """

    # TODO: set this up from the parameter file rather then parsing the
    # directories of the study? Won't allow to filter out the failed
    # variations though.
    dirs = parameter_study_directories(studyName)
    dirs = completed_runs(dirs)
    dirs.sort()

    header = assemble_dataframe_header(dirs[0])

    df = pd.DataFrame(columns=header)

    for directory in dirs:

        # Extract parameters
        parameterVector = read_parameter_vector(directory)

        # Append data to container
        df = df.append(parameterVector, ignore_index=True)

    return df


def reduced_parameter_space(studyName):
    """Similar to 'parameter_space(...)', but removes constant parameters"""

    pSpace = parameter_space(studyName)
    varyingColumns = columns_with_varying_data(pSpace)
    reducedPSpace = only_keep_given_columns(pSpace, varyingColumns)

    return reducedPSpace


def multiindex_from_successful_variations(parameterSpace, iterationList):
    """Setup a multiindex from a dataframe 'parameterSpace'which contains the
       parameter names as columns and each row represents the parameter values
       of a successful variation. The 'iterationList' accounts for the number
       of runs in each variation.
    """

    tupleList = []

    for index,row in parameterSpace.iterrows():
        parameterVector = list(row)
        iterations = iterationList[index]

        for line in iterations:
            lineIndex = list(parameterVector)
            lineIndex.append(line)
            tupleList.append(tuple(lineIndex))
    
    columnNames = list(parameterSpace.columns.values)
    columnNames.append("iteration")

    index = pd.MultiIndex.from_tuples(tupleList, names=columnNames)

    return index


def add_variation_number(df, directoryName):
    """Add the variation number as a column to the dataframe"""

    varNum = extract_variation_number_from_path(directoryName)

    df['variation number'] = pd.Series(varNum, index=df.index)

    return df


def agglomerate_data(studyName, includeVarNum=False):
    """Agglomerates all the data frames from a parameter study into a single
       dataframe.
       Note: also returns a dataframe containing the reduced parameter space.
    """

    dirs = parameter_study_directories(studyName)
    dirs = completed_runs(dirs)
    dirs.sort()
     
    # Data frames picked up from CSV files stored in simulation directories.
    dfs = []
    pSpace = reduced_parameter_space(studyName)

    # Iteration number is required to address the single lines of a variations
    # results file
    iterationsLists = list()

    for directory in dirs: 
        csvFileName = os.path.join(os.curdir, directory, result_file_name(directory))
        csvDf = dataframe_from_csv(csvFileName)
        iterationsLists.append(list(csvDf.index))

        if includeVarNum:
            csvDf = add_variation_number(csvDf, directory)

        dfs.append(csvDf)
         
    # Store all data in a df with a multidimensional index
    mIndex = multiindex_from_successful_variations(pSpace, iterationsLists)
     
    agglDf = pd.concat(dfs)
    agglDf.index = mIndex
     
    return pSpace,agglDf


#------------------------------------------------------------------------------
#   functions related to crash reporting
#------------------------------------------------------------------------------
def determine_run_status(directory):
    """Check if a variation chrashed during the random / perturbed iterations. 
       This function exploits that the CSV file containing the results is
       only written after completing all iterations.
    """

    if result_file_name(directory):
        return "ok"
    else:
        return "crashed"


def crash_report_minimal_dataframe(dataFrame):
    """Only keep a minimal set of information in the dataframe for the
       crash report.
    """

    columns = columns_with_varying_data(dataFrame)

    # Ensure that the status column is included even if it has a constant value
    if "status" not in columns:
        columns.insert(1, "status")

    return only_keep_given_columns(dataFrame, columns)


def assemble_crash_report(parameterStudyName):
    """Assemble the crash report for the given parameter study."""

    # list of all directories belonging to the parameter study
    directoryList = parameter_study_directories(parameterStudyName)

    # Skip if there are no directories belonging to a certain parameter study,
    # e.g. it has not been run
    if not directoryList:
        return pd.DataFrame()

    # container holding information: variation number, completed/crashed, parameters
    df_header = assemble_dataframe_header(directoryList[0], ["variation number","status"])
    df = pd.DataFrame(columns=df_header)
    
    # for each directory
    for directory in directoryList:

        # Extract parameters
        parameterVector = read_parameter_vector(directory)

        # Determine variation number
        parameterVector["variation number"] = extract_variation_number_from_path(
                                                    directory)
        
        # Determine completed / crash
        parameterVector["status"] = determine_run_status(directory)
        
        # Append data to container
        df = df.append(parameterVector, ignore_index=True)

    # Compile report from container
    rf = crash_report_minimal_dataframe(df)
    fdf = rf[rf['status'] == 'crashed']

    return fdf

#------------------------------------------------------------------------------
#   functions related to hard condition violation reporting
#------------------------------------------------------------------------------
def columns_are_present(studyName, columnList):
    """Check if the csv files of a study contain the given columns."""

    studyDirectories = parameter_study_directories(studyName)

    for directory in studyDirectories:
        dataFileName = result_file_name(directory)
        if dataFileName:
            df = dataframe_from_csv(os.path.join(directory,dataFileName))
            metricNames = list(df.columns.values)

            for key in columnList:
                if key in metricNames:
                    return True

            return False


def assemble_hard_condition_violation_dataframe(parameterStudyName,
        conditionNames):
    """Assemble a dataframe containing relevant information for 
       reporting the violation of hard conditions.
    """

    reducedPSpace,df = agglomerate_data(parameterStudyName, True)

    header = ["variation number"]
    for colName in conditionNames:
        header.append(colName)

    for colName in list(reducedPSpace.columns):
        header.append(colName)

    reportDf = pd.DataFrame(columns=header)

    for index,row in df.iterrows():

        for condition in conditionNames:
            if row[condition] != 1.0:
                # This is ugly...
                newRow = dict()
                newRow[header[0]] = row[header[0]]
                for con in conditionNames:
                    newRow[con] = row[con]
                offset = len(conditionNames) + 1
                for I in range(len(index)-1):
                    newRow[header[I+offset]] = index[I]
                reportDf = reportDf.append(newRow,ignore_index=True)
                break
    
    reportDf.drop_duplicates(inplace=True)

    return reportDf


def compile_hard_condition_violation_tex_table(hardConditionColumns):
    """Write a tex table reporting the variations for which the hard
       conditions given in 'hardConditionColumns' are violated.
       The assumption is that that value of 0 indicates a violation of
       the condition. If the specified columns are not found in the
       results, no report is generated.
    """

    outputPath = result_path("TESTRESULTS", "table")

    # --- loop over all studies ---
    parameterFiles = list_of_parameter_files()

    for study in parameterFiles:
        # An application using this function may use it in a test case
        # not matching the application, e.g. when running the complete
        # test suite and executing each evaluation function in each test case.
        # In this case, simply skip the parameter studies.
        if not columns_are_present(study,hardConditionColumns):
            continue

        report = assemble_hard_condition_violation_dataframe(study,
                        hardConditionColumns)

        # TODO: add nice formatters?
        reportType = "-hard_condition_violation_report.tex"
        outputFileName = os.path.join(outputPath, study + reportType)
        report.to_latex(outputFileName, index=False)


#------------------------------------------------------------------------------
#   miscellaneous auxillary functions
#------------------------------------------------------------------------------
def convert_values_to_float_if_valid(dictionary):
    for key in dictionary:
        try:
            dictionary[key] = float(dictionary[key])
        except ValueError:
            dictionary[key] = dictionary[key]

    return dictionary

