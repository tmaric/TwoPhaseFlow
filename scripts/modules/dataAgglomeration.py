# dataAgglomeration.py

import pandas as pd
import numpy as np
import itertools
import re
import os
import sys

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import BoolProxy

import dataframeWithMetadata as dwm

class data_collector:
    """
    Collect data from text files in directories matching a user-prescribed pattern.

    This class uses a regular expression pattern to find all directories belonging to a study.
    From each directory a text file containing the data is read.
    The primary
    responsibility of this class is to create a dict of Pandas DataFrames where
    each dataframe represents data from a file. The variation ID (from PyFoam)
    is used as dict key. This is accompanied
    by some additional information which is required to create a
    Pandas multiindex for the collected data.

    Data collected by this class:
        - study_dataframes: dictionary of Pandas DataFrames containing the data read from files.
                                The variation ID is used as dictionary key.
        - valid_variations: list of variation numbers. Only those with valid data are considered.
        - invalid_variations: dictionary of all failed variations that maps the variation number
                                to the reason why it is considered failed.
        - datapoints_per_variant_: dictionary mapping the IDs of valid variations to their
                                    corresponding number of datapoints.
                                    Required for constructing the Pandas Multiindex.
    """

    def __init__(self, directory_pattern, path_to_file, subdirectory=""):
        """
        Instantiate a data_collector with a directory- and file pattern.

        This constructor requires two arguments:
            - directory_pattern: string containing a regular expression. All directories
                matching this pattern are considered for data collection.
            - path_to_file: path to the text file within a study directory containing the
                data to be collected, e.g. 'postProcessing/minMaxU/0/fieldMinMax.dat'.

            Keyword arguments:
            - subdirectory: name of the subdirectory where the studydirectories are located. (default: current directory)
        """
        self.subdirectory = subdirectory

        # Regular expressions patterns for finding directories
        # associated with a study and the files containing the data
        self.directory_pattern = re.compile(directory_pattern)
        self.path_to_file = path_to_file
        self.variation_number_pattern = re.compile("[0-9]{5}")

        # Collected data
        self.study_dataframes = dict()
        self.valid_variations = list()
        self.invalid_variations = dict()
        self.datapoints_per_variant_ = dict()

        # Avoid duplicated colection of data
        self.already_collected = False

    def set_directory_pattern(self, pattern):
        """Set the pattern of the directories which shall be considered for collection."""
        self.directory_pattern = re.compile(pattern)

    def study_directories(self):
        """Find all directories in the current folder matching the pattern and return their names as a list."""
        study_folders = []
        study_directory = os.path.join(os.getcwd(), self.subdirectory)

        for file in os.listdir(study_directory):
            belongs_to_study = self.directory_pattern.match(file)
            if belongs_to_study and os.path.isdir(os.path.join(self.subdirectory,file)):
                study_folders.append(os.path.join(study_directory, file))

        # Ensure that the folders are sorted so the data can correctly
        # be matched with the multiindex later (TT)
        study_folders.sort()

        return study_folders

    def extract_variation_number(self, directory):
        """Extract the PyFoam variation number from the directory name and return it as an integer."""
        # Only consider the last folder of the given path
        variation_folder = directory.split('/')[-1]
        return int(self.variation_number_pattern.search(variation_folder).group())

    def has_valid_results(self, directory):
        """
        Search the given directory for a valid data file.

        Checks if the given directory contains a data file as specified by
        the user. Returns False if either the file does not exist
        or it is empty.
        """
        # Check if a data file exists
        variation_number = self.extract_variation_number(directory)
        full_file_path = os.path.join(directory, self.path_to_file)

        if not os.path.isfile(full_file_path):
            self.invalid_variations[variation_number] = "No data file"
            return False

        if os.stat(full_file_path).st_size == 0:
            self.invalid_variations[variation_number] = "File empty"
            return False
        else:
            return True


    def read_dataframe_from_file(self, file_name):
        """
        Read Pandas DataFrame from text file, ignore columns with empty header.

        Distinguish two formats:
        - *.csv: written by function objects/applications of the argo project. Use ','
            as separator.
        - *.dat: written by OpenFOAM function object. Uses tab '\t' as separator.

        """
        # Remember: by default, pandas names a headerless column 'Unnamed N' 
        if file_name.endswith('.csv'):
            df = pd.read_csv(file_name, usecols=lambda x: not ('Unnamed' in x), comment='#')
        elif file_name.endswith('.dat'):
            df = pd.read_table(file_name, header=1)
            df.columns = [column.strip('#') for column in df.columns]
            df.columns = [column.strip(' ') for column in df.columns]
            df.columns = [column.rstrip(' ') for column in df.columns]
        else:
            sys.exit("Error: unknown file extension", file_name.rsplit('.')[0],
                     ". Valid extensions are .csv and .dat.")

        return df

        
    def collect_data(self):
        """Performs the actual data collection process."""
        if self.already_collected:
            return
        # Create study associated folder list
        directories = self.study_directories()

        # Iterate folders:
        for directory in directories:
            if self.has_valid_results(directory):
                full_file_path = os.path.join(directory, self.path_to_file)
                variation_dataframe = self.read_dataframe_from_file(full_file_path)
                variation_id = self.extract_variation_number(directory)
                self.valid_variations.append(variation_id)
                self.datapoints_per_variant_[variation_id] = len(variation_dataframe.index)
                self.study_dataframes[variation_id] = variation_dataframe

        self.already_collected = True

    #--------------------------------------------------------------------------
    # Interface member functions
    #--------------------------------------------------------------------------
    def agglomerated_study_data(self):
        """Returns collected data as a dictionary containing Pandas DataFrames as values."""
        self.collect_data()

        return self.study_dataframes

    def existing_variations(self):
        """Returns all valid variation numbers as a list of integers."""
        self.collect_data()

        return self.valid_variations

    def datapoints_per_variant(self):
        """Returns a dictionary mapping variation IDs to the number of datapoints."""
        self.collect_data()

        return self.datapoints_per_variant_

    def failed_variations(self):
        """
        Returns information about failed variations as a dictionary.
        
        The dictionary uses the variation number of the failed variations as keys
        and maps these to the reason why the variation is considered failed.
        NOTE: missing variation directories are not considered at all, meaning
        neither successful or failed.
        """
        self.collect_data()

        return self.invalid_variations


class multiindex_assembler:
    """
    Assembles a Pandas Multiindex for a PyFoam parameter study from a parameter file.

    Uses ParsedParameterFile from PyFoam to read the values of the study parameters
    from the given parameter file. The different parameter vectors are constructed
    using a Cartesian product. Parameter vectors of failed/missing variations
    are removed.
    Operates on a DataFrame before constructing the final mutliindex.
    """
    #--------------------------------------------------------------------------
    # Constructor
    #--------------------------------------------------------------------------
    def __init__(self, parameter_file_name):
        """Instantiate the multiindex_assembler given the name of a parameter file."""
        self.parameter_file_name = parameter_file_name
        self.parameter_vector_frame = None

    #--------------------------------------------------------------------------
    # Processing member functions
    #--------------------------------------------------------------------------
    def convert_nonhashable_types(self, value_list):
        """Convert non-hashable types to appropriate hashable ones"""
        for idx in range(len(value_list)):
            if type(value_list[idx]) == BoolProxy:
                value_list[idx] = str(value_list[idx])

    def compute_complete_index(self):
        """Compute the full set of parameter vectors and store them internally as a DataFrame."""
        variation_data=ParsedParameterFile(self.parameter_file_name,
                                          noHeader=True,
                                          noVectorOrTensor=True).getValueDict()
        parameter_names = []
        parameter_values = []
        n_variations = 1

        # Read the parameters and their corresponding set of values
        for parameter in variation_data["values"]:
            # Do not keep the solver. The study has not necessarily been executed
            #   with the one mentioned in the parameter file (TT)
            if parameter == "solver":
                continue

            values = variation_data["values"][parameter]
            # Keep all parameter values even if they are constant (TT)
            #if len(values) > 1:
            # Type check for values: PyFoam may read some of the parameter values,
            # e.g. logical switches (on/off), as non-hashable types which will
            # cause an error later when creating the dataframe.
            # These types are converted here
            self.convert_nonhashable_types(values)

            # It is necessary to prepend the values rather than to append
            # to obtain the variations in the same order as in PyFoam
            parameter_names.insert(0, parameter)
            parameter_values.insert(0, values)
            n_variations = n_variations*len(values)

        self.parameter_vector_frame = pd.DataFrame(columns=parameter_names,
                                                    index=range(0, n_variations-1))

        # The loop below does the magic of the Cartesian product
        index = 0
        for element in itertools.product(*parameter_values):
            self.parameter_vector_frame.loc[index] = list(element)
            index = index + 1

        # Add variation number as additional column to double check that
        # the mapping between parameter vector and data is correct later on
        self.parameter_vector_frame.loc[:,'variation_number'] = pd.Series(range(0,
                        n_variations), index=self.parameter_vector_frame.index)

    def remove_missing_variations(self, found_variations):
        """Drop variations which are not present in 'found_variations'."""
        self.parameter_vector_frame = self.parameter_vector_frame[
                    self.parameter_vector_frame['variation_number'].
                        isin(found_variations)]

    def has_constant_value(self, column):
        """Check if the given column of a dataframe contains a single, constant value."""
        ref_value = column[0]
        for value in column:
            if value != ref_value:
                return False

        return True

    def remove_constant_parameters(self):
        """
        Drop parameters which have a constant a.k.a. non-varying value.

        This function is relevant when only a variation subset of a study is
        considered. In this case, parameters which assume varying values
        according to the parameter file, might be constant (e.g. only variants
        of one resolution are simulated.) 
        """
        for column in self.parameter_vector_frame.columns:
            if self.has_constant_value(self.parameter_vector_frame.loc[:,column]):
                self.parameter_vector_frame = self.parameter_vector_frame.drop(column, axis=1)
        

    def assemble_multiindex(self, datapoints_per_variant):
        """
        Assemble the actual multiindex.
        
        Builds the multiindex from the dataframe by incorporating the number of data points
        per variant. This is required to setup up a unique mapping since a parameter vector
        maps to the entire data of a variant.
        """
        # Variation number is not required any more
        self.parameter_vector_frame = self.parameter_vector_frame.drop('variation_number', axis=1)

        # Ensure ascending order with respect to variation number
        self.parameter_vector_frame.sort_index()

        tuple_list = []

        for index,row in self.parameter_vector_frame.iterrows():
            parameter_vector = list(row)
            steps = range(datapoints_per_variant[index])

            for number in steps:
                step_index = list(parameter_vector)
                step_index.append(number)
                tuple_list.append(tuple(step_index))

        column_names = list(self.parameter_vector_frame.columns.values)
        column_names.append('step')

        return pd.MultiIndex.from_tuples(tuple_list, names=column_names)


    #--------------------------------------------------------------------------
    # Interface member functions
    #--------------------------------------------------------------------------
    def construct_multiindex(self, found_variations, datapoints_per_variant):
        """Construct and return a Pandas MultiIndex."""
        self.compute_complete_index()
        self.remove_missing_variations(found_variations)
        # Keep all parameter values even if they are constant (TT)
        #self.remove_constant_parameters()

        return self.assemble_multiindex(datapoints_per_variant)


class data_agglomerator:
    """
        Agglomerate the data of a parameter study in a multiindexed Pandas DataFrame.

        Intended to simplify further postprocessing of data and the generation
        of plots and tables from the data.
    """

    def __init__(self, parameter_file_name, path_to_file, directory_pattern=""):
        """
        Constructor of the data_agglomerator requiring the name of a parameter file.

        Keyword arguments:
        directory_pattern -- regular expression for finding the study directories (default: parameter file name)
        """
        self.parameter_file_name = parameter_file_name
        self.multiindex_assembler = multiindex_assembler(parameter_file_name)

        study_subdirectory, pattern = split_path_from_pattern(directory_pattern)
        if not directory_pattern:
            directory_pattern = parameter_file_name.replace('.', r"\.")
            directory_pattern = directory_pattern + r"_[0-9]{5}"
        else:
            directory_pattern = regex_from_study_directory_name(pattern)

        self.data_collector = data_collector(directory_pattern, path_to_file, subdirectory=study_subdirectory) 
        self.dataframe = pd.DataFrame()
        self.already_computed = False

    def assemble_dataframe(self, data, index):
        """Build the multiindexed dataframe from a dataframe and a multiindex."""
        self.dataframe = pd.concat(data)
        self.dataframe.index = index

    def compute_dataframe(self):
        """Toplevel function for dataframe assembly calling a set of lower-level member functions."""
        if self.already_computed:
            return

        # Get agglomerated data
        study_data = self.data_collector.agglomerated_study_data()
        existing_variations = self.data_collector.existing_variations()

        # Return an empty DataFrame for an empty parameter study. 
        if (len(existing_variations) == 0):
            raise ValueError("The number of existing variations equals 0, is the parameter study empty?")

        datapoints_per_variant = self.data_collector.datapoints_per_variant()

        # Create a mutliindex for the collected data
        study_index = self.multiindex_assembler.construct_multiindex(existing_variations, datapoints_per_variant)

        # Create multiindexed dataframe
        self.assemble_dataframe(study_data, study_index)

        self.already_computed = True

    def assemble_output_name(self, file_name, path):
        """
        Assemble the output name for the multiindexed dataframe File name suffix is .csv.
        
        Uses the name of the parameter file if an empty file name is passed.
        """
        if not file_name:
            file_name = self.parameter_file_name.split('.')[0]

        if 'csv' not in file_name.split('.'):
            file_name = file_name + '.csv'

        return os.path.join(path, file_name)

    #--------------------------------------------------------------------------
    # Interface member functions
    #--------------------------------------------------------------------------
    def study_dataframe(self):
        """Return the multiindexed dataframe of the study."""
        self.compute_dataframe()

        return self.dataframe

    def write_agglomerated_study_data(self, file_name="", path=""):
        """
        Write the multiindexed dataframe to a CSV file.

        Without arguments, the file is written to the current folder using the name of the
        parameter study.
        Keyword arguments:
        file_name -- file name for the dataframe (default "")
        path      -- path for the dataframe file (default current_working_directory)
        """
        self.compute_dataframe()
        path_and_name = self.assemble_output_name(file_name, path)
        dwm.write_dataframe_with_metadata(path_and_name, self.dataframe)

    def show_failed_variations(self):
        """Show a list of all variations present as directories but not containing valid data."""
        failed_variations = self.data_collector.failed_variations()
        print("Variants without valid data:")
        print("----------------------------")
        print("#Variation | Reason")
        print("----------------------------")

        for key, value in failed_variations.items():
            print(key, '\t', value)

def regex_from_study_directory_name(directory_name):
    """
    Generate a regular expression for the directory names of a study
    from a given example directory.
    """
    number_part = re.compile(r"_[0-9]{5}_")
    match = number_part.search(directory_name)
    directory_name = directory_name.replace(match.group(), r"_[0-9]{5}_")
    directory_name = directory_name.replace(".", r"\.")

    return r"^" + directory_name + r"$"

def split_path_from_pattern(directory_name):
    """
    Split the example directory name from the path to it.
    """
    if "/" not in directory_name:
        return ("", directory_name)
    else:
        return directory_name.rsplit("/", 1)
