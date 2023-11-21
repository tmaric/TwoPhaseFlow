# dataframeWithMetadata.py

import pandas as pd

def variation_vectors(multiindex, parameters=[], drop_step=True):
    """
    Return a list of dictionaries where each dictionary corresponds to a parameter variation.

    Expects a Pandas multiindex as argument.

    Keyword arguments:
    parameters -- only consider levels with a name given in this list. Other
                  levels are dropped (default: empty list)
    drop_step  -- boolean if the level named 'step' shall be dropped from index
                  (default True)
    """
    if drop_step and "step" in multiindex.names:
        reduced_index = multiindex.droplevel(level="step").unique()
    else:
        reduced_index = multiindex

    # Keep only those parameters specified in the parameter argument
    if parameters:
        for parameter in reduced_index.names:
            if parameter not in parameters:
                reduced_index = reduced_index.droplevel(level=parameter).unique()

    vectors = []
    parameter_names = reduced_index.names

    for entry in reduced_index:
        # We need to ensure that entry is of type 'tuple' or a list-like object
        # Otherwise, zip will not return a valid key-value tuple.
        if not isinstance(entry, tuple):
            tmp = []
            tmp.append(entry)
            entry = tmp
        vectors.append(dict(zip(parameter_names, entry)))

    return vectors

class metadata_parser:
    """
    Parser for the metadata stored along with a dataframe.

    The current, supported format is as follows:
        # key_name (separator) value (newline)
    where 'value' is not allowed to contain a newline character.
    """

    def __init__(self, file_name, metadata_token='#', separator=':'):
        """
        Initialize a metadata_parser object given a file name.

        Keyword arguments:
        metadata_token -- character to indicate a metadata line in the CSV file (default '#')
        separator -- character separating key and value of metadata (default ':')
        """
        self.file_name = file_name
        self.metadata_token = metadata_token
        self.separator = separator
        self.metadata = {}

    def read_raw_data(self):
        """Read raw data and return it as a list of strings."""
        raw_data = []
        with open(self.file_name) as file:
            for line in file:
                if line[0] == self.metadata_token:
                    raw_data.append(line.lstrip('#'))
                else:
                    break
        return raw_data

    def parse_metadata(self):
        """Top level function calling the different parsing steps."""
        if {}:
            return

        raw_data = self.read_raw_data()

        for line in raw_data:
            key_separator_value = line.partition(self.separator)
            key = key_separator_value[0].strip()
            value = key_separator_value[2].strip()
            self.metadata[key] = value

    #--------------------------------------------------------------------------
    # Interface member functions
    #--------------------------------------------------------------------------
    def dataframe_metadata(self):
        self.parse_metadata()
        return self.metadata


class metadated_dataframe_reader:
    """
    Read a Pandas DataFrame from a CSV file which stores metadata as comments using '#'.
    """

    def __init__(self, file_name, multiindexed=True, metadata_token='#'):
        """
        Initialize an metadated_dataframe_reader object given a file name of a dataframe.
        
        Keyword arguments:
        multiindexed -- whether the dataframe stores a multiindex or not (default True)
        metadata_token -- character to indicate metadata/comment lines (default '#')
        """
        self.file_name = file_name
        self.has_multiindex = multiindexed
        self.metadata_token = metadata_token
        self.metadata_parser = metadata_parser(file_name, metadata_token)

    def read_dataframe(self):
        """Read and return the dataframe."""
        if self.has_multiindex:
            metadata = self.metadata_parser.dataframe_metadata()
            n_index_columns = int(metadata['n_index_columns'])
            return pd.read_csv(self.file_name, comment=self.metadata_token, index_col=list(range(n_index_columns)))
        else:
            return pd.read_csv(self.file_name, comment=self.metadata_token)
    
    def read_metadata(self):
        return self.metadata_parser.dataframe_metadata()

def write_dataframe_with_metadata(file_name, dataframe, metadata={}, has_multiindex=True):
    """
    Write a dataframe together with its metadata as comments.

    Metadata must be given as a dictionary. If the 'has_multiindex' flag is set,
    the number of index columns is added to the metadata.
    Keyword arguments:
    metadata -- dictionary containing the metadata (default: empty dict)
    has_multiindex -- flag if the dataframe uses a multiindex (default: True)
    """
    # Move index levels with constant value to metadata
    print(dataframe.index.levels)
    print(dataframe.index.names)
    levels_to_drop = []

    for idx in range(len(dataframe.index.names)):
        if len(dataframe.index.levels[idx]) == 1:
            levels_to_drop.append(idx)
            metadata[dataframe.index.names[idx]] = dataframe.index.levels[idx][0]

    dataframe.reset_index(level=levels_to_drop, drop=True, inplace=True)


    # If dataframe has a multiindex, write the number of index columns as
    # metadata
    if has_multiindex:
        metadata['n_index_columns'] = str(dataframe.index.nlevels)

    # Ensure proper file name suffix
    if not (file_name.endswith('.csv') or file_name.endswith('.CSV')):
        file_name = file_name + '.csv'

    # Write metadata
    dataframe_file = open(file_name, 'w')
    metadata_string = ""
    for key, value in metadata.items():
        metadata_string += "# " + str(key) + " : " + str(value) + "\n"
    dataframe_file.write(metadata_string)
    dataframe_file.close()

    # Append data
    dataframe.to_csv(file_name, mode='a')
