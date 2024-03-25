"""
multiple_candidates_pythia.py

This code reads in the outputted .csv file from runpythia.cpp. The .csv file contains the simulated event information, specifically the event number, PID number and the associated values of the transverse momentum (pT) and rapidity(y).
The script applys selection criteria so that only simulated events that have 0 < pT < 10 GeV/c and 0 < y < 6 are accepted. In addition events that are multiple candidates (have the same event number) are removed.

Author: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 15th February 2024
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import os
import numpy as np
import argparse
from collections import Counter
import csv
import pandas as pd

# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:

    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    Returns the parsed arguments.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--input",
        type=file_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found"
    )
    parser.add_argument(
        "--file_suffix",
        type=int,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found"
    )
    return parser.parse_args()

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string) 
    
def chunk_preprocessing(chunk):
    if len(chunk) > 0:
        selected_data = selection_criteria(chunk) # Removes events that do not have pT and rapidity within the required range.
        cut_data = remove_multiple_candidates(selected_data) # select at random multiple candidates to remove.
        print(f"The number of simulated events after multiple candidates are cut is: {len(cut_data)}")
        D0_DATA, D0BAR_DATA, DATA = split_meson(cut_data) # split data based on meson
        print(f"The number of simulated D0 mesons are {len(D0_DATA)}")
        print(f"The number of simulated D0bar mesons are {len(D0BAR_DATA)}")
        print("")
    
    return cut_data
    
def read_from_file():
    """
    Opens a .csv file and reads the values of the event information.
    Returns the data in the .csv file into an array.

    Returns:
        data (array): An array containging the simulated events information.
    """
    df_chunk = pd.read_csv(f'{args.input}', chunksize=1000000, usecols=columns_to_read_save)

    chunk_list = []  # append each chunk df here 

    # Each chunk is in df format
    for chunk in df_chunk:  
        # perform data filtering 
        chunk_filter = chunk_preprocessing(chunk)
    
        # Once the data filtering is done, append the chunk to list
        chunk_list.append(chunk_filter)
    
    # concat the list into dataframe 
    df_concat = pd.concat(chunk_list)

    return df_concat

def remove_multiple_candidates(data):
    """
    This function removes events that contain multiple candidates, which is determined when candidates have the same event number.
    Unlike the multiple candidates script for the real data, all multiple candidate events are relinquished. This is because we are not limited by the amount of data since we can just increase the simulated sample size.
    
    Args:
        data (DataFrame): DataFrame containing the pythia simulated data that has passed the selection criteria.

    Returns:
        no_multiple_candidates (DataFrame): DataFrame containing the simulated events that do not have multiple candidates.
    """
    # Count occurrences of each unique event number
    counts = data['Event '].value_counts()

    # Filter rows where count of event number is 1
    no_multiple_candidates = data[data['Event '].map(counts) == 1]

    return no_multiple_candidates


def selection_criteria(raw_data):
    """
    This function imposes the pT and rapidity (y) selection criteria on the simulated events.
    The requirements are that: 1.25 < pT < 10 GeV/c,
                               2 < y < 5.
                               
    Args:
        raw_data (DataFrame): DataFrame containing the raw pythia simulated data.

    Returns:
        filtered_data (DataFrame): The new DataFrame containing only the events that meet the selection criteria.
    """
    # Apply the selection criteria
    filtered_data = raw_data[(raw_data[' PT '] > 1.25) & (raw_data[' PT '] < 10) & (raw_data[' Y '] > 0) & (raw_data[' Y '] < 6)]
    
    print(f"The number of simulated events after selection criteria is: {len(filtered_data)}" )
    
    return filtered_data
    
def file_path(string):
    '''
    Checks if a given string is the path to a file.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(f"No such file: '{string}'")

def save_file(D0, D0bar, both, path, file_suffix):
    """
    Saves the selected simulated data into three .csv files: one each for D0 and D0bar mesons, and one for both mesons.

    Args:
        D0 (DataFrame): The simulated D0 mesons that meet all of the criteria.
        D0bar (DataFrame): The simulated D0bar mesons that meet all of the criteria.
        both (DataFrame): DataFrame containing the pythia simulated data that has passed the selection and multiple candidate requirements.
        path (string): The file path parsed in the path parser.
        file_suffix (string): Suffix to be appended to the file names.
    """
    filenames = [f'{path}/D0_clean_pythia_data_{file_suffix}.csv', f'{path}/D0bar_clean_pythia_data_{file_suffix}.csv', f'{path}/clean_pythia_data_{file_suffix}.csv']
    data_list = [D0, D0bar, both]
    
    # Writing to CSV files
    for filename, data in zip(filenames, data_list):
        data[columns_to_read_save].to_csv(filename, index=False)

def split_meson(data):
    """
    It takes data from events and splits them into three datasets, depending on the origin meson
    of the decay (D0, D0bar, any).

    Args:
        data (DataFrame): DataFrame containing the pythia simulated data that has passed the selection and multiple candidate requirements.

    Returns:
        D0_data (DataFrame): The simulated D0 mesons that meet all of the criteria.
        D0bar_data (DataFrame): The simulated D0bar mesons that meet all of the criteria.
        data (DataFrame): DataFrame containing the pythia simulated data that has passed the selection and multiple candidate requirements.
    """
    # Define masks for D0 and D0bar mesons
    mask_D0 = data[' PID '] == 421
    mask_D0bar = data[' PID '] == -421

    # Filter data based on masks
    D0_data = data[mask_D0]
    D0bar_data = data[mask_D0bar]

    print('checkpoint: Data has been split into D0 and D0bar subsets')

    # Return split data
    return D0_data, D0bar_data, data


# - - - - - - - MAIN BODY - - - - - - - #

args=parse_arguments()
columns_to_read_save = ['Event ', ' PID ', ' PT ', ' Y ']

# Import data from Pythia.
pythia_data = read_from_file()
num_simulated_events = len(pythia_data)
formatted_num_events = "{:,}".format(num_simulated_events)
print(f"The number of simulated events in the Pythia data is {formatted_num_events}.")

if len(pythia_data) > 0:
    D0_DATA, D0BAR_DATA, DATA = split_meson(pythia_data) # split data based on meson
    print(f"The number of simulated D0 mesons are {len(D0_DATA)}")
    print(f"The number of simulated D0bar mesons are {len(D0BAR_DATA)}")
    #saveall(args.year, args.size) # output data
    save_file(D0_DATA, D0BAR_DATA, DATA, args.path, args.file_suffix)
