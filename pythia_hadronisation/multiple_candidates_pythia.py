"""
multiple_candidates.py

This code is used to chech which events have multiple candidates, and for those that do one of the candidates is selected at random, while the others are removed.
The year of interest, size of the data and polarity to be analysed must be specified using the required flags --year --size --polarity. There is a fourth flag --path, which is not required. This one is used to specify the directory where the input data is located, and where the output file should be written. By default it is set to be the current working directory.
It outputs the data for each into 3 seperate root files, one containing only D0 events, another containing only D0bar, and the third containing all events.
This code is heavily inspired on the work of Camille Jarvis-Stiggants and Michael England. The minor modifications to the original code have simply added flexibility to it by Marc Oriol Pérez

Author: Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 15th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import os
import numpy as np
import argparse
from collections import Counter
import csv

# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:

    --path      Used to specify the directory in which the output files should be written. It is not required,
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
    
    return parser.parse_args()

def read_from_file():
    '''
    Opens a .txt file and reads the values of the event information.
    
    Returns the data in the .txt file into an array.
    '''
    data = np.genfromtxt('/afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/pythia_hadronisation/pythia_hadronisation68.csv', delimiter = ',', skip_header = 1)
    
    return data

def remove_multiple_candidates(data):
    """
    
    Args:

    Returns:
        _type_: _description_
    """
    # Count occurrences of each unique event number in the first column of the simulation data.
    event_number = [row[0] for row in data]
    counts = Counter(event_number)
    
    # Filter rows where count of event number is 1.
    no_multiple_candidates = [row for row in data if counts[row[0]] == 1]
    # Converts the tuple into a numpy array.
    no_multiple_candidates = np.array(no_multiple_candidates)

    return no_multiple_candidates

def selection_criteria(raw_data):
    """
    This function imposes the pT and rapidity (y) selection criteria on the simulated events.
    The requirements are that: 0 < pT < 10 GeV/c,
                               0 < y < 6.
    
    The selection requirements are imposed by using mask arrays. The combined mask is applied to the original data array to extract only the rows that meet both conditions. 
    This is done using array indexing, where only the rows corresponding to True values in the mask are selected, effectively filtering out the rows that don't satisfy the conditions.
    
    Args:
        raw_data (array): Array containing the raw pythia simulated data.
        pT (array): A subset of raw_data including only the 3rd column, which is the data for pT.
        rapidity (array): A subset of raw_data including only the 4th column, which is the data for rapidity.
        mask_pT (boolean array): A boolean array the same length as pT where the corresponding elements have values of True or False depending on if they meet the pT selection criteria.
        mask_rapidity (boolean array): A boolean array the same length as rapidity where the corresponding elements have values of True or False depending on if they meet the rapidity selection criteria.
        combined_mask (boolean array): A boolean array with the value of True only if there are corresponding values of True in the mask_pT and mask_rapidity arrays.
        filtered_data (array): The new event data array where events are only appended to this array if they have met the selection criteria.
    """
    
    # Extract columns of interest (pT and rapidity).
    pT = raw_data[:, 2]
    rapidity = raw_data[:, 3]
    # Create boolean masks for filtering.
    mask_pT = (pT > 0) & (pT < 10)
    mask_rapidity = (rapidity > 0) & (rapidity < 6)

    # Combine masks.
    combined_mask = mask_pT & mask_rapidity

    # Apply the combined mask to your data.
    filtered_data = raw_data[combined_mask]
    
    print(f"The number of simulated events after selection criteria is: {len(filtered_data)}" )
    
    return filtered_data
    
def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def save_file2(D0, D0bar, both, path):
    """_summary_
    """
    filenames = [f'{path}/D0_clean_pythia_data.csv', f'{path}/D0bar_clean_pythia_data.csv', f'{path}/clean_pythia_data.csv']
    data_list = [D0, D0bar, both]
    # Writing to CSV files
    for i in range(len(filenames)):
        with open(filenames[i], 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data_list[i])

def split_meson(data):
    """
    It takes data from events and splits them in 3 datasets, depending on the origin meson
    of the decay (D0, D0bar, any).

    Args:
        data (_type_): _description_

    Returns:
        _type_: _description_
    """
    length = len(data)
    PID = data[:, 1]
    
    # Define a mask array of same length as the data with all ones.
    mask = np.ones(length)
    
    # Returns True for an element in the D0 and D0bar mask arrays if the PID number is equal to 421 for D0 and -421 for D0bar.
    mask_D0 = np.logical_and(mask, PID==421)
    mask_D0bar = np.logical_and(mask, PID==-421)
    
    # split the simulated particles into mesons and antimesons
    D0_data = data[mask_D0]
    D0bar_data = data[mask_D0bar]

    print('checkpoint: Data has been split into D0 and D0bar subsets')
        
    return D0_data, D0bar_data, data

# - - - - - - - MAIN BODY - - - - - - - #

args=parse_arguments()

np.random.seed(482022) # implemented 04/08/2022

# Import data from Pythia.
pythia_data = read_from_file()
print(f"The number of simulated events before all cuts is: {len(pythia_data)}")
if len(pythia_data) > 0:
    selected_data = selection_criteria(pythia_data) # Removes events that do not have pT and rapidity within the required range.
    cut_data = remove_multiple_candidates(selected_data) # select at random multiple candidates to remove.
    print(f"The number of simulated events after multiple candidates are cut is: {len(cut_data)}")
    D0_DATA, D0BAR_DATA, DATA = split_meson(cut_data) # split data based on meson
    print(f"The number of simulated D0 mesons are {len(D0_DATA)}")
    print(f"The number of simulated D0bar mesons are {len(D0BAR_DATA)}")
#     save_all(args.year, args.size) # output data
    save_file2(D0_DATA, D0BAR_DATA, cut_data, args.path)
