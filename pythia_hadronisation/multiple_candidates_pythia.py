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
    """
    Opens a .csv file and reads the values of the event information.
    Returns the data in the .csv file into an array.

    Returns:
        data (array): An array containging the simulated events information.
    """
    data = np.genfromtxt('/afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/pythia_hadronisation/pythia_hadronisation68.csv', delimiter = ',', skip_header = 1)
    
    return data

def remove_multiple_candidates(data):
    """
    This function removes events that contain multiple candidates, which is determined when candidates have the same event number.
    Unlike the multiple candidates script for the real data, all multiple candidate events are relinquished. This is because we are not limited by the amount of data since we can just increase the simulated sample size.
    
    Args:
        data (array): Array containing the pythia simulated data that has passed the selection criteria.

    Returns:
        no_multiple_candidates (array): Array containing the simulated events that do not have multiple candidates.
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
    
    pT (array): A subset of raw_data including only the 3rd column, which is the data for pT.
    rapidity (array): A subset of raw_data including only the 4th column, which is the data for rapidity.
    mask_pT (boolean array): A boolean array the same length as pT where the corresponding elements have values of True or False depending on if they meet the pT selection criteria.
    mask_rapidity (boolean array): A boolean array the same length as rapidity where the corresponding elements have values of True or False depending on if they meet the rapidity selection criteria.
    combined_mask (boolean array): A boolean array with the value of True only if there are corresponding values of True in the mask_pT and mask_rapidity arrays.
    
    Args:
        raw_data (array): Array containing the raw pythia simulated data.

    Returns:
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

def save_file(D0, D0bar, both, path):
    """
    Saves the selected simulated data into three .csv files: one each for D0 and D0bar mesons, and one for both mesons.

    Args:
        D0 (array): The simulated D0 mesons that meet all of the criteria.
        D0bar (array): The simulated D0bar mesons that meet all of the criteria.
        both (array): Array containing the pythia simulated data that has passed the selection and multiple candidate requirements.
        path (string): The file path parsed in the path parser.
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
        data (array): Array containing the pythia simulated data that has passed the selection and multiple candidate requirements.

    Returns:
        D0_data (array): The simulated D0 mesons that meet all of the criteria.
        D0bar_data (array): The simulated D0bar mesons that meet all of the criteria.
        data (array): Array containing the pythia simulated data that has passed the selection and multiple candidate requirements.
    """
    length = len(data)
    PID = data[:, 1]
    
    # Define a mask array of same length as the data with all ones.
    mask = np.ones(length)
    
    # Returns True for an element in the D0 and D0bar mask arrays if the PID number is equal to 421 for D0 and -421 for D0bar.
    mask_D0 = np.logical_and(mask, PID==421)
    mask_D0bar = np.logical_and(mask, PID==-421)
    
    # split the simulated particles into mesons and antimesons.
    D0_data = data[mask_D0]
    D0bar_data = data[mask_D0bar]

    print('checkpoint: Data has been split into D0 and D0bar subsets')
        
    return D0_data, D0bar_data, data

# - - - - - - - MAIN BODY - - - - - - - #

args=parse_arguments()

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
    save_file(D0_DATA, D0BAR_DATA, cut_data, args.path)
