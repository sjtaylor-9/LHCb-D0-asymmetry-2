
"""
Aprod_pythia.py

This code is used to process the signal normalization yields and obtain the production asymmetries of the average bins and global integrated. It finally outputs the results obtained both to the secreen and to a .txt file.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size --scheme --blind. There also are the flags --input --seedval --path which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.
This code is  inspired on the work of Camille Jarvis-Stiggants and Michael England and Marc Oriol Pérez. The code has been completely rewritten and reorganised, and some features have been added to add flexibility to the code, but some of the original functions have been used here as well.

Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 17th February 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import os
import argparse
import numpy as np

# - - - - - - - FUNCTIONS - - - - - - - #
def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    --path      Used to specify the directory in which the root files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --scheme    
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
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found"
    )
    parser.add_argument(
    "--scheme",
    type=str,
    required=True,
    choices=["pT", "rapidity", "local"],
    default=os.getcwd(),
    help="flag to deduce what binning scheme is to be used"
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
        
def read_from_file(scheme, bin_num):
    '''
    Opens a .csv files and reads the values of the signal normalization constant and its uncertainty.
    
    Returns these two values.
    '''
    data = np.genfromtxt(f'{args.input}/both_{scheme}_bin{bin_num}.csv', delimiter = ',')
    PID = data[:, 1]
    
    return data, PID

def split_meson(data, PID):
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
    # Define a mask array of same length as the data with all ones.
    mask = np.ones(length)
    
    # Returns True for an element in the D0 and D0bar mask arrays if the PID number is equal to 421 for D0 and -421 for D0bar.
    mask_D0 = np.logical_and(mask, PID==421)
    mask_D0bar = np.logical_and(mask, PID==-421)
    
    # split the simulated particles into mesons and antimesons.
    D0_data = data[mask_D0]
    D0bar_data = data[mask_D0bar]
        
    return D0_data, D0bar_data

def production_asymmetry(D0, D0bar):
    
    # Number of D0 and D0bar mesons generated are the lengths of D0 and D0bar respectively.
    n_D0 = len(D0)
    n_D0bar = len(D0bar)
    
    Aprod = (n_D0 - n_D0bar) / (n_D0 + n_D0bar)
    # The error in the number of simulated events is determined by Poissonian statistics, sqrt(N).
    n_D0_err = (np.sqrt(n_D0))
    n_D0bar_err = (np.sqrt(n_D0bar))
    
    A_err = 2*(((n_D0bar**2)*(n_D0_err**2) + (n_D0**2)*(n_D0bar_err**2))**0.5)*((n_D0 + n_D0bar)**(-2))

    return 100*Aprod, 100*A_err

def integrated_asym(val, err):

    weighted_mean, sum_weights = np.average(np.array(val), weights=np.power(np.array(err), -2), returned=True)
    uncertainty = np.power(sum_weights, -0.5)
    return weighted_mean, uncertainty
# - - - - - - - MAIN CODE - - - - - - - #

args = parse_arguments()
scheme = args.scheme

local_production_asyms = []
local_production_errs = []

if scheme == 'local':
    for j in range(0,10):
        for i in range (0,10):
            bin_num = str(j)+str(i)
            
            data, PID = read_from_file(scheme, bin_num)
            D0, D0bar = split_meson(data, PID)
            # Calculate the local production asymmetries across the (pT, eta) phase space.
            Aprod_local = production_asymmetry(D0, D0bar)
            
            local_production_asyms.append(Aprod_local[0])
            local_production_errs.append(Aprod_local[1])


    global_Aprod = integrated_asym(local_production_asyms, local_production_errs)
    print(np.mean(local_production_asyms))
    print(f"The global bin-integrated production asymmetry for the simulated data is: ", round(global_Aprod[0], 3), "% +/-", round(global_Aprod[1], 3), '%')

    # Aprod, Aprod error -- both from weighted mean
    array = np.array([global_Aprod[0], global_Aprod[1]])
    np.savetxt(f"{args.path}/final_asymmetries_{scheme}.txt", array)

    file_path = f"{args.path}/final_text_asymmetries_{scheme}.txt"

    with open(file_path, "w") as file:
        text = (
            f'A_prod_intergrated: {round(global_Aprod[0], 3)}\n'
            f'A_prod_err_intergrated: {round(global_Aprod[1], 3)}\n'
        )
        file.write(text)

    print(100*((149596-150530)/(149596+150530)))
    
elif scheme == 'pT' or scheme == 'rapidity':
    for j in range(0,10):
        bin_num = str(j)
        
        data, PID = read_from_file(scheme, bin_num)
        D0, D0bar = split_meson(data, PID)
        # Calculate the local production asymmetries across the pT or eta phase space.
        Aprod_local = production_asymmetry(D0, D0bar)
        local_production_asyms.append(Aprod_local[0])
        local_production_errs.append(Aprod_local[1])
    
    print(np.mean(local_production_asyms))
