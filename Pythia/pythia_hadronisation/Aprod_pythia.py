
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
import csv
import pandas as pd
import glob

# - - - - - - - FUNCTIONS - - - - - - - #
def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    --path      Used to specify the directory in which the root files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input_bins     Used to specify the directory in which the input data should be found. It is not required,
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
        help="flag to set the path where the output files of individual results should be written to"
    )
    parser.add_argument(
        "--results_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files of final results should be written to"
    )
    parser.add_argument(
        "--input_bins",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data of each bin should be found"
    )
    parser.add_argument(
        "--input_global",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data of total data should be found"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        required=True,
        choices=["pT", "eta", "local"],
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
        
# Code to calculate asymmetry
    
def read_from_file(meson, bin_num, scheme):
    '''
    Opens a .txt file and reads the values of the signal normalization constant and its uncertainty.
    
    Returns these two values.
    '''
    bin_num = int(bin_num)
    filename = f"{args.input_bins}/number_of_events_{scheme}_{meson}.txt"
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Get the N value based on bin number
        N_value = float(lines[bin_num])  # Subtract 1 as Python uses 0-based indexing
        # Assuming you want to calculate uncertainty as the square root of N
        N_err_value = np.sqrt(N_value)
    return N_value, N_err_value
    
def get_yield(scheme, bin_num):

    N_D0, N_D0_err = read_from_file("D0", bin_num, scheme)
    N_D0bar, N_D0bar_err = read_from_file("D0bar", bin_num, scheme)


    return N_D0, N_D0_err, N_D0bar, N_D0bar_err

def count_rows(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        row_count = sum(1 for row in reader)
    return row_count

def production_asymm(N_D0, N_D0_err, N_D0bar, N_D0bar_err):

    print(f'N_D0 {N_D0}')
    print(f'N_D0bar {N_D0bar}')
    print(f'N_D0_err {N_D0_err}')
    print(f'N_D0bar_err {N_D0bar_err}')
    A_prod = (N_D0 - N_D0bar)/(N_D0 + N_D0bar)
    A_prod_err = np.sqrt((2 * N_D0bar / (N_D0 + N_D0bar)**2 * N_D0_err)**2 + (-2 * N_D0 / (N_D0 + N_D0bar)**2 * N_D0bar_err)**2)

    return 100*A_prod, 100*A_prod_err

def average_asymm(val, err):

    weighted_mean, sum_weights = np.average(np.array(val), weights=np.power(np.array(err), -2), returned=True)
    uncertainty = np.power(sum_weights, -0.5)
    return weighted_mean, uncertainty

def output_results(A_prod, A_prod_err, bin_num):
    '''
    This function takes all the necessary values and outputs them to the screen in a nicely formatted way.
    It also outputs them to a .txt file, written in the directory established by the user.
    '''

    asymmetry = str(round(A_prod, 10)) + '% +/- ' + str(round(A_prod_err, 10)) + '% (stat)'
    print(f'The binning scheme results in a prod asymmetry of bin {bin_num} is:', asymmetry)
    print("------------------------------")
    
    array = np.array([A_prod, A_prod_err])
    np.savetxt(f"{args.path}/asymmetries_pythia_{scheme}_bin{bin_num}.txt", array)

    return

def A_prod_global():
    """
    Finds asymmetry of all data, wihtout splitting into bins 'global'.
    """

    N_D0_global, N_D0_err_global, = read_from_file_global("D0")
    N_D0bar_global, N_D0bar_err_global, = read_from_file_global("D0bar")
    
    A_prod_bin_integrated, A_prod_err_bin_integrated= production_asymm(N_D0_global, N_D0_err_global, N_D0bar_global, N_D0bar_err_global)

    return A_prod_bin_integrated, A_prod_err_bin_integrated,

def read_from_file_global(meson):
    N_global = 0
    N_err_global = 0
    filename_pattern = f"{args.input_global}/{meson}_clean_pythia_data_*.csv"
    filenames = glob.glob(filename_pattern)    
    for filename in filenames:
            with open(filename):        
                N_global = count_rows(filename)
                N_err_global = np.sqrt(N_global)
                
    return N_global, N_err_global

# - - - - - - - MAIN CODE - - - - - - - #
args = parse_arguments()
scheme = args.scheme

# Initialize as empty NumPy arrays
A_prod_list = np.array([])
A_prod_err_list = np.array([])

A_prod_bin_integrated, A_prod_err_bin_integrated = A_prod_global()

asymmetry_global_string = str(round(A_prod_bin_integrated, 10)) + '% +/- ' + str(round(A_prod_err_bin_integrated, 10)) + '% (stat)'
print(f'The simulated data results in a bin integrated prod asymmetry of: ', asymmetry_global_string)

    
if scheme == "pT" or scheme == "eta":
    for j in range(0,10):
        bin_num = j
        
        N_D0, N_D0_err, N_D0bar, N_D0bar_err = get_yield(scheme, bin_num)

        # Calculate the local production asymmetries in each bin across the pT or eta phase space.
        A_prod, A_prod_err = production_asymm(N_D0, N_D0_err, N_D0bar, N_D0bar_err)
        A_prod_list = np.concatenate([A_prod_list, [A_prod]])
        A_prod_err_list = np.concatenate([A_prod_err_list, [A_prod_err]])

        output_results(A_prod, A_prod_err, bin_num)

if scheme == 'local':
    for j in range(0,10):
        for i in range (0,10):
            bin_num = str(j)+str(i)

            N_D0, N_D0_err, N_D0bar, N_D0bar_err = get_yield(scheme, bin_num)

            # Calculate the local production asymmetries in each bin across the pT or eta phase space.
            A_prod, A_prod_err = production_asymm(N_D0, N_D0_err, N_D0bar, N_D0bar_err)
            A_prod_list = np.concatenate([A_prod_list, [A_prod]])
            A_prod_err_list = np.concatenate([A_prod_err_list, [A_prod_err]])

            output_results(A_prod, A_prod_err, bin_num)

if scheme == "pt_eta":
    for j in range(0,10):
        bin_num = str(j)
        
        N_D0, N_D0_err, N_D0bar, N_D0bar_err = get_yield(scheme, bin_num)

        # Calculate the local production asymmetries in each bin across the pT or eta phase space.
        A_prod, A_prod_err = production_asymm(N_D0, N_D0_err, N_D0bar, N_D0bar_err)
        A_prod_list = np.concatenate([A_prod_list, [A_prod]])
        A_prod_err_list = np.concatenate([A_prod_err_list, [A_prod_err]])

        output_results(A_prod, A_prod_err, bin_num)

average_asymmetry_across_bins, uncertainty_in_average_asymmetry_across_bins = average_asymm(A_prod_list, A_prod_err_list)

file_path = f"{args.results_path}/final_pythia_Text_asymmetries_{args.scheme}.txt"

with open(file_path, "w") as file:
    text = (
        f'A_prod_average_across_bins: {round(average_asymmetry_across_bins, 10)}\n'
        f'A_prod_err_average_across_bins: {round(uncertainty_in_average_asymmetry_across_bins, 10)}\n'
        f'A_prod_bin_integrated: {round(A_prod_bin_integrated, 10)}\n'
        f'A_prod_err_bin_integrated: {round(A_prod_err_bin_integrated, 10)}\n'
    )
    file.write(text)

print(f"Completed asymmetry calcultions for scheme: {scheme}")