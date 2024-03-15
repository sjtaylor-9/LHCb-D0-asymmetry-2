"""
detection_asym.py

This code calculates the Kpi detection asymmetry for a given year and polarity in each of the bins in the pT, eta and (pT, eta) binning schemes. The Kpi detection asymmetry is defined as A_det(Kpi) = A_raw(D->K-pi+pi+) - A_raw(D->Ks0pi+) - A_det(Ks0).
The raw Kpipi and Kspi asymmetries are read in from .txt files saved on the eos, which are outputs from the HTCondor jobs in the kinematic reweighting procedure.
The year of interest, polarity, binning scheme and bin number to be analysed must be specified using the required flags --year --polarity --scheme --bin. There is also the --input flag to set the eos path for where the Kpipi and Kspi raw asymmetries are saved. The --path flag sets the output directory which is not required.

Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)/
Last edited: 15th March 2024
"""
# - - - - - - IMPORT STATEMENTS - - - - - - #
import os
import argparse
import numpy as np
# - - - - - - - FUNCTIONS - - - - - - - - - #
def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --polarity  Used to specify the polarity of the magnet the user is interested in.
                The argument must be one of: [up, down].
    --path      Used to specify the directory in which the root files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --scheme    Used to specify the binning scheme to be analysed.
                The argument must be one of: [pT, eta, local].
    --bin       Used to specify the bin number to be analysed.
                For the pT, eta schemes this is [0, 9] and for the local scheme this is [0, 99].
    
    Returns the parsed arguments.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--year",
        type=int,
        choices=[16,17,18],
        required=True,
        help="flag to set the data taking year."
    )
    parser.add_argument(
        "--polarity",
        type=str,
        choices=["up","down"],
        required=True,
        help="flag to set the data taking polarity."
    )
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
        required=True,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        required=True,
        choices = ['eta', 'local', 'pT', 'global'],
        default=os.getcwd(),
        help="flag to set the binning scheme to be used"
    )
    parser.add_argument(
        "--bin",
        type=int,
        required=False,
        default=os.getcwd(),
        help="flag to set the bin number"
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

def read_from_file():
    """
    Opens .txt files and reads the values of the raw asymmetries and their uncertainties.
    If args.scheme is set as global then the file path is set to be the global one.
    Returns these two values for both control modes.

    Returns:
        A_raw_dict: A dictionary containing the raw asymmetries for both control modes. The specific asymmetries are defined as A_raw_dict['kpipi'] and A_raw_dict['A_raw_kspi'].
    """
    # Initialize dictionaries to store Araw and Araw_err for each mode
    Araw_dict = {}
    Araw_err_dict = {}
    
    # Opens the .txt files for both modes and saves the asymmetry and its uncertainty to their respective dictionaries. 
    for mode in ['kpipi', 'kspi']:
        # Determines the file path based on if the code is being ran for the global or local detection asymmetries
        if args.scheme == 'global':
            with open(f'{args.input}/20{args.year}/mag{args.polarity}/{mode}/output_{mode}.txt') as f:
                for line in f:
                    # The output of the .txt file is in the format Araw +/- Araw_err % so Araw is element 0  and Araw_err is element 2
                    currentline = line.split()
                    Araw = float(currentline[0])
                    Araw_err = float(currentline[2])
                    Araw_dict[mode] = Araw
                    Araw_err_dict[mode] = Araw_err
        else:
            with open(f'{args.input}/{args.scheme}/20{args.year}/mag{args.polarity}/{args.bin}/{mode}/output_{mode}.txt') as f:
                for line in f:
                    # The output of the .txt file is in the format Araw +/- Araw_err % so Araw is element 0  and Araw_err is element 2
                    currentline = line.split()
                    Araw = float(currentline[0])
                    Araw_err = float(currentline[2])
                    Araw_dict[mode] = Araw
                    Araw_err_dict[mode] = Araw_err
            
        
    return Araw_dict, Araw_err_dict

def calculate_detection_asym(A_kpipi, A_kpipi_err, A_kspi, A_kspi_err):
    """
    This function calculates the Kpi detection asymmetry and its uncertainty for a given year, polarity and bin number within the chosen binning scheme.
    
    Args:
        A_kpipi (float): The raw asymmetry for the Kpipi control mode.
        A_kpipi_err (float): The raw uncertainty associated with the Kpipi raw asymmetry.
        A_kspi (float): The raw asymmetry for the Kspi control mode.
        A_kspi_err (float): The raw uncertainty associated with the Kspi raw asymmetry.

    Returns:
        Adet (float): The Kpi detection asymmetry.
        Adet_err (float): The detection asymmetry uncertainty asssociated with the Kpi detection asymmetry.
    """
    # Detection asymmetry for K0 found from a paper: https://arxiv.org/pdf/1405.2797.pdf
    Adet_k0 = 0.054
    Adet_err_k0 = 0.014

    # Calculate the kpi detection asymmetry and the associated error in quadrature
    Adet = A_kpipi - A_kspi - Adet_k0
    Adet_err = np.sqrt(((A_kpipi_err)**2+(A_kspi_err)**2+(Adet_err_k0)**2))
    
    return Adet, Adet_err

def output_results(Adet, uncertainty):
    """
    This function takes all the necessary values and outputs them to the screen in a nicely formatted way.
    It also outputs them to a .txt file, written in the directory established by the user.

    Args:
        Adet (float): The Kpi detection asymmetry for a given year, polarity and bin number.
        uncertainty (float): The detection asymmetry uncertainty associated with Adet.
    """
    if args.scheme == 'global':
        # Print the detection asymmetry rounded to 3 sf
        asymmetry = str(round(Adet, 3)) + ' +/- ' + str(round(uncertainty, 3))
        print(f'The global 20{args.year} Mag{args.polarity} detection asymmetry is:', asymmetry)
        # Save the detection asymmetry and its error to a .txt file
        array = np.array([Adet, uncertainty])
        np.savetxt(f"{args.path}/detection_asym_{args.year}_{args.polarity}.txt", array)
    else:
        # Print the detection asymmetry rounded to 3 sf
        asymmetry = str(round(Adet, 3)) + ' +/- ' + str(round(uncertainty, 3))
        print(f'The 20{args.year} Mag{args.polarity} detection asymmetry of bin {args.bin} in the {args.scheme} binning scheme is:', asymmetry)
        # Save the detection asymmetry and its error to a .txt file
        array = np.array([Adet, uncertainty])
        np.savetxt(f"{args.path}/detection_asym_{args.year}_{args.polarity}_bin{args.bin}.txt", array)
# - - - - - - - MAIN BODY - - - - - - - - - #

args = parse_arguments()

A_raw, Araw_err = read_from_file()
# Assign the dictionaries to individual variables
Araw_kpipi = A_raw['kpipi']
Araw_err_kpipi = Araw_err['kpipi']
Araw_kspi = A_raw['kspi']
Araw_err_kspi = Araw_err['kspi']

# Calculaete the kpi dection asymmetry
Adet, Adet_err = calculate_detection_asym(Araw_kpipi, Araw_err_kpipi, Araw_kspi, Araw_err_kspi)

output_results(Adet, Adet_err)