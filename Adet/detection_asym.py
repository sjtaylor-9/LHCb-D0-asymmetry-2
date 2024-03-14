"""
detection_asym.py


Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)/
Last edited: 
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
        choices = ['eta', 'local', 'pT'],
        default=os.getcwd(),
        help="flag to set the binning scheme to be used"
    )
    parser.add_argument(
        "--bin",
        type=int,
        required=True,
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
    '''
    Opens .txt files and reads the values of the raw asymmetries and their uncertainties.
    
    Returns these two values for both control modes.
    '''
    # Initialize dictionaries to store Araw and Araw_err for each mode
    Araw_dict = {}
    Araw_err_dict = {}
    for mode in ['kpipi', 'kspi']:
        with open(f'{args.input}/{args.scheme}/20{args.year}/mag{args.polarity}/{args.bin}/{mode}/output_{mode}.txt') as f:
            for line in f:
                # The output of the txt file is in the format Araw +/- Araw_err % so Araw is element 0  and Araw_err is element 2
                currentline = line.split()
                Araw = float(currentline[0])
                Araw_err = float(currentline[2])
                Araw_dict[mode] = Araw
                Araw_err_dict[mode] = Araw_err
    
    return Araw_dict, Araw_err_dict

def calculate_detection_asym(A_kpipi, A_kpipi_err, A_kspi, A_kspi_err):
    """

    Args:
        A_kpipi (_type_): _description_
        A_kpipi_err (_type_): _description_
        A_kspi (_type_): _description_
        A_kspi_err (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Detection asymmetry for K0 found from a paper: https://arxiv.org/pdf/1405.2797.pdf
    Adet_k0 = 0.054
    Adet_err_k0 = 0.014

    # Calculate the kpi detection asymmetry and the associated error
    Adet = A_kpipi - A_kspi - Adet_k0
    Adet_err = np.sqrt(((A_kpipi_err)**2+(A_kspi_err)**2+(Adet_err_k0)**2))
    
    return Adet, Adet_err

def output_results(Adet, uncertainty):
    '''
    This function takes all the necessary values and outputs them to the screen in a nicely formatted way.
    It also outputs them to a .txt file, written in the directory established by the user.
    '''

    asymmetry = str(round(Adet, 3)) + ' +/- ' + str(round(uncertainty, 3))
    print(f'The 20{args.year} Mag{args.polarity} detection asymmetry of bin {args.bin} in the {args.mode} phase-space is:', asymmetry)
    print("------------------------------")
    
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

print("Kpipi:", Araw_kpipi, Araw_err_kpipi)
print("Kspi:", Araw_kspi, Araw_err_kspi)

# Calculaete the kpi dection asymmetry
Adet, Adet_err = calculate_detection_asym(Araw_kpipi, Araw_err_kpipi, Araw_kspi, Araw_err_kspi)

output_results(Adet, Adet_err)