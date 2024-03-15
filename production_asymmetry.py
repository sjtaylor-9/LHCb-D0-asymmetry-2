
"""
production_asymmetry.py

This code is used to process the signal normalization yields and obtain the production asymmetries of the average bins and global integrated. It finally outputs the results obtained both to the secreen and to a .txt file.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size --scheme --blind. There also are the flags --input --seedval --path which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.
This code is  inspired on the work of Camille Jarvis-Stiggants and Michael England and Marc Oriol Pérez. The code has been completely rewritten and reorganised, and some features have been added to add flexibility to the code, but some of the original functions have been used here as well.

Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)
Last edited: 15th March 2024
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #
import random
import os
import argparse
import numpy as np
import seaborn as sns
# - - - - - - - FUNCTIONS - - - - - - - - - #
def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year          Used to specify the year at which the data was taken the user is interested in.
                    The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size          Used to specify the amount of events the user is interested in analysing.
                    The argument must be one of: [1-800]. The interger must be divisible by 10. The integers specify the number of root
                    files to be read in.
    --polarity      Used to specify the polarity of the magnet the user is interested in.
                    The argument must be one of: [up, down].
    --path          Used to specify the directory in which the output files should be written. It is not required,
                    in the case it is not specified, the default path is the current working directory.
    --model_input   Used to specify the directory in which the input data for the raw asymmetries should be found.
    --detection_input   Used to specify the directory in which the input dectection asymmetries should be found.
    --scheme        Used to specify the binning scheme to be analysed.
                    The argument must be one of: [pT, eta, local].
    --bin           Used to specify the bin number to be analysed.
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
        "--size",
        type=size_argument,
        required=True,
        help="flag to set the size of the input data."
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
        help="flag to set the path where the output files should be written to."
    )
    parser.add_argument(
        "--model_input",
        type=dir_path,
        required=True,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found."
    )
    parser.add_argument(
        "--detection_input",
        type=dir_path,
        required=True,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found."
    )
    parser.add_argument(
        "--scheme",
        type=str,
        required=True,
        choices = ['eta', 'local', 'pT'],
        default=os.getcwd(),
        help="flag to set the binning scheme to be used."
    )
    parser.add_argument(
        "--bin",
        type=int,
        required=True,
        default=os.getcwd(),
        help="flag to set the bin number."
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

def size_argument(value):
    if value.isdigit():
        # If the input is a digit, treat it as an integer
        int_value = int(value)
        if 1 <= int_value <= 800 and int_value % 10 == 0: 
            return int_value
        else:
            raise argparse.ArgumentTypeError("Integer value must be between 1 and 800 and be divisible by 10.")
    else:
        raise argparse.ArgumentTypeError("Invalid value. Choose between 'small', 'medium', 'large', or an integer between 1 and 800 that is divisible by 10.")

def read_from_file(polarity, bin_num=None, scheme=None, meson=None):
    """
    Opens a .txt files and reads the values of the asymmetry and its uncertainty.
    Returns these two values.

    Args:
        polarity (_type_): _description_
        bin_num (_type_, optional): _description_. Defaults to None.
        scheme (_type_, optional): _description_. Defaults to None.
        meson (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    if bin_num is not None:
        if meson is not None:
            with open(f'{args.model_input}/{scheme}/{bin_num}/yields_{meson}_{polarity}_{args.year}_{args.size}_bin{bin_num}.txt') as f:
                for line in f:
                    # The output of the .txt file is in the format signal_yield, signal_yield_err % so the signal yield is element 0  and uncertainty is element 1
                    currentline = line.split(",")
                    Output = float(currentline[0])
                    Output_err = float(currentline[1])
                f.close()
        else:
            with open(f'{args.detection_input}/{scheme}/{args.year}/{polarity}/detection_asym_{args.year}_{polarity}_bin{bin_num}.txt') as f:
                for line in f:
                    # The output of the .txt file is in the format Adet +/- Adet_err % so Adet is element 0  and Adet_err is element 2
                    currentline = line.split()
                    Output = float(currentline[0])
                    Output_err = float(currentline[2])
                f.close()
    else:
        if meson is not None:
            with open(f'{args.model_input}/global/yields_{meson}_{polarity}_{args.year}_{args.size}.txt') as f:
                for line in f:
                    # The output of the .txt file is in the format signal_yield, signal_yield_err % so the signal yield is element 0  and uncertainty is element 1
                    currentline = line.split(",")
                    Output = float(currentline[0])
                    Output_err = float(currentline[1])
                f.close()
        else:
            with open(f'{args.detection_input}/global/{args.year}/{polarity}/detection_asym_{args.year}_{args.polarity}.txt"') as f:
                for line in f:
                    # The output of the .txt file is in the format Adet +/- Adet_err % so Adet is element 0  and Adet_err is element 2
                    currentline = line.split()
                    Output = float(currentline[0])
                    Output_err = float(currentline[2])
                f.close()
    return Output, Output_err

def get_yield(bin_num, scheme):
    """
    Gets all the normalization yields, and their uncertainties, necessary to calculate the raw asymmetries.
    This takes into account both D0 and D0bar and both magnet polarities.
    
    Returns all the signal normalization constants, together with their uncertainties.

    Args:
        bin_num (_type_): _description_
        scheme (_type_): _description_

    Returns:
        _type_: _description_
    """
    yield_D0_up = read_from_file('up', bin_num, scheme, 'D0')
    yield_D0bar_up = read_from_file('up', bin_num, scheme, 'D0bar')
    yield_D0_down = read_from_file('down', bin_num, scheme, 'D0')
    yield_D0bar_down = read_from_file('down', bin_num, scheme, 'D0bar')
   
    return yield_D0_up[0], yield_D0_up[1], yield_D0bar_up[0], yield_D0bar_up[1], yield_D0_down[0], yield_D0_down[1], yield_D0bar_down[0], yield_D0bar_down[1]


def A_Det(bin_num, scheme, A_det_up, A_det_up_err, A_det_down, A_det_down_err):
    """
    The detection asymmetry for Kpi is defined as A_detKpi) = A_raw(D->Kpipi) - A_raw(D->Ks0pi) - A_det(Ks0)
    
    Args:
        bin_num (int): The bin number of the binning scheme being used.
        scheme (string): The binning scheme being used.
        Adet_global (boolean): Deduces if the detection asymmetry is the global one or not.
    Returns:
        : _description_
    """
    
    if bin_num is not None:
        # Read in the detection asymmetries from the .txt files outputted from detection_asym.py
        A_det_up, A_det_up_err = read_from_file('up', bin_num, scheme)
        A_det_down, A_det_down_err = read_from_file('down', bin_num, scheme)

    # Total detection asymmetry is the arithmetic mean of the up/down asymmetries and the total uncertainty is the up/down uncertainties in quadrature
    A_det = (A_det_up + A_det_down) / 2
    A_det_error = np.sqrt(((A_det_up_err)**2+(A_det_down_err)**2)) / 2

    return A_det, A_det_error

def calculate_raw_asymmetry(yeild_D0, yield_D0bar, bin_width, N_D0_err, N_D0bar_err):
    """
    It takes the signal yields for D0 and D0bar as arguments and then calculates the raw
    asymmetries from these. It also propagates the uncertainties.
    
    Returns both the asymmetry and its uncertainty as a percentage.

    Args:
        yeild_D0 (_type_): _description_
        yield_D0bar (_type_): _description_
        bin_width (_type_): _description_
        N_D0_err (_type_): _description_
        N_D0bar_err (_type_): _description_

    Returns:
        _type_: _description_
    """
    N_D0 = abs(yeild_D0)/abs(bin_width)
    N_D0bar = abs(yield_D0bar)/abs(bin_width)
    
    A = (N_D0 - N_D0bar)/(N_D0 + N_D0bar)
    A_err = 2*(((N_D0bar**2)*(N_D0_err**2) + (N_D0**2)*(N_D0bar_err**2))**0.5)*((N_D0 + N_D0bar)**(-2))
          
    return 100*A, 100*A_err

def output_results(A_raw, A_raw_err, bin_num, A_prod, A_prod_err):
    """
    This function takes all the necessary values and outputs them to the screen in a nicely formatted way.
    It also outputs them to a .txt file, written in the directory established by the user.

    Args:
        A_raw (_type_): _description_
        A_raw_err (_type_): _description_
        bin_num (_type_): _description_
        A_prod (_type_): _description_
        A_prod_err (_type_): _description_
    """
    asymmetry = str(round(A_raw, 3)) + ' +/- ' + str(round(A_raw_err, 3)) + ' (stat) +/- '
    print(f'The 20{args.year} raw asymmetry of bin {bin_num} is:', asymmetry)
    print("------------------------------")
    
    array = np.array([A_prod, A_prod_err, A_raw, A_raw_err])
    np.savetxt(f"{args.path}/asymmetries_{args.year}_bin{bin_num}.txt", array)

def production_asymm(raw_asym, raw_error, detection_asym, detection_error):
    """

    Args:
        raw_asym (_type_): _description_
        raw_error (_type_): _description_
        detection_asym (_type_): _description_
        detection_error (_type_): _description_

    Returns:
        _type_: _description_
    """


    A_prod = raw_asym - detection_asym
    A_prod_err = np.sqrt(((raw_error)**2 + (detection_error)**2))


    return A_prod, A_prod_err


def weighted_mean(val, err):
    """

    Args:
        val (_type_): _description_
        err (_type_): _description_

    Returns:
        _type_: _description_
    """
    weighted_mean, sum_weights = np.average(np.array(val), weights=np.power(np.array(err), -2), returned=True)
    uncertainty = np.power(sum_weights, -0.5)
    return weighted_mean, uncertainty

def A_prod_unbinned():
    """_summary_

    Returns:
        _type_: _description_
    """
    # Read in the signal yields
    yield_D0_up = read_from_file('up', 'D0')
    yield_D0bar_up = read_from_file('up', 'D0bar')
    yield_D0_down = read_from_file('down', 'D0')
    yield_D0bar_down = read_from_file('down', 'D0bar')
    
    # Read in the global detection asymmetry for magup and magdown
    A_det_up_global, A_det_up_err_global = read_from_file('up')
    A_det_down_global, A_det_down_err_global = read_from_file('down')
    
    # Calculate the global detection asymmetry.
    A_det_global, A_det_err_global = A_Det(A_det_up = A_det_up_global, 
                                           A_det_up_err = A_det_up_err_global, 
                                           A_det_down = A_det_down_global, 
                                           A_det_down_err = A_det_down_err_global, 
                                           bin_num = None, 
                                           scheme = None
                                           )

    # Calculate the global up/down raw asymmetries from the signal yields
    A_raw_up_global, A_raw_up_err_global = calculate_raw_asymmetry(yield_D0_up[0], yield_D0bar_up[0], 1, yield_D0_up[1], yield_D0bar_up[1])
    A_raw_down_global, A_raw_down_err_global = calculate_raw_asymmetry(yield_D0_down[0], yield_D0bar_down[0], 1, yield_D0_down[1], yield_D0bar_down[1])
    
    # Calculating the total raw asymmetry of the bin as the arithmetic mean of the up/down asymmetries
    A_raw_global = (A_raw_up_global + A_raw_down_global) / 2
    A_raw_err_global = np.sqrt(((A_raw_up_err_global)**2 + (A_raw_down_err_global)**2)) / 2

    
    # Calculate the global production asymmetry
    A_prod_global= production_asymm(A_raw_global, A_raw_err_global, A_det_global, A_det_err_global)

    return A_prod_global
# - - - - - - - MAIN CODE - - - - - - - - - #

args = parse_arguments()
scheme = args.scheme

A_prod_list = []
A_prod_err_list = []

Aprod_unbinned = A_prod_unbinned()

if scheme == 'local':
    for j in range(0,10):
        for i in range (0,10):
            scheme = 'local'
            bin_num = str(j)+str(i)
            
            # Get signal yield from desired model 
            N_D0_up, N_D0_up_err, N_D0bar_up, N_D0bar_up_err, N_D0_down, N_D0_down_err, N_D0bar_down, N_D0bar_down_err = get_yield(bin_num,scheme)

            # Get raw asymmetries for main model
            A_raw_up, A_raw_up_err = calculate_raw_asymmetry(N_D0_up, N_D0bar_up, 1, N_D0_up_err, N_D0bar_up_err)
            A_raw_down, A_raw_down_err = calculate_raw_asymmetry(N_D0_down, N_D0bar_down, 1, N_D0_down_err, N_D0bar_down_err)
            # Calculating the total raw asymmetry of the bin as the arithmetic mean of the up/down asymmetries
            A_raw = (A_raw_up + A_raw_down) / 2
            A_raw_err = np.sqrt(((A_raw_up_err)**2 + (A_raw_down_err)**2)) / 2

            # Calculate the detection asymmetry of the bin
            A_det, A_det_err = A_Det(bin_num=bin_num, 
                                     scheme=scheme,
                                     A_det_up = None, 
                                     A_det_up_err = None, 
                                     A_det_down = None, 
                                     A_det_down_err = None
                                    )

            # Calculate the production asymmetry of the bin
            A_prod_bin = production_asymm(A_raw, A_raw_err, A_det, A_det_err)

            # Output the raw asymmetry results
            output_results(A_raw, A_raw_err, bin_num, A_prod_bin[0], A_prod_bin[1])
            # Appending the production asymmetry of the bin to an array
            A_prod_list.append(A_prod_bin[0])
            A_prod_err_list.append(A_prod_bin[1])

elif scheme == 'pT' or scheme == 'eta':
    for j in range(0,10):
        bin_num = str(j)
        # Get signal yield from desired model 
        N_D0_up, N_D0_up_err, N_D0bar_up, N_D0bar_up_err, N_D0_down, N_D0_down_err, N_D0bar_down, N_D0bar_down_err = get_yield(bin_num,scheme)

        # Get raw asymmetries for main model
        A_raw_up, A_raw_up_err = calculate_raw_asymmetry(N_D0_up, N_D0bar_up, 1, N_D0_up_err, N_D0bar_up_err)
        A_raw_down, A_raw_down_err = calculate_raw_asymmetry(N_D0_down, N_D0bar_down, 1, N_D0_down_err, N_D0bar_down_err)
        # Calculating the total raw asymmetry of the bin as the arithmetic mean of the up/down asymmetries
        A_raw = (A_raw_up + A_raw_down) / 2
        A_raw_err = np.sqrt(((A_raw_up_err)**2 + (A_raw_down_err)**2)) /2

        # Calculate the detection asymmetry of the bin
        A_det, A_det_err = A_Det(bin_num, scheme)
        
        # Calculate the production asymmetry of the bin
        A_prod_bin = production_asymm(A_raw, A_raw_err, A_det, A_det_err)

        # Output the raw asymmetry results
        output_results(A_raw, A_raw_err, bin_num, A_prod_bin[0], A_prod_bin[1])
        # Appending the production asymmetry of the bin to an array
        A_prod_list.append(A_prod_bin[0])
        A_prod_err_list.append(A_prod_bin[1])


# Calculates the global production asymmetry as the weighted mean across the bins
AProd = weighted_mean(A_prod_list, A_prod_err_list)
print(f"The 20{args.year} global production asymmetry found from the weighted mean is: ", round(AProd[0],3), "% +/-", round(AProd[1], 3), '%')

print(f"The 20{args.year} global bin-integrated production asymmetry is: ", round(Aprod_unbinned[0],3), "% +/-", round(Aprod_unbinned[1], 3), '%')

# Saving Aprod calculated from weighted average and Aprod calculated from unbinned average of D0_up, D0_down, D0bar up, D0bar down
# Aprod, Aprod error -- both from weighted mean, Aprod from unbinned average, Aprod err from unbinned average
array = np.array([AProd[0],
                  AProd[1],
                  Aprod_unbinned[0], 
                  Aprod_unbinned[1]])
np.savetxt(f"{args.results_path}/final_asymmetries_{args.scheme}_{args.year}_{args.size}.txt", array)

file_path = f"{args.results_path}/final_text_asymmetries_{args.scheme}_{args.year}_{args.size}.txt"

with open(file_path, "w") as file:
    text = (
        f'A_prod_bin: {round(AProd[0], 3)}\n'
        f'A_prod_err_bin: {round(AProd[1], 3)}\n'
        f'A_prod_integrated: {round(Aprod_unbinned[0], 3)}\n'
        f'A_prod_err_integrated: {round(Aprod_unbinned[1], 3)}\n'
    )
    file.write(text)