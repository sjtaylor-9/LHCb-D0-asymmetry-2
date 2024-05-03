import os
import argparse
import numpy as np
from lhcbstyle import LHCbStyle
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import ListedColormap
from matplotlib import colormaps
import awkward as ak
import sys
import matplotlib.pyplot as plt
import mplhep; mplhep.style.use("LHCb2")
from scipy.stats import chi2
from matplotlib.patches import Rectangle


# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size      Used to specify the amount of events the user is interested in analysing.
                The argument must be one of: [1-800]. The interger must be divisible by 10. The integers specify the number of root
                files to be read in.
    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --bin_path  Used to specify the directory in which the binning scheme should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --scheme    
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
        "--input",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the detection asymmetries should be found"
    )
    parser.add_argument(
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--bin_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the binning scheme should be found"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        choices=["pT","eta"],
        required=True,
        help="flag to set whether a binned or an unbinned should be performed (y/n)"
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

def read_asymmetry_values(polarity):
    with open(f'{args.input}/global/{args.year}/{polarity}/detection_asym_{args.year}_{polarity}.txt') as f:
        for line in f:
            currentline = line.split()
            unbinned_asymm = float(currentline[0])
            unbinned_asymm_error = float(currentline[2])
        f.close()
    return unbinned_asymm, unbinned_asymm_error

def detection_asym(A_det_up, A_det_down, A_det_up_err, A_det_down_err):
    """
    This function calculates the total detection asymmetry and it's uncertainty from the up and down components that are outputted from detection_asym.py.
    The detection asymmetry for Kpi is defined as A_detKpi) = A_raw(D->Kpipi) - A_raw(D->Ks0pi) - A_det(Ks0).

    Args:
        A_det_up (float): The MagUp detection asymmetry outputted from detection_asym.py.
        A_det_up_err (float): The MagUp detection asymmetry uncertainty outputted from detection_asym.py.
        A_det_down (float):The MagDown detection asymmetry outputted from detection_asym.py.
        A_det_down_err (float): The MagDown detection asymmetry outputted from detection_asym.py.

    Returns:
        A_det: The polarity-integrated detection asymmetry.
        A_det_err: The total uncertainty associated with the detection asymmetry.
    """
    
    # Total detection asymmetry is the arithmetic mean of the up/down asymmetries and the total uncertainty is the up/down uncertainties in quadrature
    A_det = (A_det_up + A_det_down) / 2
    A_det_error = np.sqrt(((A_det_up_err)**2+(A_det_down_err)**2)) / 2

    return A_det, A_det_error
# - - - - - - - MAIN BODY - - - - - - - - - #
args = parse_arguments()

unbinned_asymm_up, unbinned_asymm_error_up = read_asymmetry_values('up')
unbinned_asymm_down, unbinned_asymm_error_down = read_asymmetry_values('down')

unbinned_asymm, unbinned_asymm_error = detection_asym(unbinned_asymm_up, 
                                                  unbinned_asymm_down, 
                                                  unbinned_asymm_error_up, 
                                                  unbinned_asymm_error_down)

Adet = np.array([])
Adet_error = np.array([])

for j in range(0,10):
    bin_num = str(j)
    with open(f'{args.input}/{args.scheme}/{args.year}/up/detection_asym_{args.year}_up_bin{bin_num}.txt') as f:
         for line in f:
            currentline = line.split()
            value_up = float(currentline[0])
            error_up = float(currentline[2])
    f.close()
    with open(f'{args.input}/{args.scheme}/{args.year}/down/detection_asym_{args.year}_down_bin{bin_num}.txt') as f:
        currentline = line.split()
        value_down = float(currentline[0])
        error_down = float(currentline[2])
    f.close()
    detection_asymmetries = detection_asym(value_up, value_down, error_up, error_down)
    Adet = np.append(Adet, detection_asymmetries[0])
    Adet_error = np.append(Adet_error, detection_asymmetries[1])
    
x_value =np.array([])
x_value_error = np.array([])

file_path = f"{args.bin_path}/{args.year}_{args.size}_{args.scheme}_bins.txt"  # Replace with the path to your file
# Open the file in read mode
with open(file_path, 'r') as file:
    # Read all lines from the file and store them in a list
    bin_lines = [float(line.strip()) for line in file.readlines()]

print(bin_lines)
for i in range(0,10):
    # Get center of bin and width of bin as error
    x_value_indivual = (bin_lines[i]+bin_lines[i+1])/2
    x_value_error_indivual = (bin_lines[i+1]) - ((bin_lines[i]+bin_lines[i+1])/2)
    x_value_error = np.append(x_value_error, x_value_error_indivual)
    x_value = np.append(x_value, x_value_indivual)


if args.scheme == 'pT':
    x_value = [x / 1000 for x in x_value] # in Gev
    x_value_error = [x / 1000 for x in x_value_error] # in Gev

bins = len(x_value)

# Plotting
fig = plt.figure()
ax1 = fig.gca()

ax1.set_ylabel(r'$A_{\mathrm{det}}$ [%]', fontsize = 50)
ax1.tick_params(axis='both', which='both', labelsize=30)
if args.scheme == 'pT':
    ax1.set_xlabel(r'$p_{T}$ [GeV$c^{-1}$]', fontsize = 50)
elif args.scheme == 'eta':
    ax1.set_xlabel(r'$\eta$', fontsize = 50)

line2 = ax1.axhline(unbinned_asymm, color='red', linestyle='solid', linewidth=5)
fill2 = ax1.axhspan(unbinned_asymm-unbinned_asymm_error, unbinned_asymm+unbinned_asymm_error, color='red', alpha=0.35, lw=0)

Data = ax1.errorbar(x_value, Adet, yerr=Adet_error,xerr=x_value_error, fmt='o', capsize=5, color = 'black', label = 'Data')

if args.scheme == 'pT':
    ax1.legend([(line2,fill2),Data],[r'Bin integrated result','Data'])#, loc='upper right')
elif args.scheme == 'eta':
    ax1.legend([(line2,fill2),Data],[r'Bin integrated result','Data'])#,='upper right')

if args.scheme == 'pT':
    plt.savefig(f'{args.path}/pT_Adet_{args.year}.pdf', bbox_inches = "tight")
elif args.scheme == 'eta':
    plt.savefig(f'{args.path}/eta_Adet_{args.year}.pdf', bbox_inches = "tight")

plt.show()