"""
create_binning_scheme.py

This code calculates the best boundary positions in  order to generete the bins. It first splits the data in 10 sets according using the transverse momentum of each event, and makes sure the number of events in each set is the same. This exact same procedure is then repeated for each of these sets but now using the pseudorapidity. Therefore, the result are 100 sets of data with an equal number of events. The boundaries calculated by this code are naturally calculated using data from both mesons and both polarities.
The year of interest and size of the data to be analysed must be specified using the required flags --year --size. There also are the flags --input and --path, which are not required. These are used to specify the directory where the input data is located and where the output file should be written, respectively. By default it is set to be the current working directory.
It outputs the values of transverse momentum and pseudorapidity at each of the boundaries of the bins. The code also calculates the best boundary position for 10 eta bins and 10 pT bins seperately.

Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)/
Last edited: 15th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #
import os
import argparse
import numpy as np
import uproot
import pandas as pd
import awkward as ak


# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size  Used to specify the amount of events the user is interested in analysing.
            The argument must be one of: [1-800]. The integer must be divisible by 10. The integers specify the number of root
            files to be read in.
    --path      Used to specify the directory in which the output files should be written. It is not required,
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
        "--size",
        type=size_argument,
        required=True,
        help="flag to set the data taking year."
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
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
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
        raise argparse.ArgumentTypeError("Invalid value.Choose an integer between 1 and 800 that is divisible by 10.")

def generate_list(size_value_local):
    result_list = []
    current_value = 10

    while not result_list or result_list[-1] < size_value_local:
        result_list.append(current_value)
        current_value += 10

    return result_list
        
        
# - - - - - - - MAIN BODY - - - - - - - #
NBINS = 10
args = parse_arguments()

# import data
tree_name = "D02Kpi_Tuple/DecayTree"
size_value = args.size
if isinstance(size_value, int):
    size_value = int(size_value)
    size_list = generate_list(size_value)
    
data = uproot.concatenate((f"{args.input}/20{args.year}/{polarity}/both/{polarity}_data_{args.year}_{size}_clean.root:{tree_name}" for polarity in ["up", "down"] for size in size_list),
  expressions=["D0_MM", "D0_PT", "D0_ETA"])

# select particles with pT below 10 GeV/c

length = len(data["D0_PT"])
print(length, 'events')
mask = np.ones(length)
mask = np.logical_and(mask, data["D0_PT"]<10000)
mask = np.logical_and(mask, data["D0_ETA"]>2)
mask = np.logical_and(mask, data["D0_ETA"]<5)

length = np.sum(mask)
data = data[mask]

# create bins
bins = np.empty((NBINS+1, NBINS+1))
pT_frame = pd.DataFrame({'Values': data["D0_PT"]})
pT_res, pT_bins = pd.qcut(pT_frame['Values'], q=NBINS, retbins=True)
bins[0] = pT_bins

length = len(data["D0_PT"])
for i in np.arange(0, NBINS):
    bin_mask = np.ones(length)
    bin_mask = np.logical_and(bin_mask, pT_bins[i]<data["D0_PT"])
    bin_mask = np.logical_and(bin_mask, pT_bins[i+1]>=data["D0_PT"])
    selected_data = data[bin_mask]
    eta_frame = pd.DataFrame({'Values': selected_data["D0_ETA"]})
    eta_frame['TimeBins'], eta_bins = pd.qcut(eta_frame['Values'], q=NBINS, retbins=True)
    timeCounts = eta_frame['TimeBins'].value_counts()
    print(timeCounts)
    bins[i+1] = eta_bins

# outputbin edges 
np.savetxt(f"{args.path}/{args.year}_{args.size}_bins.txt", bins, delimiter=',')

#output bin edges for 10 pT bins
np.savetxt(f"{args.path}/{args.year}_{args.size}_pT_bins.txt", bins[0], delimiter=',')


#Creating bin edges for 10 eta bins

length_eta = len(data["D0_ETA"])
mask_eta = np.ones(length_eta)
mask_eta = np.logical_and(mask_eta, data["D0_ETA"]<5)
length_eta = np.sum(mask_eta)
data = data[mask_eta]

# create bins
bins_eta = np.empty(NBINS+1)
eta_frame = pd.DataFrame({'Values': data["D0_ETA"]})
eta_res, eta_bins = pd.qcut(eta_frame['Values'], q=NBINS, retbins=True)

# outputbin eta edges 
np.savetxt(f"{args.path}/{args.year}_{args.size}_eta_bins.txt", eta_bins, delimiter=',')