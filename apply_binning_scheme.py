"""
apply_binning_scheme.py

This code applies the given binning scheme to a set of data, and generates root files cotaining the events in each bin.
The year of interest, size of the data, polarity and meson to be analysed must be specified using the required flags --year --size --polarity --meson. There also are the flags --input --path and --bin_path, which are not required. These are used to specify the directory where the input data is located, where the binning scheme can be found and where the output file should be written, respectively. By default it is set to be the current working directory.
It outputs the root files with the events in each individual bin, as well as a txt file with the number of events in each bin.

Author: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (lseelan@student.manchester.ac.uk)
Last edited: 25th March 2024
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
    --size      Used to specify the amount of events the user is interested in analysing.
                The argument must be one of: [1-800]. The integer must be divisible by 10. The integers specify the number of root
                files to be read in.
    --polarity  Used to specify the polarity of the magnet the user is interested in.
                The argument must be one of: [up, down].
    --meson     Used to specify the meson the user is interested in.
                The argument must be one of: [D0, D0bar, both].
    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --bin_path  Used to specify the directory in which the binning scheme should be found. It is not required,
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
        "--polarity",
        type=str,
        choices=["up","down"],
        required=True,
        help="flag to set the data taking polarity."
    )
    parser.add_argument(
        "--meson",
        type=str,
        choices=["D0","D0bar", "both"],
        required=True,
        help="flag to set the D0 meson flavour."
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
        help="flag to set the path where the input data should be found"
    )
    parser.add_argument(
        "--bin_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the binning scheme should be found"
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

args = parse_arguments()

# import data
tree_name = "D02Kpi_Tuple/DecayTree"
size_value = args.size
if isinstance(size_value, int):
    size_value = int(size_value)
    size_list = generate_list(size_value)

if args.meson=="both":
    data = uproot.concatenate((f"{args.input}/20{args.year}/{args.polarity}/both/{args.polarity}_data_{args.year}_{size}_clean.root:{tree_name}" for size in size_list),
    expressions=["D0_PT", "D0_ETA", "P1_ETA", "P2_ETA", "P1_PT", "P2_PT", "P1_PHI", "P2_PHI"])
else:
    data = uproot.concatenate((f"{args.input}/20{args.year}/{args.polarity}/{args.meson}/{args.meson}_{args.polarity}_data_{args.year}_{size}_clean.root:{tree_name}" for size in size_list),
    expressions=["D0_MM", "D0_PT", "D0_ETA"])

bins = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_bins.txt", delimiter=',')
bins_pT = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_pT_bins.txt", delimiter=',')
bins_eta = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_eta_bins.txt", delimiter=',')

# select particles with pT below 10 GeV/c
nevents=np.empty(0)
nevents_pT =np.empty(0)
nevents_eta =np.empty(0)

length = len(data["D0_PT"])

print('Number of events', length)
# iterate through all bins
for i in np.arange(0,10):
    pT_mask = np.ones(length)
    pT_mask = np.logical_and(pT_mask, data["D0_PT"]>bins[0,i])
    pT_mask = np.logical_and(pT_mask, data["D0_PT"]<=bins[0,i+1])
    for j in np.arange(0,10):
        eta_mask = pT_mask
        eta_mask = np.logical_and(eta_mask, data["D0_ETA"]>bins[i+1,j])
        eta_mask = np.logical_and(eta_mask, data["D0_ETA"]<=bins[i+1,j+1])
        selected_data = data[eta_mask]
        nevents = np.append(nevents, len(selected_data["D0_PT"]))
        # Write out bin
        if args.meson=="both":
            out_file_name = f"{args.path}/local/both/{args.polarity}_{args.year}_{args.size}_bin{j}{i}.root"
        else:
            out_file_name = f"{args.path}/local/{args.meson}_{args.polarity}_{args.year}_{args.size}_bin{j}{i}.root"
        out_tree = "D02Kpi_Tuple/DecayTree"
        print(f"Writing to {out_file_name}...")
        out_file = uproot.recreate(out_file_name)
        branches = {column: ak.type(selected_data[column]) for column in selected_data.fields}
        out_file.mktree(out_tree, branches)
        out_file[out_tree].extend({branch: selected_data[branch] for branch in branches.keys()})
        out_file.close()

#Creating root files for pT
for i in np.arange(0, 10):
    pT_mask = np.ones(length)
    pT_mask = np.logical_and(pT_mask, data["D0_PT"]>bins_pT[i])
    pT_mask = np.logical_and(pT_mask, data["D0_PT"]<=bins_pT[i+1])
    selected_data = data[pT_mask]
    nevents_pT = np.append(nevents_pT, len(selected_data["D0_PT"]))

    if args.meson=="both":
        out_file_name = f"{args.path}/pT/both/{args.polarity}_{args.year}_{args.size}_bin{i}.root"
    else:
        out_file_name = f"{args.path}/pT/{args.meson}_{args.polarity}_{args.year}_{args.size}_bin{i}.root"
    out_tree = "D02Kpi_Tuple/DecayTree"
    print(f"Writing to {out_file_name}...")
    out_file = uproot.recreate(out_file_name)
    branches = {column: ak.type(selected_data[column]) for column in selected_data.fields}
    out_file.mktree(out_tree, branches)
    out_file[out_tree].extend({branch: selected_data[branch] for branch in branches.keys()})
    out_file.close()

# # Create root files for eta 
for i in np.arange(0, 10):
    eta_mask = np.ones(length)
    eta_mask = np.logical_and(eta_mask, data["D0_ETA"]>bins_eta[i])
    eta_mask = np.logical_and(eta_mask, data["D0_ETA"]<=bins_eta[i+1])
    selected_data = data[eta_mask]
    nevents_eta = np.append(nevents_eta, len(selected_data["D0_ETA"]))

    if args.meson=="both":
        out_file_name = f"{args.path}/eta/both/{args.polarity}_{args.year}_{args.size}_bin{i}.root"
    else:
        out_file_name = f"{args.path}/eta/{args.meson}_{args.polarity}_{args.year}_{args.size}_bin{i}.root"
    out_tree = "D02Kpi_Tuple/DecayTree"
    print(f"Writing to {out_file_name}...")
    out_file = uproot.recreate(out_file_name)
    branches = {column: ak.type(selected_data[column]) for column in selected_data.fields}
    out_file.mktree(out_tree, branches)
    out_file[out_tree].extend({branch: selected_data[branch] for branch in branches.keys()})
    out_file.close()

# write out number of events in each bin
np.savetxt(f"{args.bin_path}/number_of_events_{args.meson}_{args.polarity}_{args.year}_{args.size}.txt", nevents, delimiter=',')
np.savetxt(f"{args.bin_path}/number_of_events_pT_{args.meson}_{args.polarity}_{args.year}_{args.size}.txt", nevents_pT, delimiter=',')
np.savetxt(f"{args.bin_path}/number_of_events_eta_{args.meson}_{args.polarity}_{args.year}_{args.size}.txt", nevents_eta, delimiter=',')