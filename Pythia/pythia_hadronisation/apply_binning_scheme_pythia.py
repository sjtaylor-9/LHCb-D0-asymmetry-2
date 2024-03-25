"""
apply_binning_scheme_pythia.py

This code applies the given binning scheme to a set of simulated data, and generates .csv files cotaining the events in each bin. There are 3 binning schemes used: pT, eta and (pT, eta).
The year of interest, size of the data, and meson to be analysed must be specified using the required flags --year --size --meson. There also are the flags --input --path and --bin_path, which are not required. These are used to specify the directory where the input data is located, where the binning scheme can be found and where the output file should be written, respectively. By default it is set to be the current working directory.
It outputs the .csv files with the events in each individual bin, as well as a .txt file with the number of events in each bin.

Authors: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)/
Last edited: 17th February 2024
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #
import os
import argparse
import numpy as np
import csv
import pandas as pd
import glob


# - - - - - - - - FUNCTIONS - - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size      Used to specify the amount of events the user is interested in analysing. The integers specify the number of root files to be read in.
    --meson     Used to specify the meson the user is interested in.
                The argument must be one of: [D0, D0bar].
    --path      Used to specify the directory in which the root files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --bin_path  Used to specify the directory in which the binning scheme should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --
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
        type=int,
        required=True,
       help="flag to set the size of the dataset that the binning scheme was created from."
    )
    parser.add_argument(
        "--meson",
        type=str,
        choices=["D0","D0bar","both"],
        required=False,
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
        
# - - - - - - - MAIN BODY - - - - - - - #
args = parse_arguments()

if args.meson == 'D0':
    files = glob.glob(f'{args.input}/D0_clean_pythia_data_*.csv')
elif args.meson == 'D0bar':
    files = glob.glob(f'{args.input}/D0bar_clean_pythia_data_*.csv')
elif args.meson == 'both':
    files = glob.glob(f'{args.input}/clean_pythia_data_*.csv')

# Define an empty list to store chunked dataframes
chunked_dfs = []

# Iterate over each matching file
for file in files:
    # Read CSV file in chunks
    chunks = pd.read_csv(file, usecols=[' PT ', ' Y '], chunksize=100000)  # Adjust chunksize as needed
    # Iterate over each chunk and store it in the list
    for chunk in chunks:
        chunked_dfs.append(chunk)

# Concatenate all chunked dataframes into a single DataFrame
data = pd.concat(chunked_dfs, ignore_index=True)

# Extract 'PT' and 'Y' columns from the DataFrame
pT = data[' PT ']
rapidity = data[' Y ']

# Loads in .txt of the binning scheme
bins = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_bins.txt", delimiter=',')
bins_pT = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_pT_bins.txt", delimiter=',')
bins_rapidity = np.loadtxt(f"{args.bin_path}/{args.year}_{args.size}_eta_bins.txt", delimiter=',')

nevents=np.empty(0)
nevents_pT =np.empty(0)
nevents_rapidity =np.empty(0)
length = len(data)
pT = 1000 * pT # Converts the transverse momentum outputted by Pythia from GeV/c to MeV/c as the bin edges are defined in terms MeV/c.

# iterate through all bins
for i in np.arange(0, 10):
    pT_mask = np.ones(length)
    # masks between the bins
    pT_mask = np.logical_and(pT_mask, pT > bins[0, i])
    pT_mask = np.logical_and(pT_mask, pT <= bins[0, i+1])
    for j in np.arange(0,10):
        rapidity_mask = pT_mask
        rapidity_mask = np.logical_and(rapidity_mask, rapidity > bins[i+1, j])
        rapidity_mask = np.logical_and(rapidity_mask, rapidity <= bins[i+1, j+1])
        selected_data = data[rapidity_mask]
        nevents = np.append(nevents, len(selected_data))
        # Write out bin
        out_file_name = f"{args.path}/local/{args.meson}_local_bin{j}{i}.csv"
        print(f"Writing to {out_file_name}...")
        with open(out_file_name, 'w', newline='') as file:
            writer = csv.writer(file)
            # Write data rows
            for row in selected_data:
                writer.writerow(row)
                
#Creating .csv files for pT
for i in np.arange(0, 10):
    pT_mask = np.ones(length)
    pT_mask = np.logical_and(pT_mask, pT > bins_pT[i])
    pT_mask = np.logical_and(pT_mask, pT <= bins_pT[i+1])
    selected_data = data[pT_mask]
    nevents_pT = np.append(nevents_pT, len(selected_data))
    
    out_file_name = f"{args.path}/pT/{args.meson}_pT_bin{i}.csv"
    print(f"Writing to {out_file_name}...")
    # Write out bin
    with open(out_file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        # Write data rows
        for row in selected_data:
            writer.writerow(row)

# Create .csv files for rapdidity
for i in np.arange(0, 10):
    rapidity_mask = np.ones(length)
    rapidity_mask = np.logical_and(rapidity_mask, rapidity > bins_rapidity[i])
    rapidity_mask = np.logical_and(rapidity_mask, rapidity <= bins_rapidity[i+1])
    selected_data = data[rapidity_mask]
    nevents_rapidity = np.append(nevents_rapidity, len(selected_data))
    
    out_file_name = f"{args.path}/eta/{args.meson}_eta_bin{i}.csv"
    print(f"Writing to {out_file_name}...")
    # Write out bin
    with open(out_file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        # Write data rows
        for row in selected_data:
            writer.writerow(row)

# write out number of events in each bin
np.savetxt(f"{args.path}/binning_scheme/number_of_events_local_{args.meson}.txt", nevents, delimiter=',')
np.savetxt(f"{args.path}/binning_scheme/number_of_events_pT_{args.meson}.txt", nevents_pT, delimiter=',')
np.savetxt(f"{args.path}/binning_scheme/number_of_events_eta_{args.meson}.txt", nevents_rapidity, delimiter=',')