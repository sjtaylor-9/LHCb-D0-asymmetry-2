import os
import numpy as np
import argparse

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:

    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
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
    parser.add_argument(
        "--input",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input data should be found"
    )
    return parser.parse_args()

def read_data(input_path):
    """
    Reads the data from the specified input path.

    Args:
        input_path (str): Path to the input data file.

    Returns:
        data (array): Array containing the dataset.
    """
    data = np.genfromtxt(input_path, delimiter=',', skip_header=1)
    return data

def save_binning_scheme(pT_bins, eta_bins, output_path):
    """
    Saves the binning scheme to files.

    Args:
        pT_bins (array): Array containing the bin edges for pT.
        eta_bins (array): Array containing the bin edges for eta.
        output_path (str): Path to save the output files.
    """
    np.savetxt(f"{output_path}/pT_bins.txt", pT_bins, delimiter='\n')
    np.savetxt(f"{output_path}/eta_bins.txt", eta_bins, delimiter='\n')

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def create_binning_scheme(data, num_bins=10):
    """
    Creates a binning scheme for the dataset that splits it into 10 bins of pT and 10 bins of eta (rapidity)
    with approximately the same number of data points in each bin.

    Args:
        data (array): Array containing the dataset.
        num_bins (int): Number of bins for each variable.

    Returns:
        pT_bins (array): Array containing the bin edges for pT.
        eta_bins (array): Array containing the bin edges for eta.
    """
    # Extract pT and eta (rapidity) columns
    pT = data[:, 2]
    eta = data[:, 3]

    # Calculate the number of data points in each bin
    num_points_per_bin = len(data) // num_bins

    # Sort the data
    pT_sorted = np.sort(pT)
    eta_sorted = np.sort(eta)

    # Calculate bin edges for pT
    pT_bins = [pT_sorted[i * num_points_per_bin] for i in range(num_bins)]
    pT_bins.append(pT_sorted[-1])

    # Calculate bin edges for eta
    eta_bins = [eta_sorted[i * num_points_per_bin] for i in range(num_bins)]
    eta_bins.append(eta_sorted[-1])

    return pT_bins, eta_bins

def count_data_in_bins(data, pT_bins, eta_bins):
    """
    Counts the number of events in each individual bin for pT and eta separately.

    Args:
        data (array): Array containing the dataset.
        pT_bins (array): Array containing the bin edges for pT.
        eta_bins (array): Array containing the bin edges for eta.

    Returns:
        pT_counts (array): Array containing the counts of events in each individual pT bin.
        eta_counts (array): Array containing the counts of events in each individual eta bin.
    """
    # Extract pT and eta (rapidity) columns
    pT = data[:, 2]
    eta = data[:, 3]

    # Initialize counts arrays
    pT_counts = np.zeros(len(pT_bins) - 1, dtype=int)
    eta_counts = np.zeros(len(eta_bins) - 1, dtype=int)

    # Loop through pT bins
    for i in range(len(pT_bins) - 1):
        # Count events falling within current pT bin range
        pT_counts[i] = np.sum((pT >= pT_bins[i]) & (pT < pT_bins[i + 1]))

    # Loop through eta bins
    for i in range(len(eta_bins) - 1):
        # Count events falling within current eta bin range
        eta_counts[i] = np.sum((eta >= eta_bins[i]) & (eta < eta_bins[i + 1]))

    return pT_counts, eta_counts

# Parse arguments
args = parse_arguments()

# Read data
input_data = read_data(f"{args.input}/clean_pythia_data.csv")

# Create binning scheme
pT_bins, eta_bins = create_binning_scheme(input_data)

pT_counts, eta_counts = count_data_in_bins(input_data, pT_bins, eta_bins)

# Print the counts for pT bins
print("Counts for pT bins:")
for i, count in enumerate(pT_counts):
    print(f"Bin {i+1}: {count}")

# Print the counts for eta bins
print("\nCounts for eta bins:")
for i, count in enumerate(eta_counts):
    print(f"Bin {i+1}: {count}")

# Save binning scheme
save_binning_scheme(pT_bins, eta_bins, args.path)

print("Binning scheme saved successfully.")