import pandas as pd
import os
import argparse

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    --path      Used to specify the directory in which the root files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --max_file  Used to obtain the amount files to combine
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
        "--max_file",
        type=int,
        required=True,
        help="flag to set the max file"
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

args = parse_arguments() 

integer_list = []
for i in range(1, args.max_file+1, +1):
    integer_list.append(i)


# Initialize combined dataframe
combined_chunks = []
chunk_counter = 0  # Counter for tracking the number of chunks processed

for int_value in integer_list:
    print(f"Adding file {int_value}")
    try:
        # Read the CSV file for the current integer in chunks
        chunks = pd.read_csv(f"{args.path}/Pythia/Pythia_Data/simulated_data/pythia_hadronisation{int_value}.csv", chunksize=1000)

        # Process each chunk
        processed_chunks = []
        for chunk in chunks:
            # Add 2 million to the event numbers in the current chunk
            chunk['Event '] += 2000000*(int_value-1)
            processed_chunks.append(chunk)

        # Store processed chunks for concatenation
        combined_chunks.extend(processed_chunks)
        chunk_counter += 1

        # Check if it's time to save the combined data to a new file
        if chunk_counter % 30 == 0:
            # Concatenate processed chunks into final dataframe
            combined = pd.concat(combined_chunks, ignore_index=True)

            # Write the combined data to a new CSV file
            output_path = f"{args.path}/Pythia/Pythia_Data/combined_simulated_data_{chunk_counter}.csv"
            print("Saving file to:", output_path)
            combined.to_csv(output_path, index=False)

            # Reset combined_chunks for the next batch of chunks
            combined_chunks = []

    except FileNotFoundError:
        print(f"File pythia_hadronisation{int_value}.csv not found.")

# Concatenate any remaining processed chunks into final dataframe
if combined_chunks:
    combined = pd.concat(combined_chunks, ignore_index=True)
    # Write the remaining combined data to a new CSV file
    output_path = f"{args.path}/Pythia/Pythia_Data/combined_simulated_data_{chunk_counter}.csv"
    print("Saving file to:", output_path)
    combined.to_csv(output_path, index=True)
    print("File Saved")