directory=$1 # Need to know where to save binned data using the binning scheme found in that directory
year=$2 # Need to know the year to access binning scheme
size=$3 # Need to know size used to create the binning scheme
seeding=$4
min_seeding=$5 # Does a loop for size in the range size to minseeding

# Need to configure pythia folder after downloading it. Then copy pythia_hadronisation into that file
mkdir "Pythia/Pythia_Data"
mkdir "Pythia/Pythia_Data/simulated_data"
mkdir $directory"/Pythia"
mkdir $directory"/Pythia/binned_data"
mkdir $directory"/Pythia/binned_data/local"
mkdir $directory"/Pythia/binned_data/pT"
mkdir $directory"/Pythia/binned_data/eta"
mkdir $directory"/Pythia/binned_data/binning_scheme"
mkdir $directory"/Pythia//results"
mkdir $directory"/Pythia/asymmetry"
mkdir $directory"/Pythia/asymmetry/local"
mkdir $directory"/Pythia/asymmetry/pT"
mkdir $directory"/Pythia/asymmetry/eta"

echo "The necessary directories have been created"

# if ! [[ "$min_seeding" =~ ^[0-9]+$ ]]; then
#   echo "WARNING: You did not select a valid option for the minsize fit"
#   echo
#   echo "The simulation will run over sizes in the array [1, ..., $min_seeding]"
#   min_seeding=1
# else
#   echo "The selection will run over sizes in the array [$seeding,...,$min_seeding]"
# fi

# echo "Starting Simulation"
# while [ $seeding -ge $min_seeding ]; do
#     echo "Inside the loop. Seeding #: $seeding"
#     Pythia/pythia_hadronisation/runpythia $seeding
#     git add ./Pythia/pythia_hadronisation/simulated_data/pythia_hadronisation$seeding.csv 
#     git commit -m "Commited $seeding file"
#     git push origin main
#     seeding=$((seeding - 1))  # For example, decrease 'size' by 1 in each iteration
# done
# echo "Finished simulation"

# python Pythia/pythia_hadronisation/combining_csv.py --path '/afs/cern.ch/work/l/lseelan/LHCb-D0-asymmetry-Semeter2' --max_file 85

# echo
# echo "Created a csv file with all the data"
# echo

# for integer in 30 60 85; do
#   echo "loop "$integer
#    python Pythia/pythia_hadronisation/multiple_candidates_pythia.py \
#     --input "Pythia/Pythia_Data/combined_simulated_data_"$integer".csv"  \
#     --path "Pythia/Pythia_Data" \
#     --file_suffix $integer
# done

echo
echo "The phase space selection criteria has been applied and multiple candidates have been removed from the simulation data"
echo

for meson in D0 D0bar; do 
  python Pythia/pythia_hadronisation/apply_binning_scheme_pythia.py \
    --year $year \
    --size $size \
    --meson $meson \
    --path $directory"/Pythia/binned_data" \
    --input "Pythia/Pythia_Data" \
    --bin_path $directory"/binned_data/binning_scheme"
done

# echo
# echo "The simulation data has been binned"
# echo

for scheme in eta pT; do
    python Pythia/pythia_hadronisation/Aprod_pythia.py \
      --path $directory"/Pythia/asymmetry/$scheme" \
      --input_bins $directory"/Pythia/binned_data/binning_scheme" \
      --scheme $scheme \
      --results_path $directory"/Pythia/results" \
      --input_global "Pythia/Pythia_Data"
done

echo "The global and local production asymmetries have been calculated"
echo