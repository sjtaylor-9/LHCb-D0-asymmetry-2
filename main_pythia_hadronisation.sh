directory=$1 # Directory in eos
year=$2
size=$3
minsize=$4 # Does a loop for size in the range size to minsize

mkdir $directory"/pythia_hadronisation"
mkdir $directory"/pythia_hadronisation/selected_data"
mkdir $directory"/pythia_hadronisation/binned_data"
mkdir $directory"/pythia_hadronisation/binned_data/binning_scheme"
mkdir $directory"/pythia_hadronisation/binned_data/rapidity"
mkdir $directory"/pythia_hadronisation/binned_data/pT"
mkdir $directory"/pythia_hadronisation/binned_data/local"
mkdir $directory"/pythia_hadronisation/results"
mkdir $directory"/pythia_hadronisation/results/local"
mkdir $directory"/pythia_hadronisation/results/pT"
mkdir $directory"/pythia_hadronisation/results/rapidity"

echo "The necessary directories have been created"
echo

python pythia_hadronisation/multiple_candidates_pythia.py --input $directory"/pythia_hadronisation"  --path $directory"/pythia_hadronisation/selected_data"

echo "The phase space selection criteria has been applied and multiple candidates have been removed from the simulation data"

python pythia_hadronisation/apply_binning_scheme_pythia.py --year $year --size $size --path $directory"/pythia_hadronisation/binned_data" --input $directory"/pythia_hadronisation/selected_data" --bin_path $directory"/pythia_hadronisation/binned_data/binning_scheme"

echo "The simulation data has been binned"

for scheme in local pT rapidity 
do
    python pythia_hadronisation/Aprod_pythia.py --path $directory"/pythia_hadronisation/results/$scheme" --input $directory"/pythia_hadronisation/binned_data/$scheme" --scheme $scheme
done

echo "The global and local production asymmetries have been calculated"
echo