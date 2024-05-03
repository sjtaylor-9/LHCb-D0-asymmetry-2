directory=$1
year=$2
size=$3

mkdir $directory"/binned_data"
mkdir $directory"/binned_data/binning_scheme"
mkdir $directory"/binned_data"
mkdir $directory"/binned_data/local"
mkdir $directory"/binned_data/pT"
mkdir $directory"/binned_data/eta"

These directories are made in the eos, which were used for Adet. They are not needed unless redoing Adet.
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/$year
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/$year/local
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/$year/pT
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/$year/eta
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/$year/local/both
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/$year/pT/both
mkdir /eos/lhcb/user/l/lseelan/Total/binned_data/$year/eta/both

# Size 10: lowest. Size 20: twice as big and includes data from Size 10

python create_binning_scheme.py \
    --year $year \
    --size $size \
    --path $directory"/binned_data/binning_scheme" \
    --input "/eos/lhcb/user/l/lseelan/Total/selected_data"

for meson in D0 D0bar; do 
    for polarity in up down; do    
        python apply_binning_scheme.py \
            --year $year \
            --size $size \
            --meson $meson \
            --polarity $polarity \
            --path $directory"/binned_data/" \
            --input "/eos/lhcb/user/l/lseelan/Total/selected_data" \
            --bin_path $directory"/binned_data/binning_scheme"
        python plot_phase_space.py \
            --year $year \
            --size $size \
            --meson $meson \
            --polarity $polarity \
            --path $directory"/binned_data/binning_scheme" \
            --input "/eos/lhcb/user/l/lseelan/Total/selected_data" \
            --bin_path $directory"/binned_data/binning_scheme"
        echo "Ploted 2D graph"
    done
done

echo "The data has been binned"
echo

    