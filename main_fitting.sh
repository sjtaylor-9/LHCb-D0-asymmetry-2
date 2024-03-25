directory=$1
year=$2
size=$3
binned=$4
model=$5

mkdir $directory"/model_fitting"
mkdir $directory"/model_fitting/global"
mkdir $directory"/model_fitting/local"
mkdir $directory"/model_fitting/pT"
mkdir $directory"/model_fitting/eta"

for ind in {0..99}; do
    index=$( printf '%02d' $ind)
    mkdir $directory"/model_fitting/local/"$index
done

for ind in {0..9}; do
    index=$( printf '%01d' $ind)
    mkdir $directory"/model_fitting/pT/"$index
    mkdir $directory"/model_fitting/eta/"$index
done

# mkdir $directory"/raw_asymmetry_outcome"
# mkdir $directory"/raw_asymmetry_outcome/chi_squared"
# mkdir $directory"/raw_asymmetry_outcome/raw_asymmetry"
# mkdir $directory"/raw_asymmetry_outcome/raw_asymmetry/pT"
# mkdir $directory"/raw_asymmetry_outcome/raw_asymmetry/eta"
# mkdir $directory"/raw_asymmetry_outcome/raw_asymmetry/local"
# mkdir $directory"/results"

echo "The necessary directories have been created"
echo

# Size 10: lowest. Size 20: twice as big and includes data from Size 10

# Perform the global fit
echo "Fitting using Model "$model
python "Models/Model"$model"_pythonfiles/fit_global_model"$model".py"  \
    --year $year \
    --size $size \
    --path $directory"/model_fitting/global" \
    --binned_fit $binned \
    --input "/eos/lhcb/user/l/lseelan/Total/selected_data" \
    --scheme "total" \
    --initial_guess_path $directory"/model_fitting"
for meson in D0 D0bar; do
    for polarity in  down up; do 
            python "Models/Model"$model"_pythonfiles/model_fitting_model"$model".py" \
                --year $year \
                --size $size \
                --polarity $polarity \
                --meson $meson \
                --path $directory"/model_fitting/global" \
                --input "/eos/lhcb/user/l/lseelan/Total/selected_data" \
                --parameters_path $directory"/model_fitting/global" \
                --scheme 'total' \
                --binned_fit $binned
    done
done

echo "Plotted the global models"

# Perform the local fits in the (pT, eta) bins
for ind in {0..99}; do
    index=$( printf '%02d' $ind)
    python fit_global.py \
        --year $year \
        --size $size \
        --path $directory"/model_fitting/local/"$index \
        --binned_fit $binned \
        --input $directory"/binned_data/local" \
        --bin $index \
        --scheme 'pT_eta' \
        --initial_guess_path $directory"/model_fitting"
    echo "Fitted Bin "$index
done

for meson in D0 D0bar; do
   for polar in up down; do 
        for ind in {0..99}; do
            index=$( printf '%02d' $ind)
            python model_fitting.py \
                --year $year \
                --size $size \
                --meson $meson \
                --polarity $polar  \
                --path $directory"/model_fitting/local/"$index \
                --input $directory"/binned_data/local" \
                --parameters_path $directory"/model_fitting/local/"$index \
                --bin $index \
                --binned_fit $binned \
                --scheme 'pT_eta'
        done
    done
done

echo "Plotted the local models in the 100 (pT, eta) bins"

# Perform the local fits in the individual 10 pT and eta bins
for ind in {0..9}; do
    index=$( printf '%01d' $ind)
    python fit_global.py \
        --year $year \
        --size $size \
        --path $directory"/model_fitting/pT/"$index \
        --binned_fit $binned \
        --input $directory"/binned_data/pT" \
        --bin $index \
        --scheme 'pT'
    python fit_global.py \
        --year $year \
        --size $size \
        --path $directory"/model_fitting/eta/"$index \
        --binned_fit $binned \
        --input $directory"/binned_data/eta" \
        --bin $index \
        --scheme 'eta'
    echo "Fitted Bin "$index
done

for meson in D0 D0bar; do
   for polar in down up; do 
        for ind in {0..9}; do
            index=$( printf '%01d' $ind)
            python model_fitting.py \
                --year $year \
                --size $size \
                --meson $meson \
                --polarity $polar  \
                --path $directory"/model_fitting/pT/"$index \
                --input $directory"/binned_data/pT" \
                --parameters_path $directory"/model_fitting/pT/"$index \
                --bin $index \
                --binned_fit $binned \
                --scheme 'pT'
            python model_fitting.py \
                --year $year \
                --size $size --meson $meson \
                --polarity $polar  \
                --path $directory"/model_fitting/eta/"$index \
                --input $directory"/binned_data/eta" \
                --parameters_path $directory"/model_fitting/eta/"$index \
                --bin $index \
                --binned_fit $binned \
                --scheme 'eta'
        done
    done
done

echo "Plotted the local models in the 10 pT and 10 eta bins"
 