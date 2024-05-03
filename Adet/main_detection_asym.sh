directory=$1
size=$2
year=$3

mkdir Adet/Outputs
mkdir Adet/Outputs/global

# Make the necessary directories to calculate the global detection asymmetries
for years in '16' '17' '18'; do
    for polarity in 'up' 'down'; do
        mkdir Adet/Outputs/global/${years}
        mkdir Adet/Outputs/global/${years}/${polarity}

        # Calculate the global detection asymmetry for each year and polarity
        python Adet/detection_asym.py \
        --year ${years} \
        --polarity ${polarity} \
        --scheme 'global' \
        --input '/eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet/global' \
        --path 'Adet/Outputs/global/'${years}/${polarity}
    done
done


# Make the necessary directories in the individual pT, eta bins and
for scheme in 'eta'; do
    for years in '16' '17' '18'; do
        for polarity in 'down' 'up'; do
            mkdir Adet/Outputs/${scheme}
            mkdir Adet/Outputs/${scheme}/${years}
            mkdir Adet/Outputs/${scheme}/${years}/${polarity}

            # Calculate the detection asymmetry in each of the bins and save the outputs to the directories
            for bin in {0..9}; do
                python Adet/detection_asym.py \
                --year ${years} \
                --polarity ${polarity} \
                --scheme ${scheme} \
                --input '/eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet' \
                --bin ${bin} \
                --path 'Adet/Outputs/'${scheme}/${years}/${polarity} 
            done
        done
    done
done


for scheme in 'pT' 'eta'; do
    python Adet/plot_Adet.py \
    --year ${year} \
    --scheme ${scheme} \
    --input 'Adet/Outputs' \
    --path 'Adet/Outputs/'${scheme}/${year} \
    --bin_path $directory"/binned_data/binning_scheme" \
    --size $size
done

