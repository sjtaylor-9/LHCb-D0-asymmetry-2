# Make the necessary directories in the individual pT, eta bins and
for scheme in 'eta' 'pT'; do
    for year in '16' '17' '18'; do
        for polarity in 'up' 'down'; do
            mkdir Outputs/${scheme}
            mkdir Outputs/${scheme}/${year}
            mkdir Outputs/${scheme}/${year}/${polarity}

            # Calculate the detection asymmetry in each of the bins and save the outputs to the directories
            for bin in {0..9}; do
                python detection_asym.py \
                --year ${year} \
                --polarity ${polarity} \
                --scheme ${scheme} \
                --input '/eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet' \
                --bin ${bin} \
                --path 'Adet/Outputs/'${scheme}/${year}/${polarity}
            done
        done
    done
done