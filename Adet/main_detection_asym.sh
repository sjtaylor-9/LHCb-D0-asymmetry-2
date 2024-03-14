# Make the necessary directories in the individual pT, eta bins
for scheme in 'eta' 'pT'; do
    for year in '16' '17' '18'; do
        for polarity in 'up' 'down'; do
            mkdir Outputs/${scheme}
            mkdir Outputs/${scheme}/${year}
            mkdir Outputs/${scheme}/${year}/${polarity}
        done
    done
done



#python Adet/detection_asym.py --year $year --polarity $polarity --scheme $scheme --input '/eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet' --bin '6' --path 'Adet/Outputs/'$scheme/$year