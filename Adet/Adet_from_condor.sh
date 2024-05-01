#!/bin/bash

# this file details what you would do at the command line, if you were to do everything locally on lxplus.
# instead write the instructions for the batch system to do it:

# Get the bin number from $6
bin=$6
# Check if bin number is less than 10
# if [ $bin -lt 10 ]; then
#     # If bin number is less than 10, prepend '0' to make it two digits
#     bin="0$bin"
# fi
# echo "Formatted input: $bin"

# first do the reweighting.
# you need to pass your input kpi data
# the kpipi and kspi data are read from a location on /eos
# so dont worry about that
#--input /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/${1}/1mil_${1:2:2}_${2:3:4}.root:D02Kpi_Tuple/DecayTree \
python reweight_control_modes.py \
--year $1 --polarity $2 --n-estimators $3 --learning-rate $4 --max-depth $5 \
--input /eos/lhcb/user/l/lseelan/Total/binned_data/${1:2:2}/pT/both/${2:3}_${1:2:2}_70_bin${bin}.root:D02Kpi_Tuple/DecayTree \
--output-directory /eos/lhcb/user/l/lseelan/Total/Adet/pT_set_params/$1/$2/${bin}/ --overwrite

# this produces some new data files in the specified output directory
# of kpipi and kspi with the generated weights (data.root and temp_**.root)


# you can remove $3/$4/$5/ from everywhere if you stick with the 500/0.1/4/ configuration
# if you decide to investigate configurations then you need to include it so that directories arent overwritten


# next plot the results using the outputs from the previous command
python plot_reweighted_control_modes.py \
--input /eos/lhcb/user/l/lseelan/Total/Adet/pT_set_params/$1/$2/${bin}/data.root \
--kpipi-input /eos/lhcb/user/l/lseelan/Total/Adet/pT_set_params/$1/$2/${bin}/temp_kpipi.root \
--kspi-input /eos/lhcb/user/l/lseelan/Total/Adet/pT_set_params/$1/$2/${bin}/temp_kspi.root \
--output-directory /eos/lhcb/user/l/lseelan/Total/Adet/pT_set_params/$1/$2/${bin}/ --overwrite
# # writes to the same location, doesn't actually overwrite anything
# # as it only generates plots


# finally use the generated weights (if you deem them satisfactory)
# to reweight and fit the m(D+/-) distributions for both kpipi and kspi modes
for mode in {kpipi,kspi}; do
  python fit_control_modes.py \
  --input /eos/lhcb/user/l/lseelan/Total/Adet/pT_set_params/$1/$2/${bin}/temp_${mode}.root \
  --control-mode ${mode} \
  --output-directory /eos/lhcb/user/l/lseelan/Total/Adet/pT_set_params/$1/$2/${bin}/${mode}/ --overwrite \
  --polarity $2 \
  --set-parameters-to /eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet/global/$1/$2/${mode}/fit_params_${mode}.txt
done

# rm /eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet/pT/$1/$2/${bin}/data.root
# rm /eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet/pT/$1/$2/${bin}/temp_kpipi.root
# rm /eos/lhcb/user/s/sjtaylor/D0_asymmetry/Adet/pT/$1/$2/${bin}/temp_kspi.root


# to submit jobs do:

#    $ conda activate my_env
#    (my_env) $ condor_submit Adet_from_condor.sub

# where my_env contains all necessary packages