#!/bin/bash

# this file details what you would do at the command line, if you were to do everything locally on lxplus.
# instead write the instructions for the batch system to do it:



# first do the reweighting.
# you need to pass your input kpi data
# the kpipi and kspi data are read from a location on /eos
# so dont worry about that
python reweight_control_modes.py \
--year $1 --polarity $2 --n-estimators $3 --learning-rate $4 --max-depth $5 $6 \
--input /eos/lhcb/user/l/lseelan/Total/selected_data/${1}/${2:3:4}/both/:D02Kpi_Tuple/DecayTree \
--output-directory /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/Adet/$1/$2/${6:2:5}/
# this produces some new data files in the specified output directory
# of kpipi and kspi with the generated weights (data.root and temp_**.root)


# you can remove $3/$4/$5/ from everywhere if you stick with the 500/0.1/4/ configuration
# if you decide to investigate configurations then you need to include it so that directories arent overwritten


# next plot the results using the outputs from the previous command
python plot_reweighted_control_modes.py \
--input /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/Adet/$1/$2/${6:2:5}/data.root \
--kpipi-input /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/Adet/$1/$2/${6:2:5}/temp_kpipi.root \
--kspi-input /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/Adet/$1/$2/${6:2:5}/temp_kspi.root \
--output-directory /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/Adet/$1/$2/${6:2:5}/ --overwrite
# writes to the same location, doesn't actually overwrite anything
# as it only generates plots


# finally use the generated weights (if you deem them satisfactory)
# to reweight and fit the m(D+/-) distributions for both kpipi and kspi modes
for mode in {kpipi,kspi}; do
  python fit_control_modes.py \
  --input /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/Adet/$1/$2/${6:2:5}/temp_${mode}.root \
  --control-mode ${mode} \
  --output-directory /afs/cern.ch/user/s/sjtaylor/WorkSpace/D0_production_asymmetry_Sem2/LHCb_D0_asymmetry_2/Adet/$1/$2/${6:2:5}/${mode}/
done




# to submit jobs do:

#    $ conda activate my_env
#    (my_env) $ condor_submit Adet_from_condor.sub

# where my_env contains all necessary packages