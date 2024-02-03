print("importing...")
import os
import gc
import uproot as ur
import awkward as ak
import numpy as np
from hep_ml.reweight import GBReweighter
from pathlib import Path
import argparse


parser = argparse.ArgumentParser()
parser.add_argument(
    "--input",
    type=str,
    nargs="+",
    # required=True,
    help="Can be .root file or directory containing only .root files, or a list of"
)
parser.add_argument(
    "--output-directory",
    type=str,
    required=True,
    help="Location to write output to (need not exist yet)"
)
parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Pass to overwrite any existing data in --output-directory"
)
parser.add_argument(
    "--year",
    required=True,
    choices=["2016", "2017", "2018"]
)
parser.add_argument(
    "--polarity",
    required=True,
    choices=["magup", "magdown"]
)
parser.add_argument(
    "--n-estimators",
    required=True,
    type=int,
)
parser.add_argument(
    "--learning-rate",
    required=True,
    type=float,
)
parser.add_argument(
    "--max-depth",
    required=True,
    type=int,
)
parser.add_argument(
    "--small-kpipi",
    action="store_true",
    help="pass this to use the smaller kpipi file for testing configurations",
)
args = parser.parse_args()


files = {}
inputs = args.input
for input in inputs:
    unpack = input.split(":")
    file = unpack[0]
    tree = unpack[1] if len(unpack)>1 else None
    if os.path.isfile(file):
        files[file] = tree
    elif os.path.isdir(file):
        for f in os.listdir(file):
            files[f"{file}/{f}"] = tree
    else:
        raise ValueError(f"\"{file}\" is neither a file nor a directory")
output_directory = Path(args.output_directory)
output_directory.mkdir(parents=True, exist_ok=args.overwrite)


##############################################################################################################################
##############################################################################################################################
#################################################       CHANGE THIS         ##################################################
##############################################################################################################################
##############################################################################################################################
print("reading signal...")
variables = ["P1_ETA", "P2_ETA", "P1_PT", "P2_PT", "P1_PHI", "P2_PHI"]
kpi = ur.concatenate(files, variables)

kpi["kminus_pT"] = kpi["P1_PT"]
kpi["kminus_eta"] = kpi["P1_ETA"]
kpi["kminus_phi"] = kpi["P1_PHI"]
kpi["kminus_px"] = kpi["kminus_pT"] * np.cos(kpi["kminus_phi"])
kpi["kminus_py"] = kpi["kminus_pT"] * np.sin(kpi["kminus_phi"])
kpi["kminus_pz"] = kpi["kminus_pT"] / (np.tan(2 * np.arctan(np.exp(-(kpi["kminus_eta"])))))
kpi["kminus_p"] = np.sqrt(np.power(kpi["kminus_px"], 2) + np.power(kpi["kminus_py"], 2) + np.power(kpi["kminus_pz"], 2))
kpi["piplus_pT"] = kpi["P2_PT"]
kpi["piplus_eta"] = kpi["P2_ETA"]
kpi["piplus_phi"] = kpi["P2_PHI"]
kpi["piplus_px"] = kpi["piplus_pT"] * np.cos(kpi["piplus_phi"])
kpi["piplus_py"] = kpi["piplus_pT"] * np.sin(kpi["piplus_phi"])
kpi["piplus_pz"] = kpi["piplus_pT"] / (np.tan(2 * np.arctan(np.exp(-(kpi["piplus_eta"])))))




print(f"writing to {str(output_directory)}...")
outfile = ur.recreate(f"{str(output_directory)}/data.root")
branches = {column: ak.type(kpi[column]) for column in kpi.fields}
outfile.mktree(f"RS_DT", branches)
outfile[f"RS_DT"].extend({branch: kpi[branch] for branch in branches.keys()})
outfile.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


kpi = ak.to_dataframe(kpi)
gc.collect()

kpipi_to_kpi_variables = ["kminus_p", "kminus_eta", "kminus_phi", "piplus_pt", "piplus_eta", "piplus_phi"]
kspi_to_kpipi_variables = ["dplus_pt", "dplus_eta", "dplus_phi", "trigger_pi_pt", "trigger_pi_eta", "trigger_pi_phi"]
all_kpipi_variables = []
all_kpipi_variables.extend(kpipi_to_kpi_variables)
all_kpipi_variables.extend(kspi_to_kpipi_variables)



print("reading kpipi...")
vars = [
    "K_P", "K_ETA", "K_PHI", "Pi2_PT", "Pi2_ETA", "Pi2_PHI",
    "D_PT", "D_ETA", "D_PHI", "Pi1_PT", "Pi1_ETA", "Pi1_PHI",
    "D_ID", "D_DTF_M_PV"
]
kpipi_cut = "(K_PIDK>8.0) & (Pi2_PIDK<-5.0) & (D_PT/1000.0 < -4.5*(D_ETA - 3) + (15.2)) & (D_PT/1000.0 < 15.5*tanh(4.0*(D_ETA-1.8))) & (D_PT/1000.0 < 98 - 19.5*D_ETA)"
kpipi = ur.concatenate(f"/eos/lhcb/wg/TrackingAlignment/KPiAsymmetry/{args.year}/{args.polarity}/kpipi{'_small' if args.small_kpipi else ''}.root", vars, cut=kpipi_cut)
kpipi = ak.to_dataframe(kpipi)
gc.collect()

for i in range(len(all_kpipi_variables)):
    kpipi.columns.values[i] = all_kpipi_variables[i]
kpipi["kminus_theta"] = 2 * np.arctan(np.exp(-1*kpipi["kminus_eta"]))
kpipi["kminus_px"] = kpipi["kminus_p"] * np.sin(kpipi["kminus_theta"]) * np.cos(kpipi["kminus_phi"])
kpipi["kminus_py"] = kpipi["kminus_p"] * np.sin(kpipi["kminus_theta"]) * np.sin(kpipi["kminus_phi"])
kpipi["kminus_pz"] = kpipi["kminus_p"] * np.cos(kpipi["kminus_theta"])
kpipi["piplus_theta"] = 2 * np.arctan(np.exp(-1*kpipi["piplus_eta"]))
kpipi["piplus_px"] = kpipi["piplus_pt"] * np.cos(kpipi["piplus_phi"])
kpipi["piplus_py"] = kpipi["piplus_pt"] * np.sin(kpipi["piplus_phi"])
kpipi["piplus_pz"] = kpipi["piplus_pt"] / np.tan(kpipi["piplus_theta"])
kpipi["dplus_theta"] = 2 * np.arctan(np.exp(-1*kpipi["dplus_eta"]))
kpipi["dplus_px"] = kpipi["dplus_pt"] * np.cos(kpipi["dplus_phi"])
kpipi["dplus_py"] = kpipi["dplus_pt"] * np.sin(kpipi["dplus_phi"])
kpipi["dplus_pz"] = kpipi["dplus_pt"] / np.tan(kpipi["dplus_theta"])
kpipi["trigger_pi_theta"] = 2 * np.arctan(np.exp(-1*kpipi["trigger_pi_eta"]))
kpipi["trigger_pi_px"] = kpipi["trigger_pi_pt"] * np.cos(kpipi["trigger_pi_phi"])
kpipi["trigger_pi_py"] = kpipi["trigger_pi_pt"] * np.sin(kpipi["trigger_pi_phi"])
kpipi["trigger_pi_pz"] = kpipi["trigger_pi_pt"] / np.tan(kpipi["trigger_pi_theta"])



print("reading kspi...")
vars = [
    "D_PT", "D_ETA", "D_PHI", "H_PT", "H_ETA", "H_PHI",
    "D_ID", "D_DTF_M_PV_CON"
]
kspi_cut = "(D_PT/1000.0 < -4.5*(D_ETA - 3) + (15.2)) & (D_PT/1000.0 < 15.5*tanh(4.0*(D_ETA-1.8))) & (D_PT/1000.0 < 98 - 19.5*D_ETA)"
kspi = ur.concatenate(f"/eos/lhcb/wg/TrackingAlignment/KPiAsymmetry/{args.year}/{args.polarity}/kspi.root", vars, cut=kspi_cut)
kspi = ak.to_dataframe(kspi)
gc.collect()

for i in range(len(kspi_to_kpipi_variables)):
    kspi.columns.values[i] = kspi_to_kpipi_variables[i]
kspi["dplus_theta"] = 2 * np.arctan(np.exp(-1*kspi["dplus_eta"]))
kspi["dplus_px"] = kspi["dplus_pt"] * np.cos(kspi["dplus_phi"])
kspi["dplus_py"] = kspi["dplus_pt"] * np.sin(kspi["dplus_phi"])
kspi["dplus_pz"] = kspi["dplus_pt"] / np.tan(kspi["dplus_theta"])
kspi["trigger_pi_theta"] = 2 * np.arctan(np.exp(-1*kspi["trigger_pi_eta"]))
kspi["trigger_pi_px"] = kspi["trigger_pi_pt"] * np.cos(kspi["trigger_pi_phi"])
kspi["trigger_pi_py"] = kspi["trigger_pi_pt"] * np.sin(kspi["trigger_pi_phi"])
kspi["trigger_pi_pz"] = kspi["trigger_pi_pt"] / np.tan(kspi["trigger_pi_theta"])



def kpipi_reweighter(kpi, kpipi, original_weights, columns):
    print(f"reweighting by {columns}...")
    reweighter = GBReweighter(n_estimators=args.n_estimators, learning_rate=args.learning_rate, max_depth=args.max_depth, min_samples_leaf=1000, loss_regularization=5.0)
    reweighter.fit(original=kpipi[columns], target=kpi[columns], original_weight=original_weights)
    weights = reweighter.predict_weights(kpipi[columns], original_weight=original_weights)

    print(f"weighted, effective number of kpipi events: {np.sum(weights):.0f}, {np.power(np.sum(weights), 2) / np.sum(np.power(weights, 2)):.0f}")

    return weights

print(f"initial number of kpipi events: {kpipi.shape[0]:.0f}")
kpipi_weights = kpipi_reweighter(kpi, kpipi, np.ones(kpipi.shape[0]), ["kminus_p", "kminus_eta", "kminus_phi", "piplus_pt", "piplus_eta", "piplus_phi"])
gc.collect()

kpipi["weight"] = kpipi_weights
out_file_name = f"{str(output_directory)}/temp_kpipi.root"
out_tree_name = "weightedTree"
print(f"Writing to {out_file_name}...")
with ur.recreate(out_file_name) as out_file:
    branches = {}
    for column in kpipi.columns:
        branches[column] = np.dtype(kpipi[column])
    out_file.mktree(out_tree_name, branches)
    out_file[out_tree_name].extend({branch: kpipi[branch] for branch in branches.keys()})
gc.collect()



def kspi_reweighter(kpipi, kspi, original_weights, target_weights, columns):
    print(f"reweighting by {columns}...")
    reweighter = GBReweighter(n_estimators=args.n_estimators, learning_rate=args.learning_rate, max_depth=args.max_depth, min_samples_leaf=1000, loss_regularization=5.0)
    reweighter.fit(original=kspi[columns], target=kpipi[columns], original_weight=original_weights, target_weight=target_weights)
    weights = reweighter.predict_weights(kspi[columns], original_weight=original_weights)

    print(f"weighted, effective number of kspi events: {np.sum(weights):.0f}, {np.power(np.sum(weights), 2) / np.sum(np.power(weights, 2)):.0f}")

    return weights

print(f"intital number of kspi events: {kspi.shape[0]:.0f}")
kspi_weights = kspi_reweighter(kpipi, kspi, np.ones(kspi.shape[0]), kpipi_weights, ["dplus_pt", "dplus_eta", "dplus_phi", "trigger_pi_pt", "trigger_pi_eta", "trigger_pi_phi"])
gc.collect()

kspi["weight"] = kspi_weights
out_file_name = f"{str(output_directory)}/temp_kspi.root"
out_tree_name = "weightedTree"
print(f"Writing to {out_file_name}...")
with ur.recreate(out_file_name) as out_file:
    branches = {}
    for column in kspi.columns:
        branches[column] = np.dtype(kspi[column])
    out_file.mktree(out_tree_name, branches)
    out_file[out_tree_name].extend({branch: kspi[branch] for branch in branches.keys()})
