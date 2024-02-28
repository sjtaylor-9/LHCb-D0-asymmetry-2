print("importing...")
import os
import string, random
import numpy as np
import gc
import ROOT as R
R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.ERROR) # mute warnings
from pathlib import Path
import argparse
from utils import plot


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
    "--control-mode",
    required=True,
    type=str,
    choices=["kpipi", "kspi"]
)
parser.add_argument(
    "--set-parameters-to",
    type=str,
    default="",
    help="Pass the fit_params_{args.control_mode}.txt location if you want to fix parameters to the result of a previous fit"
)
parser.add_argument(
    "--polarity",
    type=str,
    default="",
    help="The magnet polarity",
    choices=["magup", "magdown"],
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

# functions taken from: https://gitlab.cern.ch/lhcb-slb/AsymmetryTools/-/blob/master/projects/KPiTwoTrackAsymmetry/KPiTwoTrackAsymmetry/utilities.py?ref_type=heads
weight_truncation_value = 80.
weight_cut_value = 500. * weight_truncation_value

def id_generator(size=6, chars=string.ascii_uppercase):
    return ''.join(random.choice(chars) for _ in range(size))

def compute_integral_of_expression_in_tree ( input_tree, tformula_expression, cuts="(1)" ):
    """
    Quickly computes the integral of an expression in a ROOT.TTree by 
    exploiting the mean of the histogram (which is always exact and
    not influenced by the binning), together with the quick access to the
    number of entries. This is tested to be equal to a loop over entries and
    an explicit sum per entry. 
    
    @return double
    """
    
    histogram_name = id_generator(10)
    
    input_tree.Draw(f"{tformula_expression}>>{histogram_name}", cuts, "GOFF")
    histogram = R.gDirectory.Get(histogram_name);
    
    return histogram.GetMean()*histogram.GetEntries()

def normalise_weights( input_tree, weight_variable ):
    normalisation_factor = compute_integral_of_expression_in_tree( input_tree, weight_variable)
    
    normalisation = compute_integral_of_expression_in_tree( input_tree, f"({weight_variable}*{weight_variable})")
    
    if normalisation > 0:
        normalisation_factor /= normalisation
    else:
        print("Got invalid result for weight*weight")
        print(normalisation)
        exit()
    
    return normalisation_factor


print("reading data...")
ttree = R.TChain("weightedTree")
for file in files:
    ttree.Add(file)

plot_type = args.polarity
plot_type_capitalised = plot_type[0].upper() + plot_type[1:3] + plot_type[3].upper() + plot_type[4:]

if args.control_mode == "kpipi":
    mass_variable = "D_DTF_M_PV"
    mass_range = (1819, 1919)
    m_obs = R.RooRealVar(mass_variable, "m(K#pi#pi) [MeVc^{-2}]", *mass_range)
    plot_type_p = plot_type_capitalised + " #it{D^{+} #rightarrow K^{-}#pi^{+}#pi^{+}}"
    plot_type_n = plot_type_capitalised + " #it{D^{-} #rightarrow K^{+}#pi^{-}#pi^{-}}"
else:
    mass_variable = "D_DTF_M_PV_CON"
    mass_range = (1800, 1935)
    m_obs = R.RooRealVar(mass_variable, "m(K_{S}^{0}#pi) [MeVc^{-2}]", *mass_range)
    plot_type_p = plot_type_capitalised + " #it{D^{+} #rightarrow K^{0}_{S}#pi^{+}}"
    plot_type_n = plot_type_capitalised + " #it{D^{-} #rightarrow K^{0}_{S}#pi^{-}}"
    
ttree.SetBranchStatus("*", 0)
ttree.SetBranchStatus(mass_variable, 1)
ttree.SetBranchStatus("D_ID", 1)
ttree.SetBranchStatus("weight", 1)

norm_factor = normalise_weights(ttree, "weight")

numbins = (mass_range[-1] - mass_range[0])*1
bins_for_plotting = np.linspace(*mass_range, numbins+1)
hist_P = R.TH1F("hist_P", "hist_P", numbins, *mass_range)
hist_M = R.TH1F("hist_M", "hist_M", numbins, *mass_range)

ttree.Draw(f"{mass_variable}>>hist_P", f"TMath::Min({weight_truncation_value}, {norm_factor}*weight) * (D_ID>0 && {norm_factor}*weight<{weight_cut_value})", "goff")
ttree.Draw(f"{mass_variable}>>hist_M", f"TMath::Min({weight_truncation_value}, {norm_factor}*weight) * (D_ID<0 && {norm_factor}*weight<{weight_cut_value})", "goff")
roo_hist_P = R.RooDataHist("roo_hist_P", "Data", R.RooArgList(m_obs), hist_P)
roo_hist_M = R.RooDataHist("roo_hist_M", "Data", R.RooArgList(m_obs), hist_M)


def MakeKPiPiModel(model_dict: dict, suffix: str="") -> tuple:
    # create sub-dictionaries to store PDFs
    if "Signals" not in model_dict.keys():
        model_dict["Signals"] = {}
    if "Backgrounds" not in model_dict.keys():
        model_dict["Backgrounds"] = {}
    # create sub-dictionary to store parameters
    if "Parameters" not in model_dict["Signals"].keys():
        model_dict["Signals"]["Parameters"] = {}
    if "Parameters" not in model_dict["Backgrounds"].keys():
        model_dict["Backgrounds"]["Parameters"] = {}

    model_dict["Signals"]["Parameters"][f"Gmu{suffix}"] = R.RooRealVar(f"Gmu{suffix}", f"Gmu{suffix}", 1870, *mass_range)
    model_dict["Signals"]["Parameters"][f"Gsig{suffix}"] = R.RooRealVar(f"Gsig{suffix}", f"Gsig{suffix}", 6.6, 5, 10)
    model_dict["Signals"][f"Gaussian{suffix}"] = R.RooGaussian(f"Gaussian{suffix}", "Gaussian", m_obs, model_dict["Signals"]["Parameters"][f"Gmu{suffix}"], model_dict["Signals"]["Parameters"][f"Gsig{suffix}"])

    # model_dict["Signals"]["Parameters"][f"Jmu{suffix}"] = R.RooRealVar(f"Jmu{suffix}", f"Jmu{suffix}", 1870, *mass_range)
    model_dict["Signals"]["Parameters"][f"Jlam{suffix}"] = R.RooRealVar(f"Jlam{suffix}", f"Jlam{suffix}", 14.2, 5, 20)
    model_dict["Signals"]["Parameters"][f"Jgam{suffix}"] = R.RooRealVar(f"Jgam{suffix}", f"Jgam{suffix}", 0.12, -1, 1)
    model_dict["Signals"]["Parameters"][f"Jdel{suffix}"] = R.RooRealVar(f"Jdel{suffix}", f"Jdel{suffix}", 1.32, 0, 1.7)
    model_dict["Signals"][f"Johnson{suffix}"] = R.RooJohnson(f"Johnson{suffix}", "Johnson", m_obs, model_dict["Signals"]["Parameters"][f"Gmu{suffix}"], model_dict["Signals"]["Parameters"][f"Jlam{suffix}"], model_dict["Signals"]["Parameters"][f"Jgam{suffix}"], model_dict["Signals"]["Parameters"][f"Jdel{suffix}"])

    model_dict["Signals"]["Parameters"][f"frac{suffix}"] = R.RooRealVar(f"frac{suffix}", f"frac{suffix}", 0.45, 0.35, 0.55)
    # model_dict["Signals"]["Parameters"][f"frac2{suffix}"] = R.RooRealVar(f"frac2{suffix}", f"frac2{suffix}", 0.3, 0, 1)
    model_dict["Signals"][f"Signal{suffix}"] = R.RooAddPdf(f"Signal{suffix}", "Signal", R.RooArgList(model_dict["Signals"][f"Gaussian{suffix}"], model_dict["Signals"][f"Johnson{suffix}"]), R.RooArgList(model_dict["Signals"]["Parameters"][f"frac{suffix}"]))

    model_dict["Backgrounds"]["Parameters"][f"Ec{suffix}"] = R.RooRealVar(f"Ec{suffix}", f"Ec{suffix}", -0.001, -0.01, -0.0001)
    model_dict["Backgrounds"][f"Exponential{suffix}"] = R.RooExponential(f"Exponential{suffix}", "Comb. background", m_obs, model_dict["Backgrounds"]["Parameters"][f"Ec{suffix}"])

    return model_dict["Signals"][f"Signal{suffix}"], model_dict["Backgrounds"][f"Exponential{suffix}"]

def MakeKSPiModel(model_dict: dict, suffix: str="") -> tuple:
    # create sub-dictionaries to store PDFs
    if "Signals" not in model_dict.keys():
        model_dict["Signals"] = {}
    if "Backgrounds" not in model_dict.keys():
        model_dict["Backgrounds"] = {}
    # create sub-dictionary to store parameters
    if "Parameters" not in model_dict["Signals"].keys():
        model_dict["Signals"]["Parameters"] = {}
    if "Parameters" not in model_dict["Backgrounds"].keys():
        model_dict["Backgrounds"]["Parameters"] = {}

    model_dict["Signals"]["Parameters"][f"Gmu{suffix}"] = R.RooRealVar(f"Gmu{suffix}", f"Gmu{suffix}", 1871, *mass_range)
    model_dict["Signals"]["Parameters"][f"Gsig1{suffix}"] = R.RooRealVar(f"Gsig1{suffix}", f"Gsig1{suffix}", 5.5, 4, 7)
    model_dict["Signals"]["Parameters"][f"Gsig2{suffix}"] = R.RooRealVar(f"Gsig2{suffix}", f"Gsig2{suffix}", 5.5, 4, 7)
    model_dict["Signals"][f"BifurGauss{suffix}"] = R.RooBifurGauss(f"BifurGauss{suffix}", "BifurGauss", m_obs, model_dict["Signals"]["Parameters"][f"Gmu{suffix}"], model_dict["Signals"]["Parameters"][f"Gsig1{suffix}"], model_dict["Signals"]["Parameters"][f"Gsig2{suffix}"])
    
    # model_dict["Signals"]["Parameters"][f"CBmu{suffix}"] = R.RooRealVar(f"CBmu{suffix}", f"CBmu{suffix}", 1871, *mass_range)
    model_dict["Signals"]["Parameters"][f"CBsig{suffix}"] = R.RooRealVar(f"CBsig{suffix}", f"CBsig{suffix}", 10, 5, 20)
    model_dict["Signals"]["Parameters"][f"CBa{suffix}"] = R.RooRealVar(f"CBa{suffix}", f"CBa{suffix}", 1.5, 0, 10)
    model_dict["Signals"]["Parameters"][f"CBn{suffix}"] = R.RooRealVar(f"CBn{suffix}", f"CBn{suffix}", 20, 0, 30)
    model_dict["Signals"][f"CrystalBall{suffix}"] = R.RooCrystalBall(f"CrystalBall{suffix}", "CrystalBall", m_obs, model_dict["Signals"]["Parameters"][f"Gmu{suffix}"], model_dict["Signals"]["Parameters"][f"CBsig{suffix}"], model_dict["Signals"]["Parameters"][f"CBa{suffix}"], model_dict["Signals"]["Parameters"][f"CBn{suffix}"])

    model_dict["Signals"]["Parameters"][f"frac{suffix}"] = R.RooRealVar(f"frac{suffix}", f"frac{suffix}", 0.22, 0, 0.55)
    # model_dict["Signals"]["Parameters"][f"frac2{suffix}"] = R.RooRealVar(f"frac2{suffix}", f"frac2{suffix}", 0.22, 0, 0.55)
    model_dict["Signals"][f"Signal{suffix}"] = R.RooAddPdf(f"Signal{suffix}", "Signal", R.RooArgList(model_dict["Signals"][f"BifurGauss{suffix}"], model_dict["Signals"][f"CrystalBall{suffix}"]), R.RooArgList(model_dict["Signals"]["Parameters"][f"frac{suffix}"]))

    model_dict["Backgrounds"]["Parameters"][f"Ec{suffix}"] = R.RooRealVar(f"Ec{suffix}", f"Ec{suffix}", -0.001, -0.01, -0.0001)
    model_dict["Backgrounds"][f"Exponential{suffix}"] = R.RooExponential(f"Exponential{suffix}", "Comb. background", m_obs, model_dict["Backgrounds"]["Parameters"][f"Ec{suffix}"])

    return model_dict["Signals"][f"Signal{suffix}"], model_dict["Backgrounds"][f"Exponential{suffix}"]

def ReadParameters(model_dict: dict, parameter_path: str):
    with open(parameter_path) as f:
        for l in f:
            line = l.strip()
            if line == "":
                continue

            colon = line.find(":")
            name = line[:colon]

            pm = line.find("+/-")
            val = float(line[colon+1:pm])
            err = float(line[pm+3:])

            if name in model_dict["Backgrounds"]["Parameters"].keys():
                model_dict["Backgrounds"]["Parameters"][name].setVal(val)
                model_dict["Backgrounds"]["Parameters"][name].setConstant(True)
            else:
                model_dict["Signals"]["Parameters"][name].setVal(val)
                model_dict["Signals"]["Parameters"][name].setConstant(True)


print("building model...")
model_dict = {}
Signal_P, Background_P = MakeKPiPiModel(model_dict, "_P") if args.control_mode == "kpipi" else MakeKSPiModel(model_dict, "_P")
Signal_M, Background_M = MakeKPiPiModel(model_dict, "_M") if args.control_mode == "kpipi" else MakeKSPiModel(model_dict, "_M")
if args.set_parameters_to != "":
    ReadParameters(model_dict, args.set_parameters_to)

######################
## simultaneous fit ##
######################
category = R.RooCategory("category", "category")
simultaneous_model = R.RooSimultaneous("simultaneous", "simultaneous", category)

N_sig = R.RooRealVar("N_sig", "N_sig", 0.9*(hist_P.GetEntries()+hist_M.GetEntries()), 0, (hist_P.GetEntries()+hist_M.GetEntries()))
N_bkg = R.RooRealVar("N_bkg", "N_sig", 0.1*(hist_P.GetEntries()+hist_M.GetEntries()), 0, (hist_P.GetEntries()+hist_M.GetEntries()))

A_sig = R.RooRealVar("A_sig", "A_sig", -1, -7.5, 7.5) # in %
A_bkg = R.RooRealVar("A_bkg", "A_bkg", 0.00, -7.5, 7.5) # in %

Nsig_P = R.RooFormulaVar("Nsig_P", "Nsig_P", "0.5*(1+@0/100.0)*@1", R.RooArgList(A_sig, N_sig))
Nsig_M = R.RooFormulaVar("Nsig_M", "Nsig_M", "0.5*(1-@0/100.0)*@1", R.RooArgList(A_sig, N_sig))
Nbkg_P = R.RooFormulaVar("Nbkg_P", "Nbkg_P", "0.5*(1+@0/100.0)*@1", R.RooArgList(A_bkg, N_bkg))
Nbkg_M = R.RooFormulaVar("Nbkg_M", "Nbkg_M", "0.5*(1-@0/100.0)*@1", R.RooArgList(A_bkg, N_bkg))

# D+
category.defineType("P")
Model_P = {
    "total": R.RooAddPdf("Model_P", "D^{+} Model", R.RooArgList(Signal_P, Background_P), R.RooArgList(Nsig_P, Nbkg_P)),
    "signals": {var.GetName(): var.GetTitle() for key, var in model_dict["Signals"].items() if key[-2:]=="_P"},
    "backgrounds": {var.GetName(): var.GetTitle() for key, var in model_dict["Backgrounds"].items() if key[-2:]=="_P"},
}
simultaneous_model.addPdf(Model_P["total"], "P")

# D-
category.defineType("M")
Model_M = {
    "total": R.RooAddPdf("Model_M", "D^{-} Model", R.RooArgList(Signal_M, Background_M), R.RooArgList(Nsig_M, Nbkg_M)),
    "signals": {var.GetName(): var.GetTitle() for key, var in model_dict["Signals"].items() if key[-2:]=="_M"},
    "backgrounds": {var.GetName(): var.GetTitle() for key, var in model_dict["Backgrounds"].items() if key[-2:]=="_M"},
}
simultaneous_model.addPdf(Model_M["total"], "M")

imports = [R.RooFit.Import("P", roo_hist_P), R.RooFit.Import("M", roo_hist_M)]
roo_hist = R.RooDataHist("roo_hist", "roo_hist", R.RooArgList(m_obs),R.RooFit.Index(category), *imports)



########################
## maximum likelihood ##
########################
simultaneous_model.fitTo(roo_hist, R.RooFit.Save(True), R.RooFit.Extended(True), R.RooFit.SumW2Error(True))

plot(m_obs, roo_hist_P, Model_P, plot_type_p, bins=bins_for_plotting, setlogy=True, save_to=f"{output_directory}/positive_likelihood_fit", unit="MeVc^{-2}")
plot(m_obs, roo_hist_M, Model_M, plot_type_n, bins=bins_for_plotting, setlogy=True, save_to=f"{output_directory}/negative_likelihood_fit", unit="MeVc^{-2}")


####################
## FINAL CHI2 FIT ##
####################
nll = R.RooChi2Var("nll", "nll", simultaneous_model, roo_hist, R.RooFit.Verbose(False),  R.RooFit.Extended(True), R.RooFit.PrintLevel(8), R.RooFit.NumCPU(1), R.RooFit.DataError(R.RooAbsData.SumW2))
m = R.RooMinuit(nll)
m.setStrategy(1)
m.setVerbose(False)#True)
m.setPrintLevel(-1)
m.migrad()
m.hesse()

plot(m_obs, roo_hist_P, Model_P, plot_type_p, bins=bins_for_plotting, setlogy=True, save_to=f"{output_directory}/positive_chi2_fit", unit="MeVc^{-2}")
plot(m_obs, roo_hist_M, Model_M, plot_type_n, bins=bins_for_plotting, setlogy=True, save_to=f"{output_directory}/negative_chi2_fit", unit="MeVc^{-2}")

print("writing result to file...")
with open(f"{output_directory}/output_{args.control_mode}.txt", "w") as f:
    f.write(f"{A_sig.getValV()} +/- {A_sig.getError()} % \n")

with open(f"{output_directory}/fit_params_{args.control_mode}.txt", "w") as f:
    for name, variable in model_dict["Signals"]["Parameters"].items():
        f.write(f"{name}: {variable.getValV():.5f} +/- {variable.getError():.5f} \n")
    for name, variable in model_dict["Backgrounds"]["Parameters"].items():
        f.write(f"{name}: {variable.getValV():.5f} +/- {variable.getError():.5f} \n")
