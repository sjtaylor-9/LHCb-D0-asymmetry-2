print("importing...")
import os
import uproot as ur
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep; mplhep.style.use("LHCb2")
plt.style.use("seaborn-v0_8-muted")
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
    "--kpipi-input",
    type=str,
    nargs="+",
    required=True,
    help="Can be .root file or directory containing only .root files, or a list of"
)
parser.add_argument(
    "--kspi-input",
    type=str,
    nargs="+",
    required=True,
    help="Can be .root file or directory containing only .root files, or a list of"
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
kpipi_files = {}
kpipi_inputs = args.kpipi_input
for kpipi_input in kpipi_inputs:
    unpack = kpipi_input.split(":")
    file = unpack[0]
    tree = unpack[1] if len(unpack)>1 else None
    if os.path.isfile(file):
        kpipi_files[file] = tree
    elif os.path.isdir(file):
        for f in os.listdir(file):
            kpipi_files[f"{file}/{f}"] = tree
    else:
        raise ValueError(f"\"{file}\" is neither a file nor a directory")
kspi_files = {}
kspi_inputs = args.kspi_input
for kspi_input in kspi_inputs:
    unpack = kspi_input.split(":")
    file = unpack[0]
    tree = unpack[1] if len(unpack)>1 else None
    if os.path.isfile(file):
        kspi_files[file] = tree
    elif os.path.isdir(file):
        for f in os.listdir(file):
            kspi_files[f"{file}/{f}"] = tree
    else:
        raise ValueError(f"\"{file}\" is neither a file nor a directory")
output_directory = Path(args.output_directory)
output_directory.mkdir(parents=True, exist_ok=args.overwrite)


print("reading signal...")
kpi = ur.concatenate(files, library="pd")

print("reading kpipi...")
kpipi = ur.concatenate(kpipi_files, library="pd")

print("reading kspi...")
kspi = ur.concatenate(kspi_files, library="pd")

kw_args = {
   "bins": 50,
   "alpha": 0.7,
   "density": True
}

def kpipi_reweighter(kpi, kpipi, kw_args):
    
    print("plotting kpipi pt, eta, phi...")
    fig = plt.figure(figsize=(45,25))
    gs = fig.add_gridspec(5, 3, height_ratios=[4, 1, 1, 4, 1], hspace=0.1) # middle row is used as a spacer and so not used

    # kaon momentum
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ran = (0,1e5)
    before, *_ = ax1.hist(kpipi["kminus_p"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$", **kw_args);
    after, *_ = ax1.hist(kpipi["kminus_p"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax1.hist(kpi["kminus_p"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})K^{\mp}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax1.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax2.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.grid(True)
    ax2.set_ylim(0.8,1.2)
    ax2.set_xlabel(r"$p$ [MeV/$c$]")
    ax2.set_ylabel(r"$K\pi/K\pi\pi$")


    # kaon pseudorapidity
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1], sharex=ax3)
    ran = (1.8,5)
    before, *_ = ax3.hist(kpipi["kminus_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$", **kw_args);
    after, *_ = ax3.hist(kpipi["kminus_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax3.hist(kpi["kminus_eta"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})K^{\mp}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax3.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax4.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.grid(True)
    ax4.set_ylim(0.8,1.2)
    ax4.set_xlabel(r"$\eta$")
    ax4.set_ylabel(r"$K\pi/K\pi\pi$")


    # kaon azimuthal angle
    ax5 = fig.add_subplot(gs[0, 2])
    ax6 = fig.add_subplot(gs[1, 2], sharex=ax5)
    ran = (-np.pi,np.pi)
    before, *_ = ax5.hist(kpipi["kminus_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$", **kw_args);
    after, *_ = ax5.hist(kpipi["kminus_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax5.hist(kpi["kminus_phi"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})K^{\mp}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax5.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax6.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.grid(True)
    ax6.set_ylim(0.8,1.2)
    ax6.set_xlabel(r"$\phi$")
    ax6.set_ylabel(r"$K\pi/K\pi\pi$")


    # pion transverse momentum
    ax7 = fig.add_subplot(gs[3, 0])
    ax8 = fig.add_subplot(gs[4, 0], sharex=ax7)
    ran = (0,1e4)
    before, *_ = ax7.hist(kpipi["piplus_pt"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax7.hist(kpipi["piplus_pt"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax7.hist(kpi["piplus_pt"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})\pi^{\pm}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax7.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax8.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.grid(True)
    ax8.set_ylim(0.8,1.2)
    ax8.set_xlabel(r"$p_{T}$ [MeV/$c$]")
    ax8.set_ylabel(r"$K\pi/K\pi\pi$")


    # pion pseudorapidity
    ax9 = fig.add_subplot(gs[3, 1])
    ax10 = fig.add_subplot(gs[4, 1], sharex=ax9)
    ran = (1.8,5)
    before, *_ = ax9.hist(kpipi["piplus_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax9.hist(kpipi["piplus_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax9.hist(kpi["piplus_eta"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})\pi^{\pm}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax9.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax10.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.grid(True)
    ax10.set_ylim(0.8,1.2)
    ax10.set_xlabel(r"$\eta$")
    ax10.set_ylabel(r"$K\pi/K\pi\pi$")
    

    # pion azimuthal angle
    ax11 = fig.add_subplot(gs[3, 2])
    ax12 = fig.add_subplot(gs[4, 2], sharex=ax11)
    ran = (-np.pi,np.pi)
    before, *_ = ax11.hist(kpipi["piplus_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax11.hist(kpipi["piplus_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax11.hist(kpi["piplus_phi"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})\pi^{\pm}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax11.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax12.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.grid(True)
    ax12.set_ylim(0.8,1.2)
    ax12.set_xlabel(r"$\phi$")
    ax12.set_ylabel(r"$K\pi/K\pi\pi$")


    for a in [ax1,ax3,ax7,ax5,ax9,ax11]:
        a.set_ylabel("Normalised entries [a.u.]")
        a.legend(loc="best")

    plt.savefig(f"{str(output_directory)}/kpipi_ptetaphi_reweightings.pdf")
    
    print("plotting kpipi px, py, pz...")
    fig = plt.figure(figsize=(45,25))
    gs = fig.add_gridspec(5, 3, height_ratios=[4, 1, 1, 4, 1], hspace=0.1) # middle row is used as a spacer and so not used

    # kaon px
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ran = (0,1e4)
    before, *_ = ax1.hist(kpipi["kminus_px"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$", **kw_args);
    after, *_ = ax1.hist(kpipi["kminus_px"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax1.hist(kpi["kminus_px"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})K^{\mp}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax1.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax2.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.grid(True)
    ax2.set_ylim(0.8,1.2)
    ax2.set_xlabel(r"$p_{x}$ [MeV/$c$]")
    ax2.set_ylabel(r"$K\pi/K\pi\pi$")


    # kaon py
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1], sharex=ax3)
    ran = (0,1e4)
    before, *_ = ax3.hist(kpipi["kminus_py"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$", **kw_args);
    after, *_ = ax3.hist(kpipi["kminus_py"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax3.hist(kpi["kminus_py"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})K^{\mp}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax3.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax4.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.grid(True)
    ax4.set_ylim(0.8,1.2)
    ax4.set_xlabel(r"$p_{y}$ [MeV/$c$]")
    ax4.set_ylabel(r"$K\pi/K\pi\pi$")


    # kaon pz
    ax5 = fig.add_subplot(gs[0, 2])
    ax6 = fig.add_subplot(gs[1, 2], sharex=ax5)
    ran = (-np.pi,np.pi)
    ran = (0,1e5)
    before, *_ = ax5.hist(kpipi["kminus_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$", **kw_args);
    after, *_ = ax5.hist(kpipi["kminus_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})K^{\mp}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax5.hist(kpi["kminus_pz"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})K^{\mp}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax5.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax6.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.grid(True)
    ax6.set_ylim(0.8,1.2)
    ax6.set_xlabel(r"$p_{z}$ [MeV/$c$]")
    ax6.set_ylabel(r"$K\pi/K\pi\pi$")


    # pion px
    ax7 = fig.add_subplot(gs[3, 0])
    ax8 = fig.add_subplot(gs[4, 0], sharex=ax7)
    ran = (0,1e4)
    before, *_ = ax7.hist(kpipi["piplus_px"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax7.hist(kpipi["piplus_px"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax7.hist(kpi["B_DTF_P2_PX"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})\pi^{\pm}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax7.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax8.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.grid(True)
    ax8.set_ylim(0.8,1.2)
    ax8.set_xlabel(r"$p_{x}$ [MeV/$c$]")
    ax8.set_ylabel(r"$K\pi/K\pi\pi$")


    # pion py
    ax9 = fig.add_subplot(gs[3, 1])
    ax10 = fig.add_subplot(gs[4, 1], sharex=ax9)
    ran = (0,1e4)
    before, *_ = ax9.hist(kpipi["piplus_py"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax9.hist(kpipi["piplus_py"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax9.hist(kpi["B_DTF_P2_PY"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})\pi^{\pm}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax9.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax10.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.grid(True)
    ax10.set_ylim(0.8,1.2)
    ax10.set_xlabel(r"$p_{y}$ [MeV/$c$]")
    ax10.set_ylabel(r"$K\pi/K\pi\pi$")
    

    # pion pz
    ax11 = fig.add_subplot(gs[3, 2])
    ax12 = fig.add_subplot(gs[4, 2], sharex=ax11)
    ran = (0,1e5)
    before, *_ = ax11.hist(kpipi["piplus_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax11.hist(kpipi["piplus_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    target, edges, *_ = ax11.hist(kpi["B_DTF_P2_PZ"], range=ran, histtype="step", linewidth=5, label=r"$(D^{0} \rightarrow K^{\mp}\pi^{\pm})\pi^{\pm}$", weights=kpi["sWeight"], **kw_args);
    plt.setp(ax11.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax12.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.grid(True)
    ax12.set_ylim(0.8,1.2)
    ax12.set_xlabel(r"$p_{z}$ [MeV/$c$]")
    ax12.set_ylabel(r"$K\pi/K\pi\pi$")


    for a in [ax1,ax3,ax7,ax5,ax9,ax11]:
        a.set_ylabel("Normalised entries [a.u.]")
        a.legend(loc="best")

    plt.savefig(f"{str(output_directory)}/kpipi_pxpypz_reweightings.pdf")
    

kpipi_reweighter(kpi, kpipi, kw_args)

def kspi_reweighter(kpipi, kspi, kw_args):

    print("plotting kspi pt, eta, phi...")
    fig = plt.figure(figsize=(45,25))
    gs = fig.add_gridspec(5, 3, height_ratios=[4, 1, 1, 4, 1], hspace=0.1) # middle row is used as a spacer and so not used

    # D transverse momentum
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ran = (0,2e4)
    before, *_ = ax1.hist(kspi["dplus_pt"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$", **kw_args);
    after, *_ = ax1.hist(kspi["dplus_pt"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax1.hist(kpipi["dplus_pt"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})D^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax1.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax2.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.grid(True)
    ax2.set_ylim(0.8,1.2)
    ax2.set_xlabel(r"$p_{T}$ [MeV/$c$]")
    ax2.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # D pseudorapidity
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1], sharex=ax3)
    ran = (1.8,5)
    before, *_ = ax3.hist(kspi["dplus_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$", **kw_args);
    after, *_ = ax3.hist(kspi["dplus_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax3.hist(kpipi["dplus_eta"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})D^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax3.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax4.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.grid(True)
    ax4.set_ylim(0.8,1.2)
    ax4.set_xlabel(r"$\eta$")
    ax4.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # D azimuthal angle
    ax5 = fig.add_subplot(gs[0, 2])
    ax6 = fig.add_subplot(gs[1, 2], sharex=ax5)
    ran = (-np.pi,np.pi)
    before, *_ = ax5.hist(kspi["dplus_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$", **kw_args);
    after, *_ = ax5.hist(kspi["dplus_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax5.hist(kpipi["dplus_phi"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})D^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax5.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax6.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.grid(True)
    ax6.set_ylim(0.8,1.2)
    ax6.set_xlabel(r"$\phi$")
    ax6.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # trigger pion transverse momentum
    ax7 = fig.add_subplot(gs[3, 0])
    ax8 = fig.add_subplot(gs[4, 0], sharex=ax7)
    ran = (0,1e4)
    before, *_ = ax7.hist(kspi["trigger_pi_pt"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax7.hist(kspi["trigger_pi_pt"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax7.hist(kpipi["trigger_pi_pt"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax7.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax8.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.grid(True)
    ax8.set_ylim(0.8,1.2)
    ax8.set_xlabel(r"$p_{T}$ [MeV/$c$]")
    ax8.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # trigger pion pseudorapidity
    ax9 = fig.add_subplot(gs[3, 1])
    ax10 = fig.add_subplot(gs[4, 1], sharex=ax9)
    ran = (1.8,5)
    before, *_ = ax9.hist(kspi["trigger_pi_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax9.hist(kspi["trigger_pi_eta"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax9.hist(kpipi["trigger_pi_eta"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax9.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax10.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.grid(True)
    ax10.set_ylim(0.8,1.2)
    ax10.set_xlabel(r"$\eta$")
    ax10.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # trigger pion azimuthal angle
    ax11 = fig.add_subplot(gs[3, 2])
    ax12 = fig.add_subplot(gs[4, 2], sharex=ax11)
    ran = (-np.pi,np.pi)
    before, *_ = ax11.hist(kspi["trigger_pi_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax11.hist(kspi["trigger_pi_phi"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax11.hist(kpipi["trigger_pi_phi"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax11.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax12.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.grid(True)
    ax12.set_ylim(0.8,1.2)
    ax12.set_xlabel(r"$\phi$")
    ax12.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    for a in [ax1,ax3,ax5,ax7,ax9,ax11]:
        a.set_ylabel("Normalised entries [a.u.]")
        a.legend(loc="best")
    
    plt.savefig(f"{str(output_directory)}/kspi_ptetaphi_reweightings.pdf")
    

    print("plotting kspi px, py, pz...")
    fig = plt.figure(figsize=(45,25))
    gs = fig.add_gridspec(5, 3, height_ratios=[4, 1, 1, 4, 1], hspace=0.1) # middle row is used as a spacer and so not used

    # D px
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ran = (0,1e4)
    before, *_ = ax1.hist(kspi["dplus_px"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$", **kw_args);
    after, *_ = ax1.hist(kspi["dplus_px"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax1.hist(kpipi["dplus_px"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})D^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax1.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax2.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax2.grid(True)
    ax2.set_ylim(0.8,1.2)
    ax2.set_xlabel(r"$p_{x}$ [MeV/$c$]")
    ax2.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # D py
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1], sharex=ax3)
    ran = (0,1e4)
    before, *_ = ax3.hist(kspi["dplus_py"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$", **kw_args);
    after, *_ = ax3.hist(kspi["dplus_py"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax3.hist(kpipi["dplus_py"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})D^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax3.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax4.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax4.grid(True)
    ax4.set_ylim(0.8,1.2)
    ax4.set_xlabel(r"$p_{y} [MeV/$c$]$")
    ax4.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # D pz
    ax5 = fig.add_subplot(gs[0, 2])
    ax6 = fig.add_subplot(gs[1, 2], sharex=ax5)
    ran = (0,2e5)
    before, *_ = ax5.hist(kspi["dplus_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$", **kw_args);
    after, *_ = ax5.hist(kspi["dplus_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})D^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax5.hist(kpipi["dplus_pz"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})D^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax5.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax6.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax6.grid(True)
    ax6.set_ylim(0.8,1.2)
    ax6.set_xlabel(r"$p_{z}$ [MeV/$c$]")
    ax6.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # trigger pion px
    ax7 = fig.add_subplot(gs[3, 0])
    ax8 = fig.add_subplot(gs[4, 0], sharex=ax7)
    ran = (0,1e4)
    before, *_ = ax7.hist(kspi["trigger_pi_px"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax7.hist(kspi["trigger_pi_px"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax7.hist(kpipi["trigger_pi_px"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax7.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax8.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax8.grid(True)
    ax8.set_ylim(0.8,1.2)
    ax8.set_xlabel(r"$p_{x}$ [MeV/$c$]")
    ax8.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # trigger pion py
    ax9 = fig.add_subplot(gs[3, 1])
    ax10 = fig.add_subplot(gs[4, 1], sharex=ax9)
    ran = (0,1e4)
    before, *_ = ax9.hist(kspi["trigger_pi_py"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax9.hist(kspi["trigger_pi_py"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax9.hist(kpipi["trigger_pi_py"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax9.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax10.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax10.grid(True)
    ax10.set_ylim(0.8,1.2)
    ax10.set_xlabel(r"$p_{y}$ [MeV/$c$]")
    ax10.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    # trigger pion pz
    ax11 = fig.add_subplot(gs[3, 2])
    ax12 = fig.add_subplot(gs[4, 2], sharex=ax11)
    ran = (0,2e5)
    before, *_ = ax11.hist(kspi["trigger_pi_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$", **kw_args);
    after, *_ = ax11.hist(kspi["trigger_pi_pz"], range=ran, label=r"$(D^{\pm} \rightarrow K_{s}^{0}\pi^{\pm})\pi^{\pm}$ w.", weights=kspi["weight"], **kw_args);
    target, edges, *_ = ax11.hist(kpipi["trigger_pi_pz"], range=ran, histtype="step", linewidth=5, label=r"$(D^{\pm} \rightarrow K^{\mp}\pi^{\pm}\pi^{\pm})\pi^{\pm}$ w.", weights=kpipi["weight"], **kw_args);
    plt.setp(ax11.get_xticklabels(), visible=False)
    target_to_before_ratio = [target[i]/before[i] for i in range(kw_args["bins"])]
    target_to_after_ratio = [target[i]/after[i] for i in range(kw_args["bins"])]
    centres = [(edges[i]+edges[i+1])/2 for i in range(kw_args["bins"])]
    widths = [(edges[i+1]-edges[i])/2 for i in range(kw_args["bins"])]
    ax12.errorbar(centres, target_to_before_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.errorbar(centres, target_to_after_ratio, xerr=widths, fmt=".", markersize=4, elinewidth=1, capthick=0)
    ax12.grid(True)
    ax12.set_ylim(0.8,1.2)
    ax12.set_xlabel(r"$p_{z}$ [MeV/$c$]")
    ax12.set_ylabel(r"$K\pi\pi/K_{S}^{0}\pi$")


    for a in [ax1,ax3,ax5,ax7,ax9,ax11]:
        a.set_ylabel("Normalised entries [a.u.]")
        a.legend(loc="best")
    
    plt.savefig(f"{str(output_directory)}/kspi_pxpypz_reweightings.pdf")


kspi_reweighter(kpipi, kspi, kw_args)
