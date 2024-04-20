"""
model_fitting_model1.py

Bifurcated Bifurcated Johnson

This code is used to fit the data in one of the bims. It then returns the relevant plots of the best fit to the data and a .txt file containing the values and errors on the normalization constant of both signal and background, the mean and standard deviation of the pull distribution and the reduced chi squared value. It can both fit the data using a binned approach or an unbinned one. The model used consists of a Gaussian function and a Crystal Ball function for the signal, and a Chevychev polynomial for the background. Some of the parameters are fixed to be the same as the best-fit values obtained during the global fit in order to obtain better convergence.
The year of interest, size of the data, meson of interest and polarity to be analysed must be specified using the required flags --year --size --meson --polarity. It is also required to specify the bin to be analyzed using the flag --bin, and if the fit should be done on the binned data or the unbinned data using the flag --binned_fit. There also are the flags --input --parameteers_path and --path, which are not required. These are used to specify the directory where the input data is located, where the global best-fit parameters can be found and where the output should be written, respectively. By default it is set to be the current working directory.
Original Code by Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk), which was updated and improved.

Author: Author: Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)

Last edited: 16th September 2023
"""

# - - - - - - IMPORT STATEMENTS - - - - - - #

import ROOT
import argparse
import os
import numpy as np
from ROOT import TChain, RooRealVar, RooDataSet, RooGaussian, RooCrystalBall, RooExponential, RooAddPdf, RooArgList, RooFit, RooArgSet, RooDataHist, RooBifurGauss, RooJohnson, RooArgusBG, RooFormulaVar
from lhcbstyle import LHCbStyle 
import gc
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# - - - - - - - FUNCTIONS - - - - - - - #
def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def size_argument(value):
    if value.isdigit():
        # If the input is a digit, treat it as an integer
        int_value = int(value)
        if 1 <= int_value <= 800 and int_value % 10 == 0: 
            return int_value
        else:
            raise argparse.ArgumentTypeError("Integer value must be between 1 and 800 and be divisible by 10.")
    else:
        raise argparse.ArgumentTypeError("Invalid value. Choose between 'small', 'medium', 'large', or an integer between 1 and 800 that is divisible by 10.")

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size      Used to specify the amount of events the user is interested in analysing.
                The argument must be one of: [1-800]. The integers specify the number of root
                files to be read in.
    --polarity  Used to specify the polarity of the magnet the user is interested in.
                The argument must be one of: [up, down].
    --meson     Used to specify the meson the user is interested in.
                The argument must be one of: [D0, D0bar, both].
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --parameters_path
                Used to specify the directory in which the global best-fit parameters should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --binned_fit
                Used to specify if the data should be binned before performing the fit or an unbinned fit should be performed.
                Type either y or Y for a binned fit. Type n or N for an unbinned fit.
                
    Returns the parsed arguments.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--year",
        type=int,
        choices=[16,17,18],
        required=True,
        help="flag to set the data taking year."
    )
    parser.add_argument(
        "--size",
        type=size_argument,
        required=True,
        help="flag to set the data taking year."
    )
    parser.add_argument(
        "--polarity",
        type=str,
        choices=["up","down"],
        required=True,
        help="flag to set the data taking polarity."
    )
    parser.add_argument(
        "--meson",
        type=str,
        choices=["D0","D0bar","both"],
        required=True,
        help="flag to set the D0 meson flavour."
    )    
    parser.add_argument(
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--parameters_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the global best fit parameters are found"
    )
    parser.add_argument(
        "--input",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the input files should be taken from"
    )
    parser.add_argument(
        "--binned_fit",
        type=str,
        choices=["y", "Y", "n", "N"],
        required=True,
        help="flag to set whether a binned or an unbinned should be performed (y/n)"
    )

    parser.add_argument(
        "--bin",
        type=str,
        required=False,
        help="flag to set whether a binned or an unbinned should be performed (y/n)"
    )

    parser.add_argument(
        "--scheme",
        type=str,
        choices=["total","pT_eta","pT","eta"],
        required=True,
        help="flag to set which binning scheme to use"
    )
    
    return parser.parse_args()
def enableBinIntegrator(func, num_bins):
    """
    Force numeric integration and do this numeric integration with the
    RooBinIntegrator, which sums the function values at the bin centers.
    """
    custom_config = ROOT.RooNumIntConfig(func.getIntegratorConfig())
    custom_config.method1D().setLabel("RooBinIntegrator")
    custom_config.getConfigSection("RooBinIntegrator").setRealValue("numBins", num_bins)
    func.setIntegratorConfig(custom_config)
    func.forceNumInt(True)

def disableBinIntegrator(func):
    """
    Reset the integrator config to disable the RooBinIntegrator.
    """
    func.setIntegratorConfig()
    func.forceNumInt(False)
def gaussian(x, A, mu, sig):
    return A * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def generate_list(size_value_local):
    result_list = []
    current_value = 10

    while not result_list or result_list[-1] < size_value_local:
        result_list.append(current_value)
        current_value += 10

    return result_list
# - - - - - - - MAIN BODY - - - - - - - #

options = parse_arguments()
numbins = 130
lower_boundary = 1820
upper_boundary = 1905
meson = options.meson
polarity = options.polarity
polarity = polarity[0].upper() + polarity[1:]
if options.year == 19:
    if options.scheme == "eta" or options.scheme == "pT" or options.scheme == "pT_eta":
        plot_type = f"Total Mag{polarity} Bin{options.bin}"
    else:
        plot_type = f"Total Mag{polarity}"
else:
    if options.scheme == "eta" or options.scheme == "pT" or options.scheme == "pT_eta":
        plot_type = f"20{options.year} Mag{polarity} Bin{options.bin}"
    else:
        plot_type = f"20{options.year} Mag{polarity}"


if options.binned_fit=="y" or options.binned_fit=="Y":
    binned = True
else:
    binned = False

# Read data
ttree = TChain("D02Kpi_Tuple/DecayTree")

size_value = options.size
size_value = int(size_value)
size_list = generate_list(size_value)
for valuesize in size_list:
    if options.scheme == "total":
        ttree.Add(f"{options.input}/20{options.year}/{options.polarity}/{options.meson}/{options.meson}_{options.polarity}_data_{options.year}_{valuesize}_clean.root")
        ttree.SetBranchStatus("*", 0)
        ttree.SetBranchStatus("D0_MM", 1)

    else:
        #######################
        # Come back to this later
        # Selects invariant mass (D0_MM) of DO
        #Total years Bins
        ttree.Add(f"{options.input}/{options.meson}_{options.size}_bin{options.bin}.root")
        ttree.SetBranchStatus("*", 0)
        ttree.SetBranchStatus("D0_MM", 1)


# Reads in the fit parameters generated by fit_global.py, these will be either for a binned/unbinned fit depending on if fit_global.py was ran as a binned fit or not
parameters = np.loadtxt(f"{options.parameters_path}/fit_parameters.txt", delimiter=',')

ttree.SetBranchStatus("*", 0)
ttree.SetBranchStatus("D0_MM", 1)
D0_M = RooRealVar("D0_MM", r" m($K\pi$) [MeV$\mathrm{c^{-2}}$]", 1815, 1910) # Data - invariant mass

# #Model 1
# Define variables for signal model, using the best fit parameters generated from fit_global.py
# Model Bifurcated Gaussian
bifurmean = RooRealVar("bifurmean", "bifurmean", parameters[25])
sigmaL = RooRealVar("sigmaL", "sigmaL", parameters[13])
sigmaR = RooRealVar("sigmaR", "sigmaR", parameters[14])
Bifurgauss = RooBifurGauss("Bifurgauss", "Bifurgauss", D0_M, bifurmean, sigmaL, sigmaR)

# Model Johnson SU Distribution
Jmu = RooRealVar("Jmu", "Jmu", parameters[21])
Jlam = RooRealVar("Jlam", "Jlam", parameters[22])
Jgam = RooRealVar("Jgam", "Jgam", parameters[23])
Jdel = RooRealVar("Jdel", "Jdel", parameters[24])
Johnson = RooJohnson("Johnson","Johnson", D0_M, Jmu, Jlam, Jgam, Jdel)

# Model Gaussian
mean = RooRealVar("mean", "mean", parameters[31])
sigma = RooRealVar("sigma", "sigma", parameters[32])
gauss = RooGaussian("Gaussian", "Gaussian", D0_M, mean, sigma)

# Model Bifurcated Gaussian
bifurmean2 = RooRealVar("bifurmeam2", "bifurmean2", parameters[26])
sigmaL2 = RooRealVar("sigmaL2", "sigmaL2", parameters[15])
sigmaR2 = RooRealVar("sigmaR2", "sigmaR2", parameters[16])
Bifurgauss2 = RooBifurGauss("Bifurgauss2", "Bifurgauss2", D0_M, bifurmean2, sigmaL2, sigmaR2)

# Model Exponential Background
a = RooRealVar("a0", "a0", parameters[0])
background = RooExponential("Exponential", "Exponential", D0_M, a)

if options.meson == "D0":
    # D0 MagDown
    if options.polarity == "down":
        frac = RooRealVar("frac_D0_down", "frac_D0_down", parameters[1])
        # frac = RooRealVar("frac_D0_down", "frac_D0_down", 0.203)

        frac2 = RooRealVar("frac_D0_down_2", "frac_D0_down_2", parameters[17])
        # frac2 = RooRealVar("frac_D0_down_2", "frac_D0_down_2", 0.2)

        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[5])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[6])
        Nsig_error = parameters[27]
    # D0 MagUp
    elif options.polarity == "up":
        frac = RooRealVar("frac_D0_up", "frac_D0_up", parameters[2])
        frac2 = RooRealVar("frac_D0_up_2", "frac_D0_up_2", parameters[18])
        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[7])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[8])
        Nsig_error = parameters[28]
elif options.meson == "D0bar":
    # D0bar MagDown
    if options.polarity == "down":
        frac = RooRealVar("frac_D0bar_down", "frac_D0bar_down", parameters[3])
        frac2 = RooRealVar("frac_D0bar_down_2", "frac_D0bar_down_2", parameters[19])
        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[9])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[10])
        Nsig_error = parameters[29]
    # D0bar MagUp
    elif options.polarity == "up":
        frac = RooRealVar("frac_D0bar_up", "frac_D0bar_up", parameters[4])
        frac2 = RooRealVar("frac_D0bar_up_2", "frac_D0bar_up_2", parameters[20])
        Nsig = RooRealVar("Nbkg_D0_up", "Nbkg_D0_up", parameters[11])
        Nbkg = RooRealVar("Nbkg_D0_down", "Nbkg_D0_down", parameters[12])
        Nsig_error = parameters[30]



signal = RooAddPdf("signal", "signal", RooArgList(Johnson, Bifurgauss, gauss), RooArgList(frac, frac2))
model = {
    "total": RooAddPdf("total", "Total", RooArgList(signal, background), RooArgList(Nsig, Nbkg)), # extended likelihood
    "signals": {
        Johnson.GetName(): Johnson.GetTitle(),
        Bifurgauss.GetName(): Bifurgauss.GetTitle(),
        gauss.GetName(): gauss.GetTitle(),
    },
    "backgrounds": {
        background.GetName(): background.GetTitle()
    }
}

# Fit data
if binned:
    with LHCbStyle():
        # Creates the histogram for the meson by converting the TTree D0_MM data inside the TChain to a TH1(base class of ROOT histograms)
        # TTree.Draw plots a histogram with name D0_Hist and given bin parameters and saves it to memory using: >>
        ttree.Draw(f"D0_MM>>D0_Hist({numbins},{lower_boundary},{upper_boundary})")
        # D0_Hist recalled from memory and saved to the local variable
        D0_Hist = ROOT.gPad.GetPrimitive("D0_Hist")
        # Creating Binned container sets using RooDataHist
        Binned_data = RooDataHist("Binned_data", "Binned Data Set", RooArgList(D0_M), D0_Hist)
        enableBinIntegrator(signal, D0_M.numBins())
        result = model["total"].fitTo(Binned_data, RooFit.Save(True), RooFit.Extended(True), IntegrateBins = 1e-03)


        frame = D0_M.frame(RooFit.Name(""))
        legend_entries = dict()

        Binned_data.plotOn(frame, ROOT.RooFit.Name("remove_me_A"))
        model["total"].plotOn(
            frame,
            RooFit.Name(model["total"].GetName()),
            RooFit.LineWidth(8),
            RooFit.LineColor(ROOT.kOrange + 1),
        )
        pull_hist = frame.pullHist()

        legend_entries[model["total"].GetName()] = {"title": model["total"].GetTitle(), "style": "l"}

        # plot signal components
        signal_colours = [ROOT.kRed, ROOT.kSpring, ROOT.kAzure + 7]
        signal_line_styles = [2, 7, 9, 10]
        i = 0
        for name, title in model["signals"].items():
            legend_name = f"S{i}"
            model["total"].plotOn(
                frame,
                ROOT.RooFit.Components(name),
                ROOT.RooFit.Name(legend_name),
                ROOT.RooFit.LineWidth(4),
                ROOT.RooFit.LineColor(signal_colours[i % len(signal_colours)]),
                ROOT.RooFit.LineStyle(signal_line_styles[i % len(signal_line_styles)]),
            )
            legend_entries[legend_name] = {"title": title, "style": "l"}
            i += 1

        # plot background components
        background_colours = [ROOT.kMagenta + 2, ROOT.kPink + 7, ROOT.kMagenta + 4]
        background_line_styles = [5, 8, 6]
        i = 0
        for name, title in model["backgrounds"].items():
            legend_name = f"B{i}"
            model["total"].plotOn(
                frame,
                ROOT.RooFit.Components(name),
                ROOT.RooFit.Name(legend_name),
                ROOT.RooFit.LineWidth(4),
                ROOT.RooFit.LineColor(background_colours[i % len(background_colours)]),
                ROOT.RooFit.LineStyle(background_line_styles[i % len(background_line_styles)]),
            )
            legend_entries[legend_name] = {"title": title, "style": "l"}
            i += 1
        
        # plot data points on top again
        Binned_data.plotOn(frame, ROOT.RooFit.Name("remove_me_B"))
        frame.remove("remove_me_A")
        frame.remove("remove_me_B")

        D0_Hist.SetMarkerStyle(20)
        D0_Hist.SetMarkerSize(0.9)
        D0_Hist.SetMarkerColor(ROOT.kBlack)
        D0_Hist.SetTitle("Data")

        # Scale up the error bars by a factor of 5 (you can adjust the scaling factor as needed)
        #if global_local == False:
            # for bin in range(1, D0_Hist.GetNbinsX() + 1):
            #     error = D0_Hist.GetBinError(bin)
            #     D0_Hist.SetBinError(bin, 10 * error)

        frame.addTH1(D0_Hist, "pe")
        legend_entries[D0_Hist.GetName()] = {"title": D0_Hist.GetTitle(), "style": "pe"}


        numbins = D0_Hist.GetNbinsX()
        mD0_bins = []
        for i in range(1, numbins+1):
            mD0_bins.append(D0_Hist.GetBinLowEdge(i))
        mD0_bins.append(D0_Hist.GetBinLowEdge(numbins) + D0_Hist.GetBinWidth(numbins))
        mD0_bins = np.array(mD0_bins, dtype=float)
        frame.SetYTitle(r"Entries [MeVc^{-2}]")
        frame.GetYaxis().SetTitleOffset(0.95)



        c = ROOT.TCanvas("fit", "fit", 900, 800)
        fit_pad = ROOT.TPad("fit_pad", "fit pad", 0, 0.2, 1.0, 1.0)
        fit_pad.SetLogy()
        fit_pad.Draw()
        fit_pad.cd()
        frame.Draw()
        
        frame.GetXaxis().SetLabelSize(0)
        frame.GetXaxis().SetTitleSize(0)



        frame.Draw()
        title_size = frame.GetYaxis().GetTitleSize() * 2.5
        label_size = frame.GetYaxis().GetLabelSize() * 2.5

        # plot_type + total + signals + backgrounds + data
        nlines = 1 + 1 + len(model["signals"]) + len(model["backgrounds"]) + 1
        xwidth = 0.4
        ywidth = 0.04 * nlines

        latex = ROOT.TLatex()
        latex.SetNDC()
        if meson == "D0":
            latex.DrawLatex(0.7, 0.8, "#it{D^{0} #rightarrow K^{-}#pi^{+}}")
        else:
            latex.DrawLatex(0.7, 0.8, "#it{#bar{D}^{0} #rightarrow K^{+}#pi^{-}}")

        # Draw the text on the canvas
        latex.Draw('same')

        legend = ROOT.TLegend(
            0.16, 0.91 - ywidth - 0.05, 0.1 + xwidth, 0.91, "#bf{#it{"+plot_type+"}}"
        )
        legend.SetFillStyle(0)
        #legend.SetFillColor(ROOT.kRed)
        legend.SetBorderSize(0)
        #legend.SetFillColorAlpha(ROOT.kRed, 0.3)
        legend.SetTextSize(label_size*0.30)
        for key, val in legend_entries.items():
            legend.AddEntry(key, val["title"], val["style"])
        legend.Draw("same")

        

        # Plots the pull distribution, where bad pulls (>5 sigma away from the fit) are made to be red
        pull_frame = D0_M.frame(ROOT.RooFit.Title(" "))
        pull_TH1 = ROOT.TH1D("pull_TH1", "pull_TH1", numbins, mD0_bins)
        bad_pull_TH1 = ROOT.TH1D("bad_pull_TH1", "bad_pull_TH1", numbins, mD0_bins)
        for i in range(pull_hist.GetN()):
            if pull_hist.GetPointY(i) > 5:
                pull_TH1.SetBinContent(i + 1, 5)
                bad_pull_TH1.SetBinContent(i + 1, 5)
            elif pull_hist.GetPointY(i) < -5:
                pull_TH1.SetBinContent(i + 1, -5)
                bad_pull_TH1.SetBinContent(i + 1, -5)
            elif pull_hist.GetPointY(i) == 0:
                pull_TH1.SetBinContent(i + 1, 0)
                bad_pull_TH1.SetBinContent(i + 1, 0)
            else:
                pull_TH1.SetBinContent(i + 1, pull_hist.GetPointY(i))
                if abs(pull_hist.GetPointY(i)) >= 3:
                    bad_pull_TH1.SetBinContent(i + 1, pull_hist.GetPointY(i))

        bad_pull_TH1.SetFillColor(ROOT.kRed)
        pull_frame.addTH1(pull_TH1, "bar min0")
        pull_frame.addTH1(bad_pull_TH1, "bar min0")

        c.cd(0)
        pull_pad = ROOT.TPad("pull_pad", "pull pad", 0.0, 0.0, 1.0, 0.31)
        pull_pad.SetBottomMargin(0.4)
        pull_pad.Draw()
        pull_pad.cd()

        # Defines the axes specifications of the pull distribution plot
        pull_frame.GetXaxis().SetLabelSize(label_size)
        pull_frame.GetXaxis().SetTitleSize(title_size)
        pull_frame.GetXaxis().SetTitleOffset(1)
        pull_frame.GetYaxis().SetRangeUser(-5, 5)
        pull_frame.GetYaxis().SetNdivisions(5)
        pull_frame.GetYaxis().SetTitle("Pull [#sigma]")
        pull_frame.GetYaxis().SetLabelSize(label_size)
        pull_frame.GetYaxis().SetTitleSize(title_size)
        pull_frame.GetYaxis().SetTitleOffset(0.35)
        if meson == "D0":
            pull_frame.GetXaxis().SetTitle(r"#it{m(K^{-}\pi^{+})} [MeVc^{-2}]")
        elif meson == "D0bar":
            pull_frame.GetXaxis().SetTitle(r"#it{m(K^{+}\pi^{-})} [MeVc^{-2}]")

        line = ROOT.TLine(D0_M.getMin(), 0, D0_M.getMax(), 0)
        pull_frame.Draw()
        line.Draw("same")

        three = ROOT.TLine(D0_M.getMin(), 3, D0_M.getMax(), 3)
        nthree = ROOT.TLine(D0_M.getMin(), -3, D0_M.getMax(), -3)
        three.SetLineColor(ROOT.kRed)
        three.SetLineStyle(9)
        nthree.SetLineColor(ROOT.kRed)
        nthree.SetLineStyle(9)
        three.Draw("same")
        nthree.Draw("same")


        #Plotting pulls
        pull = ROOT.TH1D("Pulls", "Pulls", 50, -5, 5)
        for i in range(pull_hist.GetN()):
            pull.Fill(pull_hist.GetPointY(i))

        pull_canvas = ROOT.TCanvas("fit_pull", "pull canvas", 600, 400)
        pull.GetXaxis().SetTitle("Distance from fit [#sigma]")
        pull.GetYaxis().SetTitle("Entries")
        pull.Draw('same')



        #Putting Bincentres and Content into Lists
        y = []
        x = []
        for i in range(pull.GetNbinsX()):
            y.append(pull.GetBinContent(i))
            x.append(pull.GetBinCenter(i))

        #Using curve_fit to find Paramters
        try:
            params, cov, *_ = curve_fit(gaussian, x, y, p0=[max(x),0,1], bounds=([0,-np.inf,0],[np.inf,np.inf,np.inf]))
            errs = np.sqrt(np.diag(cov))
            rounded_pull_mean = round(params[1],2)
            rounded_pull_std = round(params[2],2)
            rounded_pull_mean_error = round(errs[1], 2)
            rounded_pull_std_error = round(errs[2], 2)
            str_pull_mean = str(rounded_pull_mean)
            str_pull_sigma = str(rounded_pull_std)
            Failed = 0
        except RuntimeError:
            print("Optimal parameters not found")
            Failed = 1
        if options.year == 19:
            if meson == "D0":
                if options.scheme == "eta" or options.scheme == "pT" or options.scheme == "pT_eta":
                    plot_type2 = f"Total D^{{0}} Mag{polarity} Bin{options.bin}"
                else:
                    plot_type2 = f"Total D^{{0}} Mag{polarity}"
            else:
                if options.scheme == "eta" or options.scheme == "pT" or options.scheme == "pT_eta":
                    plot_type2 = f"Total \\bar{{D}}^{{0}} Mag{polarity} Bin{options.bin}"
                else:
                    plot_type2 = f"Total \\bar{{D}}^{{0}} Mag{polarity}"
        else:
            if meson == "D0":
                if options.scheme == "eta" or options.scheme == "pT" or options.scheme == "pT_eta":
                    plot_type2 = f"20{options.year} D^{{0}} Mag{polarity} Bin{options.bin}"
                else:
                    plot_type2 = f"20{options.year} D^{{0}} Mag{polarity}"
            else:
                if options.scheme == "eta" or options.scheme == "pT" or options.scheme == "pT_eta":
                    plot_type2 = f"20{options.year} \\bar{{D}}^{{0}} Mag{polarity} Bin{options.bin}"
                else:
                    plot_type2 = f"20{options.year} \\bar{{D}}^{{0}} Mag{polarity}"


        if Failed == 0:
            gaussian_fit = ROOT.TF1("gaussian_fit", "gaus", -5, 5)
            gaussian_fit.SetLineColor(4)
            gaussian_fit.SetParameters(params[0],params[1],params[2])
            gaussian_fit.Draw('same')

        # legend2 = ROOT.TLegend(
        #     0.67, 0.67,0.8,0.90, "#bf{#it{"+plot_type2+"}}"
        # )
        if options.scheme == "total":
            legend2 = ROOT.TLegend(
                0.67, 0.74,0.79,0.89, "#bf{#it{"+plot_type2+"}}"
            )
        else:
            legend2 = ROOT.TLegend(
                0.64, 0.74,0.79,0.89, "#bf{#it{"+plot_type2+"}}"
            )


        legend2.SetFillStyle(0)
        legend2.SetBorderSize(0)
        # text1 = f"\\mu: {str(rounded_pull_mean)} \\pm {str(rounded_pull_mean_error)}"
        # text2 = f"\\sigma: {str(rounded_pull_std)} \\pm {str(rounded_pull_std_error)}"
        if options.scheme == "total": 
            legend2.SetTextSize(0.056)
        else:
            legend2.SetTextSize(0.047)

        legend2.AddEntry('Data', 'Data', "l")
        if Failed == 0:
            legend2.AddEntry(gaussian_fit, "Gaussian Fit", "l")
            # legend2.AddEntry(gaussian_fit, text1, "l")
            # legend2.AddEntry(gaussian_fit, text2, "l")
        legend2.Draw("same")


        if Failed == 0:
            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.056)
            if options.scheme == "total":
                latex.DrawLatex(0.67 ,0.70 , 'Fit \mu:  ' + str(rounded_pull_mean) + ' \pm ' + str(rounded_pull_mean_error))
                latex.DrawLatex(0.67,0.65 , 'Fit \sigma:  ' + str(rounded_pull_std) + ' \pm ' + str(rounded_pull_std_error))
            else:
                latex.DrawLatex(0.64 ,0.70 , 'Fit \mu:  ' + str(rounded_pull_mean) + ' \pm ' + str(rounded_pull_mean_error))
                latex.DrawLatex(0.64,0.65 , 'Fit \sigma:  ' + str(rounded_pull_std) + ' \pm ' + str(rounded_pull_std_error))


    


        pull_canvas.Show()

        # nparams = model["total"].getParameters(Binned_data).selectByAttrib("Constant", False).getSize()
        # chi2 = frame.chiSquare(model["total"].GetName(), 'Binned Data Set', nparams - 1)

        print("Saving plots")
        if options.scheme != "total":
            c.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}_bin{options.bin}_fit_ANA.pdf")
            if Failed == 0:
                pull_canvas.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}_bin{options.bin}_fit_pulls.pdf")
            file = open(f"{options.path}/yields_{options.meson}_{options.polarity}_{options.year}_{options.size}_bin{options.bin}.txt", "w+")
            if Failed == 0:
                text = str(Nsig.getValV()) + ', ' + str(Nsig_error) + ', ' + str(Nbkg.getValV()) + ', ' + str(Nbkg.getError()) + ', ' + str_pull_mean + ', ' + str_pull_sigma
            else:
                text = str(Nsig.getValV()) + ', ' + str(Nsig_error) + ', ' + str(Nbkg.getValV()) + ', ' + str(Nbkg.getError())
            file.write(text)
            file.close()
        else:
            c.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}_fit_ANA.pdf")
            if Failed == 0:
                pull_canvas.SaveAs(f"{options.path}/{options.meson}_{options.polarity}_{options.year}_{options.size}_fit_pulls.pdf")
            file = open(f"{options.path}/yields_{options.meson}_{options.polarity}_{options.year}_{options.size}.txt", "w+")
            # Nsig Nsig_Err NBkg NBkg error pull_mean pull_sigma
            if Failed == 0:
                text = str(Nsig.getValV()) + ', ' + str(Nsig_error) + ', ' + str(Nbkg.getValV()) + ', ' + str(Nbkg.getError()) + ', ' + str_pull_mean + ', ' + str_pull_sigma
            else:
                text = str(Nsig.getValV()) + ', ' + str(Nsig_error) + ', ' + str(Nbkg.getValV()) + ', ' + str(Nbkg.getError())                
            file.write(text)
            file.close()
        disableBinIntegrator(signal)
        
        

#print(ttree.GetEntries())
gc.collect()
exit()