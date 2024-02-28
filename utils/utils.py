import logging
logging.basicConfig(level=logging.NOTSET)
logger = logging.getLogger("utils")

import argparse
def get_config() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--year",
        type=str,
        choices=["16", "17", "18"],
        required=True
    )
    parser.add_argument(
        "--polarity",
        type=str,
        choices=["u", "d"],
        required=True
    )
    parser.add_argument(
        "--sign",
        type=str,
        choices=["RS", "WS"],
        required=True
    )
    parser.add_argument(
        "--background",
        action="store_true"
    )
    parser.add_argument(
        "--sideband",
        action="store_true"
    )
    parser.add_argument(
        "--mc",
        action="store_true"
    )

    args=parser.parse_args()

    return args

def get_parser() -> argparse.ArgumentParser:
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
    return parser

def gaussian(x, A, mu, sig):
    return A * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

import matplotlib.pyplot as plt

def setup_mpl():
    # must annoying warnings
    plt.set_loglevel(level = "warning")
    import logging
    logging.getLogger("matplotlib.pyplot").disabled = True
    logging.getLogger("matplotlib.font_manager").disabled = True
    logging.getLogger("matplotlib.ticker").disabled = True
    logging.getLogger("matplotlib.backends.backend_pdf").disabled = True
    logging.getLogger("matplotlib.dviread").disabled = True
    logging.getLogger("matplotlib.texmanager").disabled = True
    logging.getLogger("matplotlib._type1font").disabled = True
    logging.getLogger("matplotlib.colorbar").disabled = True
    logging.getLogger("PIL.PngImagePlugin").disabled = True
    logging.getLogger("progressbar.utils").disabled = True

def multiplot(nrows: int, ncols: int, ):
    fig = plt.figure(figsize=((ncols + (ncols-1)*0.35)*plt.rcParams["figure.figsize"][0],(nrows + (nrows-1)*0.35)*plt.rcParams["figure.figsize"][1]))
    gs = fig.add_gridspec(nrows, ncols, height_ratios=[1 for i in range(nrows)], hspace=0.35, width_ratios=[1 for i in range(ncols)], wspace=0.35)

    axs = np.ndarray([nrows, ncols], dtype=plt.Axes)
    for i in range(nrows):
        for j in range(ncols):
            axs[i][j] = fig.add_subplot(gs[i,j], autoscale_on=True)
    
    return fig, axs

import numpy as np
def plot(
    observable,#: R.RooRealVar,
    data,#: R.RooAbsData,
    model: dict,
    plot_type: str,
    bins: int|np.ndarray = 50,
    ran: str = "FULL",
    setlogy: bool = False,
    rescaley: None | tuple[float,float] = None,
    save_to: None | str = None,
    unit: str = "",
) -> None:
    import ROOT as R
    from ROOT import gROOT, RooGaussian
    from lhcbstyle import LHCbStyle
    from scipy.optimize import curve_fit
    
    # everything below is now in LHCb sytle
    with LHCbStyle():
        observable.setRange("FULL", observable.getMin(), observable.getMax())

        # Draw fit with pulls beneath
        fit_frame = observable.frame(R.RooFit.Title(" "))
        legend_entries = dict()

        # plot data points
        if type(bins) == int:
            nbins = bins
            # data.plotOn(fit_frame, R.RooFit.Name(data.GetName()), R.RooFit.Binning(bins), R.RooFit.CutRange(ran))
            data.plotOn(fit_frame, R.RooFit.Name(data.GetName()), R.RooFit.Binning(bins))
        else:
            nbins = len(bins) - 1
            binning = R.RooBinning(nbins, np.array(bins))

            data.plotOn(
                fit_frame,
                R.RooFit.Name(data.GetName()),
                R.RooFit.Binning(binning),
                R.RooFit.DataError(R.RooAbsData.SumW2),
                R.RooFit.Range(ran)
            )

        # nData = data.sumEntries("", ran)
        # Make clear that the target normalisation is nData. The enumerator NumEvent
        # is needed to switch between relative and absolute scaling.
        # model.plotOn(frame01, Normalization(nData, RooAbsReal::NumEvent),
        # ProjectionRange("SB1"));
        # plot model
        model["total"].plotOn(
            fit_frame,
            R.RooFit.Name(model["total"].GetName()),
            R.RooFit.LineWidth(5),
            R.RooFit.LineColor(R.kAzure),
            R.RooFit.NormRange(ran),
            # R.RooFit.Range("FULL")
            R.RooFit.Range(ran),
            # R.RooFit.Normalization(nData, R.RooAbsReal.NumEvent),
            # R.RooFit.ProjectionRange(ran)
        )
        legend_entries[model["total"].GetName()] = model["total"].GetTitle()
        pull_hist = fit_frame.pullHist()

        # plot signal components
        signal_colours = [R.kRed, R.kSpring, R.kAzure + 7, R.kOrange + 7]
        signal_line_styles = [2, 7, 9, 10]
        i = 0
        for name, title in model["signals"].items():
            legend_name = f"S{i}"
            model["total"].plotOn(
                fit_frame,
                R.RooFit.Components(name),
                R.RooFit.Name(legend_name),
                R.RooFit.LineWidth(4),
                R.RooFit.LineColor(signal_colours[i % len(signal_colours)]),
                R.RooFit.LineStyle(signal_line_styles[i % len(signal_line_styles)]),
                R.RooFit.NormRange(ran),
                R.RooFit.Range(ran)
                # R.RooFit.Range("FULL")
                # R.RooFit.Normalization(nData, R.RooAbsReal.NumEvent),
                # R.RooFit.ProjectionRange(ran)
            )
            legend_entries[legend_name] = title
            i += 1

        # plot background components
        background_colours = [R.kMagenta + 2, R.kPink + 7, R.kMagenta + 4]
        background_line_styles = [5, 8, 6]
        i = 0
        for name, title in model["backgrounds"].items():
            legend_name = f"B{i}"
            model["total"].plotOn(
                fit_frame,
                R.RooFit.Components(name),
                R.RooFit.Name(legend_name),
                R.RooFit.LineWidth(4),
                R.RooFit.LineColor(background_colours[i % len(background_colours)]),
                R.RooFit.LineStyle(background_line_styles[i % len(background_line_styles)]),
                R.RooFit.NormRange(ran),
                R.RooFit.Range(ran)
                # R.RooFit.Range("FULL")
                # R.RooFit.Normalization(nData, R.RooAbsReal.NumEvent),
                # R.RooFit.ProjectionRange(ran)
            )
            legend_entries[legend_name] = title
            i += 1

        # have to repeat it down here else the curves disappear...
        if type(bins) == int:
            plotted_data = data.plotOn(fit_frame, R.RooFit.Name(data.GetName()), R.RooFit.Binning(bins))
            # plotted_data = data.plotOn(fit_frame, R.RooFit.Name(data.GetName()), R.RooFit.Binning(bins), R.RooFit.CutRange(ran))
        else:
            plotted_data = data.plotOn(
                fit_frame,
                R.RooFit.Name(data.GetName()),
                R.RooFit.Binning(binning),
                R.RooFit.DataError(R.RooAbsData.SumW2),
            )

        fit_frame.SetYTitle(
            f"Candidates / ({(observable.getRange('Full').second-observable.getRange('Full').first)/nbins:.2f} {unit})"
        )

        fit_canvas = R.TCanvas("fit", "fit canvas", 1125, 800)
        fit_canvas.Draw()
        fit_pad = R.TPad("fit_pad", "fit pad", 0, 0, 0.8, 1.0)
        fit_pad.Draw()
        if setlogy:
            fit_pad.SetLogy()
        if rescaley is not None:
            fit_frame.GetYaxis().SetRangeUser(*rescaley)

        fit_pad.cd()

        # make plot for paper first
        fit_frame.GetYaxis().SetMaxDigits(len(str(int(plotted_data.GetMaximum()/1000)))+1)  # TODO : a better fix than this?
        if plotted_data.GetMaximum()<20 and "dM_fit_bin_" in save_to:
            # this should pick up the case of sideband binned swapped dM fits
            # where nentries is very low
            # goal is to not have legend overlapping with data point error bars
            fit_frame.GetYaxis().SetRangeUser(0, 20)
        fit_frame.Draw()
        # the pull_pad object is not 1x1 but 0.31x1 so the text gets shrunk
        # 2.5 seems to multiple the lhcbstyle font size back to approx the right size
        title_size = fit_frame.GetYaxis().GetTitleSize() * 2.5
        label_size = fit_frame.GetYaxis().GetLabelSize() * 2.5

        # plot_type + total + signals + backgrounds
        nlines = 1 + 1 + len(model["signals"]) + len(model["backgrounds"])
        yheight = 0.04 * nlines
        legend = R.TLegend(0.2, 0.9 - yheight, 0.3, 0.9, "#bf{#it{"+plot_type+"}}")
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.037)
        for key, val in legend_entries.items():
            legend.AddEntry(key, val, "l")
        legend.AddEntry(data.GetName(), data.GetTitle(), "PE")
        # Moves the legend to outside the plot
        #fit_canvas.cd()
        # legend_pad = R.TPad("legend_pad", "legend pad", 0.8, 0.2, 1.0, 1.0)
        # legend_pad.Draw()
        # if setlogy:
        #     legend_pad.SetLogy()
        # legend_pad.cd()
        legend.Draw("same")

        if save_to is not None:
            logger.info(f"Saving plots to {save_to}_PAPER.root")
            fit_canvas.SaveAs(save_to + "_PAPER.root")
            fit_canvas.SaveAs(save_to + "_PAPER.pdf")
            fit_canvas.SaveAs(save_to + "_PAPER.jpg")
            fit_canvas.Print(save_to + "_PAPER.tex")

        # make plot for ana note next first
        fit_frame.GetXaxis().SetLabelSize(0)
        fit_frame.GetXaxis().SetTitleSize(0)
        fit_canvas.cd()
        fit_pad.cd()
        fit_pad.SetPad(0, 0.2, 0.8, 1.0)
        fit_frame.Draw()
        #legend_pad.cd()
        legend.Draw("same")

        # now pulls beneath
        # pull_hist = fit_frame.pullHist(data.GetName(), model["total"].GetName())
        pull_frame = observable.frame(R.RooFit.Title(" "))
        # pull_frame.addPlotable(pull_hist, "P") # this makes points with error bars - we want bar chart
        # go about it in an awkward way:
        
        if type(bins) == int:
            bin_width = pull_hist.GetPointX(1) - pull_hist.GetPointX(0)
            minimum = pull_hist.GetPointX(0) - bin_width / 2
            maximum = pull_hist.GetPointX(pull_hist.GetN() - 1) + bin_width / 2
            pull_TH1 = R.TH1D("pull_TH1", "pull_TH1", pull_hist.GetN(), minimum, maximum)
            bad_pull_TH1 = R.TH1D("bad_pull_TH1", "bad_pull_TH1", pull_hist.GetN(), minimum, maximum)
        else:
            pull_TH1 = R.TH1D("pull_TH1", "pull_TH1", nbins, np.array(bins))
            bad_pull_TH1 = R.TH1D("bad_pull_TH1", "bad_pull_TH1", nbins, np.array(bins))

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
            
        bad_pull_TH1.SetFillColor(R.kRed)
        pull_frame.addTH1(pull_TH1, "bar min0")
        pull_frame.addTH1(bad_pull_TH1, "bar min0")

        fit_canvas.cd(0)
        pull_pad = R.TPad("pull_pad", "pull pad", 0.0, 0.0, 0.8, 0.31)
        pull_pad.SetBottomMargin(0.4)
        pull_pad.Draw()
        pull_pad.cd()

        pull_frame.GetXaxis().SetLabelSize(label_size)
        pull_frame.GetXaxis().SetTitleSize(title_size)
        pull_frame.GetXaxis().SetTitleOffset(1)
        pull_frame.GetYaxis().SetRangeUser(-5, 5)
        pull_frame.GetYaxis().SetNdivisions(5)
        pull_frame.GetYaxis().SetTitle("Pull [#sigma]")
        pull_frame.GetYaxis().SetLabelSize(label_size)
        pull_frame.GetYaxis().SetTitleSize(title_size)
        pull_frame.GetYaxis().SetTitleOffset(0.39)
        
        line = R.TLine(observable.getMin(), 0, observable.getMax(), 0)
        pull_frame.Draw()
        line.Draw("same")

        three = R.TLine(observable.getMin(), 3, observable.getMax(), 3)
        nthree = R.TLine(observable.getMin(), -3, observable.getMax(), -3)
        three.SetLineColor(R.kRed)
        three.SetLineStyle(9)
        nthree.SetLineColor(R.kRed)
        nthree.SetLineStyle(9)
        three.Draw("same")
        nthree.Draw("same")

        pull = R.TH1D("pull", "pull", 50, -5, 5)
        for i in range(pull_hist.GetN()):
            pull.Fill(pull_hist.GetPointY(i))

        pull_canvas = R.TCanvas("fit_pull", "pull canvas", 600, 400)
        pull.GetXaxis().SetTitle("Distance from fit [#sigma]")
        pull.Draw()
        height = pull.GetMaximum()
        print(height)
        
        #Putting Bincentres and Content into Lists
        y = []
        x = []
        for i in range(pull.GetNbinsX()):
            y.append(pull.GetBinContent(i))
            x.append(pull.GetBinCenter(i))
        
        #Using curve_fit to find Paramters
        params, cov, *_ = curve_fit(gaussian, x, y, p0=[max(x),0,1], bounds=([0,-np.inf,0],[np.inf,np.inf,np.inf]))
        errs = np.sqrt(np.diag(cov))
        
        rounded_pull_mean = round(params[1],2)
        rounded_pull_std = round(params[2],2)
        rounded_pull_mean_error = round(errs[1], 2)
        rounded_pull_std_error = round(errs[2], 2)


        gaussian_fit = R.TF1("gaussian_fit", "gaus", -5, 5)
        gaussian_fit.SetLineColor(4)
        gaussian_fit.SetParameters(params[0],params[1],params[2])
        gaussian_fit.Draw('same')

        legend2 = R.TLegend(
            0.68, 0.78,0.8,0.90, "#bf{#it{"+plot_type+"}}"
        )

        legend2.SetFillStyle(0)
        legend2.SetBorderSize(0)

        legend2.SetTextSize(0.04)
        legend2.AddEntry(data.GetName(), data.GetTitle(), "l")
        legend2.AddEntry(gaussian_fit, "Gaussian Fit", "l")
        legend2.Draw("same")
        latex = R.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.68 ,0.75 , 'pull \mu: ' + str(rounded_pull_mean) + ' \pm ' + str(rounded_pull_mean_error))
        latex.DrawLatex(0.68 ,0.71 , 'pull \sigma: ' + str(rounded_pull_std) + ' \pm ' + str(rounded_pull_std_error))
        
        if save_to is not None:
            logger.info(f"Saving plots to {save_to}_ANA.root")
            fit_canvas.SaveAs(save_to + "_ANA.root")
            fit_canvas.SaveAs(save_to + "_ANA.pdf")
            fit_canvas.SaveAs(save_to + "_ANA.jpg")
            fit_canvas.Print(save_to + "_ANA.tex")
            pull_canvas.SaveAs(save_to + "_pulls.root")
            pull_canvas.SaveAs(save_to + "_pulls.pdf")
            pull_canvas.SaveAs(save_to + "_pulls.jpg")
            pull_canvas.Print(save_to + "_pulls.tex")
        
        nparams = model["total"].getParameters(data).selectByAttrib("Constant", False).getSize()
        red_chi2 = fit_frame.chiSquare(model["total"].GetName(), data.GetName(), nparams - 1)
        logger.info(f"Reduced chi-squared for {nbins} bins and {nparams} free parameters: {red_chi2}")
