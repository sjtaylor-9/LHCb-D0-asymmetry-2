import os
import argparse
import numpy as np
from lhcbstyle import LHCbStyle
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.colors import ListedColormap
from matplotlib import colormaps
import awkward as ak
import sys
from iminuit import Minuit
import matplotlib.pyplot as plt
import mplhep; mplhep.style.use("LHCb2")
from scipy.stats import chi2, ks_2samp
from matplotlib.patches import Rectangle

# - - - - - - - FUNCTIONS - - - - - - - #

def parse_arguments():
    '''
    Parses the arguments needed along the code. Arguments:
    
    --year      Used to specify the year at which the data was taken the user is interested in.
                The argument must be one of: [16, 17, 18]. These referr to 2016, 2017 & 2018, respectively.
    --size      Used to specify the amount of events the user is interested in analysing.
                The argument must be one of: [large, small, medium, 1-8]. The integers specify the number of root
                files to be read in. Large is equivalent to 8. Medium is equivalent to 4. Small takes 200000 events.
    --polarity  Used to specify the polarity of the magnet the user is interested in.
                The argument must be one of: [up, down].
    --meson     Used to specify the meson the user is interested in.
                The argument must be one of: [D0, D0bar, both].
    --path      Used to specify the directory in which the output files should be written. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --input     Used to specify the directory in which the input data should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    --bin_path  Used to specify the directory in which the binning scheme should be found. It is not required,
                in the case it is not specified, the default path is the current working directory.
    
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
        type=str,
        required=True,
        help="flag to set the data taking size."
    )
    parser.add_argument(
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--bin_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the binning scheme should be found"
    )
    parser.add_argument(
        "--sim_bin_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the binning scheme should be found"
    )
    parser.add_argument(
        "--asymm_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the production asymmetry for each bin should be found"
    )
    parser.add_argument(
        "--sim_asymm_path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the production asymmetry for each bin should be found"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        choices=["pT","eta"],
        required=True,
        help="flag to set whether a binned or an unbinned should be performed (y/n)"
    )
    return parser.parse_args()

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def read_asymmetry_values():
    with open(f'{args.path}/final_asymmetries_{args.scheme}_{args.year}_{args.size}.txt') as f:
        lines = f.readlines()
        unbinned_asymm = float(lines[0])
        unbinned_asymm_error = float(lines[1])
        f.close()
    return unbinned_asymm, unbinned_asymm_error
def getCumY(yy, yKS):
    # sort the data and build cumulative distribution on the set of y axis points
    yy.sort()
    inc = 1. / len(yy)
    cumY = 0
    n = 0
    cKS = []
    for y in yy:
        while n < len(yKS) and y > yKS[n]:
            cKS.append(cumY)
            n += 1
        cumY += inc
    for nn in range(nSamples+1-len(cKS)):
        cKS.append(cumY)
    return cKS
# - - - - - - - MAIN BODY - - - - - - - #
args = parse_arguments()
scheme = args.scheme
unbinned_asymm, unbinned_asymm_error = read_asymmetry_values()

asymmetry = np.array([])
asymmetry_error = np.array([])
for j in range(0,10):
    bin_num = str(j)
    with open(f'{args.asymm_path}/asymmetries_{args.year}_{args.size}_bin{bin_num}_detection_scheme_{scheme}.txt') as f:
        lines = f.readlines()
        A_prod = float(lines[0])
        A_prod_err = float(lines[1])
    f.close()
    asymmetry = np.append(asymmetry, A_prod)
    asymmetry_error = np.append(asymmetry_error, A_prod_err)

asymmetry_global_detection = np.array([])
asymmetry_error_global_detection = np.array([])
for j in range(0,10):
    bin_num = str(j)
    with open(f'{args.asymm_path}/asymmetries_{args.year}_{args.size}_bin{bin_num}_detection_scheme_global.txt') as f:
        lines = f.readlines()
        A_prod = float(lines[0])
        A_prod_err = float(lines[1])
    f.close()
    asymmetry_global_detection = np.append(asymmetry_global_detection, A_prod)
    asymmetry_error_global_detection = np.append(asymmetry_error_global_detection, A_prod_err)

simulated_asymmetry = np.array([])
simulated_asymmetry_error = np.array([])
for j in range(0,10):
    bin_num = str(j)
    with open(f'{args.sim_asymm_path}/asymmetries_pythia_{args.scheme}_bin{bin_num}.txt') as f:
        lines = f.readlines()
        sim_A_prod = float(lines[0])
        sim_A_prod_err = float(lines[1])
    f.close()
    simulated_asymmetry = np.append(simulated_asymmetry, sim_A_prod)
    simulated_asymmetry_error = np.append(simulated_asymmetry_error, sim_A_prod_err)


x_value =np.array([])
x_value_error = np.array([])

file_path = f"{args.bin_path}/{args.year}_{args.size}_{args.scheme}_bins.txt"  # Replace with the path to your file
# Open the file in read mode
with open(file_path, 'r') as file:
    # Read all lines from the file and store them in a list
    bin_lines = [float(line.strip()) for line in file.readlines()]

print(bin_lines)
for i in range(0,10):
    # Get center of bin and width of bin as error
    x_value_indivual = (bin_lines[i]+bin_lines[i+1])/2
    x_value_error_indivual = (bin_lines[i+1]) - ((bin_lines[i]+bin_lines[i+1])/2)
    x_value_error = np.append(x_value_error, x_value_error_indivual)
    x_value = np.append(x_value, x_value_indivual)


if args.scheme == 'pT':
    x_value = [x / 1000 for x in x_value] # in Gev
    x_value_error = [x / 1000 for x in x_value_error] # in Gev


##########
bins = len(x_value)

# Plotting
fig = plt.figure()
ax1 = fig.gca()

ax1.set_ylabel(r'$A_{\mathrm{prod}}$ [%]', fontsize = 50)
ax1.tick_params(axis='both', which='both', labelsize=30)

line2 = ax1.axhline(unbinned_asymm, color='red', linestyle='solid', linewidth=5)
fill2 = ax1.axhspan(unbinned_asymm-unbinned_asymm_error, unbinned_asymm+unbinned_asymm_error, color='red', alpha=0.35, lw=0)

Data = ax1.errorbar(x_value, asymmetry, yerr=asymmetry_error,xerr=x_value_error, fmt='o', capsize=5, color = 'black', label = r'Data (Inc. local $A_{\mathrm{det}}$)', alpha=0.75)
Data2 = ax1.errorbar(x_value, asymmetry_global_detection, yerr=asymmetry_error_global_detection, xerr=x_value_error, fmt='o', capsize=5, color = 'blue', label = r'Data (Inc. global $A_{\mathrm{det}}$)', alpha=0.6)

# Fill between simulated error bars
for i in range(len(x_value)):
    ax1.fill_between([x_value[i] - x_value_error[i], x_value[i] + x_value_error[i]],
                    simulated_asymmetry[i] - simulated_asymmetry_error[i],
                    simulated_asymmetry[i] + simulated_asymmetry_error[i],
                    color='green', alpha=0.4, linewidth=0)
    
extra = Rectangle((0, 0), 1, 1, fc="green", fill=True, edgecolor='none', linewidth=0, alpha=0.4)

# Fit = ax.axhline(values[0], color='purple', linestyle=':', linewidth=5)

if args.scheme == 'pT':
    ax1.legend([(line2,fill2),Data, Data2, extra],[r'Bin integrated result', r'Data (Inc. local $A_{\mathrm{det}}$)', r'Data (Inc. global $A_{\mathrm{det}}$)', 'Pythia'], loc='upper left')
elif args.scheme == 'eta':
    ax1.legend([(line2,fill2),Data, Data2, extra],[r'Bin integrated result', r'Data (Inc. local $A_{\mathrm{det}}$)', r'Data (Inc. global $A_{\mathrm{det}}$)','Pythia'], loc='best')

if args.scheme == 'pT':
    ax1.set_xlabel(r'$p_{T}$ [GeV$c^{-1}$]', fontsize = 50)
elif args.scheme == 'eta':
    ax1.set_xlabel(r'$\eta$', fontsize = 50)

if args.scheme == 'pT':
    plt.savefig(f'{args.path}/pT_Asymm_{args.year}_{args.size}.pdf', bbox_inches = "tight")
elif args.scheme == 'eta':
    plt.savefig(f'{args.path}/eta_Asymm_{args.year}_{args.size}.pdf', bbox_inches = "tight")

plt.show()

min_value_asymmetry = min(asymmetry.min(), simulated_asymmetry.min())
max_value_asymmetry = max(asymmetry.max(), simulated_asymmetry.max())
print(min_value_asymmetry)
print(max_value_asymmetry)
# y axis points for plotting and evaluation of the difference
nSamples = 10000
yKS = np.linspace(min_value_asymmetry, max_value_asymmetry, nSamples+1)


cKS1 = getCumY(asymmetry, yKS)
cKS2 = getCumY(simulated_asymmetry, yKS)

# calculate difference dataset, the maximum of which is the KS test statistic
dKS = [abs(x-y) for x,y in zip(cKS1, cKS2)]
# plotting just cdf

fig,ax = plt.subplots(2,1,figsize=(16, 8))

ax[0].plot(yKS,cKS1, label = 'Measured Data')
ax[0].plot(yKS,cKS2, label = 'Simulated Data')
ax[1].plot(yKS,dKS, color = "green", label = r'$\Delta$ CDF')
# Add more ticks on the y-axis
ax[0].set_yticks(np.arange(0, 1 + 0.2, 0.25))
ax[1].set_yticks(np.arange(0, max(dKS) + 0.2, 0.25))
# Add more ticks on the x-axis
if max_value_asymmetry > 1.5 or min_value_asymmetry < 1.5:
        ax[0].set_xticks(np.arange(round(min_value_asymmetry,1)-0.1, round(max_value_asymmetry,1) + 0.1, 0.4))
else:
    ax[0].set_xticks(np.arange(round(min_value_asymmetry,1)-0.1, round(max_value_asymmetry,1) + 0.1, 0.2))
ax[1].set_xticks(ax[0].get_xticks())
ax[0].set_ylabel(r'CDF',fontsize=50)
ax[1].set_xlabel(r'Production Asymmetry [%]',fontsize=50)
ax[1].set_ylabel(r'|$\Delta$ CDF| ',fontsize=50)
# Adding legends
y,x=max(zip(dKS,yKS))
ax[1].plot(x,y, 'o', color = 'red', label = r'Maximum $\Delta$ CDF value')
ax[0].legend(fontsize=25, loc='best')
ax[1].legend(fontsize=25, loc='best')
if args.scheme == 'pT':
    plt.savefig(f'{args.path}/pT_KS-Test CDF and DeltaCDF.pdf', bbox_inches = "tight")
elif args.scheme == 'eta':
    plt.savefig(f'{args.path}/eta_KS-Test CDF and DeltaCDF.pdf.pdf', bbox_inches = "tight")

# print output
n1 = len(asymmetry)
n2 = len(simulated_asymmetry)
# Calculate KS test statistic and p-value
ks_statistic, p_value = ks_2samp(asymmetry, simulated_asymmetry)

# Open a text file in write mode
if args.scheme == 'pT':
    filename_text_output = f'{args.path}/eta_KS-Test CDF and DeltaCDF.txt'
elif args.scheme == 'eta':
    filename_text_output = f'{args.path}/pT_KS-Test CDF and DeltaCDF.txt'
    
with open(filename_text_output, "w") as f:
    # Write the test statistic and p-value to the file
    f.write("KS-Test Statistic from scipy: {}\n".format(ks_statistic))
    f.write("p-value from scipy: {}\n".format(p_value))

    # Print the test statistic and p-value to console
    print("KS-Test Statistic from scipy:", ks_statistic)
    print("p-value from scipy:", p_value)

    # Print the KS test output to console and write it to the file
    output = "The KS test output for {0} and {1} entries is D={2:.3f} and d={3:.3f}.\n".format(n1, n2, max(dKS), max(dKS)*(n1*n2/(n1+n2))**0.5)
    print(output)
    f.write(output)