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

unbinned_asymm, unbinned_asymm_error = read_asymmetry_values()

asymmetry = np.array([])
asymmetry_error = np.array([])
for j in range(0,10):
    bin_num = str(j)
    with open(f'{args.asymm_path}/asymmetries_{args.year}_{args.size}_bin{bin_num}.txt') as f:
        lines = f.readlines()
        A_prod = float(lines[0])
        A_prod_err = float(lines[1])
    f.close()
    asymmetry = np.append(asymmetry, A_prod)
    asymmetry_error = np.append(asymmetry_error, A_prod_err)

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

# initialise globally so can print after
store = {}
def get_chi2(c: float):
    chi2 = 0.
    ndof = -1 # 1 for c

    for i in range(bins):
        if asymmetry_error[i] == 0:
            continue

        pol0 = c
        residual = (asymmetry[i] - pol0) / asymmetry_error[i]

        chi2 += residual * residual
        ndof += 1
    
    store["chi2"] = chi2
    store["ndof"] = ndof

    return chi2



m = Minuit(get_chi2, c=0) # 0 is initial guess
result = m.migrad()
values = list(m.values)
errors = list(m.errors)
prob = 1 - chi2.cdf(store["chi2"], store["ndof"])

result_string = "\n".join((
    rf"$c = {values[0]:.3f} \pm {errors[0]:.3f}, m=0$",
    rf"$\chi^{{2}} / $nDOF$ = {store['chi2']:.2f} / {store['ndof']}$",
    rf"$p = {prob*100:.1f} \%$"
))

_, ax = plt.subplots()

ax.errorbar(x_value, asymmetry, yerr=asymmetry_error, fmt="k.", label="data")
ax.plot(x_value,[values[0] for _ in range(bins)],color="b", label="fit")
# ax.text(0.55, 0.1, result_string, transform=ax.transAxes, fontsize=30)
ax.set_xlabel(r"$p_T$ [GeV$/c$]")
ax.set_ylabel("Production asymmetry")
ax.legend(loc="best")
plt.savefig("asym_flatness.pdf", bbox_inches="tight")


###################

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

ax1.set_ylabel(r'$A_{\mathrm{prod}}$ [%]', fontsize = 50)
ax1.tick_params(axis='both', which='both', labelsize=30)

line2 = ax1.axhline(unbinned_asymm, color='red', linestyle='solid', linewidth=5)
fill2 = ax1.axhspan(unbinned_asymm-unbinned_asymm_error, unbinned_asymm+unbinned_asymm_error, color='red', alpha=0.35, lw=0)

Data = ax1.errorbar(x_value, asymmetry, yerr=asymmetry_error,xerr=x_value_error, fmt='o', capsize=5, color = 'black', label = 'Data')

# Fill between simulated error bars
for i in range(len(x_value)):
    ax1.fill_between([x_value[i] - x_value_error[i], x_value[i] + x_value_error[i]],
                    simulated_asymmetry[i] - simulated_asymmetry_error[i],
                    simulated_asymmetry[i] + simulated_asymmetry_error[i],
                    color='green', alpha=0.4, linewidth=0)
    
extra = Rectangle((0, 0), 1, 1, fc="green", fill=True, edgecolor='none', linewidth=0, alpha=0.4)

# Fit = ax.axhline(values[0], color='purple', linestyle=':', linewidth=5)

if args.scheme == 'pT':
    ax1.legend([(line2,fill2),Data, extra],[r'Bin integrated result','Data','Pythia'])#, loc='upper right')
elif args.scheme == 'eta':
    ax1.legend([(line2,fill2),Data, extra],[r'Bin integrated result','Data','Pythia'])#,='upper right')

file_path = f"{args.path}/result_of_fit.txt"
with open(file_path, "w") as file:
  file.write(f"{result_string}")

ratio = (simulated_asymmetry / asymmetry)
ratio_error = np.sqrt((simulated_asymmetry_error/asymmetry)**2 + ((simulated_asymmetry*asymmetry_error)/(asymmetry)**2)**2)

# Plot residuals against observed data with error bars
for i, (x, y, xerr, yerr) in enumerate(zip(x_value, ratio, x_value_error, ratio_error)):
    color = 'black' if (y + yerr < 1) and (y - yerr < 1) else 'black'
    ax2.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color=color)

# ax2.errorbar(x_value, ratio, xerr=x_value, yerr=ratio_error, fmt='o', color='black')

ax2.axhline(y=1, color='grey', linestyle='--') 
# ax2.axhline(y=0, color='grey', linestyle='--')  
# ax2.axhline(y=-3, color='red', linestyle='--')  

if args.scheme == 'pT':
    ax2.set_xlabel(r'$p_{T}$ [GeV$c^{-1}$]', fontsize = 50)
elif args.scheme == 'eta':
    ax2.set_xlabel(r'$\eta$', fontsize = 50)
ax2.set_ylabel(r'Ratio [Pythia/asymmetry]', fontsize = 50)



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
ax[0].plot(yKS,cKS1)
ax[0].plot(yKS,cKS2)
ax[1].plot(yKS,dKS)
# Add more ticks on the y-axis
ax[0].set_yticks(np.arange(0, 1 + 0.1, 0.2))
ax[1].set_yticks(np.arange(0, max(dKS) + 0.1, 0.1))
# Add more ticks on the x-axis
ax[0].set_xticks(np.arange(round(min_value_asymmetry,1)-0.1, round(max_value_asymmetry,1) + 0.1, 0.2))
ax[1].set_xticks(ax[0].get_xticks())
ax[0].set_ylabel('Cumulative Distribution Frequency')
ax[1].set_xlabel('Production Asymmetry')
ax[1].set_ylabel(r'{/delta} CDF ')
y,x=max(zip(dKS,yKS))
ax[1].plot(x,y, 'o')
plt.savefig("KS-Test_y_axis.pdf",bbox_inches="tight" )

# print output
n1 = len(asymmetry)
n2 = len(simulated_asymmetry)
print('The KS test output for {0} and {1} entries is D={2:.3f} and d={3:.3f}.'.format(n1, n2, max(dKS), max(dKS)*(n1*n2/(n1+n2))**0.5))

ks_statistic, p_value = ks_2samp(asymmetry, simulated_asymmetry)

# Print the test statistic and p-value
print("Test Statistic:", ks_statistic)
print("p-value:", p_value)