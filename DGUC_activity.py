import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
from matplotlib.transforms import ScaledTranslation
import matplotlib.ticker as ticker

# Set plot font globally to Times New Roman and set globally font size to 20
# locally set fontsize will overwrite globally defined setting!
plt.rcParams.update({'font.family': 'Arial'})

# open csv file of interest but skip the first two rows
data = pd.read_csv(
    r'C:\Users\hellmold\Nextcloud\Experiments\Activity_Assay_Photometric\Dehalogenase\20231220_DGUC_Rebecca\20231220_BV_MV.csv',
    skiprows=0, delimiter=',')  # falls mit tabs getrennt '\t' oder';'

# parameters for plot/ figure
figure_name = "20231220_DGUC_BV_MV"
fraction_size = 1
activity_test = 'Dehalogenase+Hydrogenase'

xmin = 0
xmax = 15
ymin = 0
ymax = 23.5

# define df for activity for hyd and dehy test (columns: a, b, c, d, e)
# define df for activity for hyd or dehy test (columns: a, b, c)
if activity_test == 'Dehalogenase+Hydrogenase':
    activity = data.iloc[:, 0:5].dropna()
    activity.columns = ["DGUC fraction", "Activity BV [%]", "STABWN BV", "Activity MV [%]", "STABWN MV"]
    # Masking zero standard deviations
    std_dev_masked_BV = [np.nan if sd == 0 else sd for sd in activity['STABWN BV']]
    std_dev_masked_MV = [np.nan if sd == 0 else sd for sd in activity['STABWN MV']]
elif activity_test == 'Dehalogenase' or activity_test == 'Hydrogenase':
    activity = data.iloc[:, 0:3].dropna()
    activity.columns = ["DGUC fraction","Activity [%]", "STABWN"]
    # Masking zero standard deviations
    std_dev_masked = [np.nan if sd == 0 else sd for sd in activity['STABWN']]
print(activity)

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title(figure_name, fontsize=22, y=1.3)

# set x and y axis + labels + define tick size + tick shift
host.set_xlabel("Collected fractions", fontsize=18)
host.set_ylabel("Relative enzymatic activity (%)", fontsize=18, color='black')
host.tick_params(axis='y', labelsize=13)
host.tick_params(axis='x', labelsize=13)

host.set_xlim(xmin, xmax)
host.set_ylim(ymin, ymax)

dx, dy = 40, 0 #12
offset = ScaledTranslation(dx / fig.dpi, dy / fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to xticklabels
for label in host.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

if activity_test == "Dehalogenase":
    host.set_ylabel('Relative dehalogenase activity (%)', fontsize=18, color='red')
if activity_test == "Hydrogenase":
    host.set_ylabel('Relative hydrogenase activity (%)', fontsize=18, color='blue')
if activity_test == "Dehalogenase+Hydrogenase":
    host.set_ylabel('Relative enzymatic activity (%)', fontsize=18, color='black')

if activity_test == 'Dehalogenase':
    bars = host.bar(activity['DGUC fraction'], activity['Activity [%]'], fraction_size,
                    align='edge', color='red', alpha=0.4,
                    edgecolor="black", yerr=std_dev_masked,
                    error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
elif activity_test == 'Hydrogenase':
    bars = host.bar(activity['DGUC fraction'], activity['Activity BV'], fraction_size,
                    align='edge', color='blue', alpha=0.3,
                    edgecolor="black", yerr=std_dev_masked,
                    error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
elif activity_test == 'Dehalogenase+Hydrogenase':
    bars1 = host.bar(activity['DGUC fraction'], activity['Activity BV [%]'], fraction_size,
                     align='edge', color='blue', alpha=0.4,
                     edgecolor="black", yerr=std_dev_masked_BV,
                     error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
    bars2 = host.bar(activity['DGUC fraction'] , activity['Activity MV [%]'], fraction_size,
                     align='edge', color='red', alpha=0.5,
                     edgecolor="black", yerr=std_dev_masked_MV,
                     error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})

# set grid
host.xaxis.grid(color='gray', linestyle='-', linewidth=0.5)
host.yaxis.grid(color='gray', linestyle='-', linewidth=0.5)

# Adjust spacings w.r.t. figsize
fig.tight_layout()

plt.savefig(r"C:\Users\hellmold\Nextcloud\Experiments\Activity_Assay_Photometric\Dehalogenase\20231220_DGUC_Rebecca\DGUC_BV_MV.svg", format='svg')
plt.show()
fig.savefig(figure_name)

