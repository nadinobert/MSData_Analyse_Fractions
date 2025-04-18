import data as data
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import Path
import pandas as pd
from scipy.interpolate import make_interp_spline
from matplotlib.transforms import ScaledTranslation
import matplotlib.ticker as ticker
import os

# TODO es ist jetzt aktuell die standardabweichung der aktivitäten in einer extra spalte berücksichtigt. Muss jetzt immer händisch in das csv file eingefügt werden
# TODO es soll aufgefordrt werden Activity werte hinzuzufügen und eine extra spalte dafür eingefügt werden (column 23)

# Set plot font globally to Times New Roman and set globally font size to 20
# locally set fontsize will overwrite globally defined setting!
plt.rcParams.update({'font.family': 'Arial'})

# open csv file of interest but skip the first two rows
directory_path = Path(r"C:\Users\hellmold\Nextcloud\Experiments\Size_exclusion_chromatography\20241031_SEC_AEX_frac") # the r enables the use of \ without interpreting it as escape sign
file_name = '20241031_SEC_AEX_frac_digitonin.csv'
data = pd.read_csv(
    os.path.join(directory_path, file_name),
    skiprows=2, delimiter=',')  # falls mit tabs getrennt '\t' oder';'

# change the headers (3 x mAU) to unique headers for UV absorbance
names = data.columns.tolist()
names[1] = 'UV1_280nm'
names[3] = 'UV2_360nm'
names[5] = 'UV3_410nm'
data.columns = names

df = pd.DataFrame(data, columns=['min', 'UV1_280nm', 'UV2_360nm', 'UV3_410nm'])
df = df.drop_duplicates(subset=['min'], keep='last').dropna()

# parameters for plot/ figure
figure_name = "20241030_SEC_MV"
flowrate = 0.4  # [ml/min]
xmin = 5
xmax = 20
stepx1 = 1  # steps on x-axis
y1min = 0
y1max = 20
stepy1 = 2
y2min = 0
y2max = 20
stepy2 = 2
fraction_size = 0.5

column_type = "big"             #small or big sec column
activity_test = "Dehalogenase"  #Dehalogenase, Hydrogenase or Dehalogenase+Hydrogenase

# define the different columns as plot and x values
# "elution" ist die von "min" in elution volumen umgerechnete spalte
# "frac" ist die spalte mit den gesammelten fraktionen in dem data sheet
elution = df['min'].dropna().apply(lambda x: x * flowrate)
frac = data.iloc[:, [20]].dropna().apply(lambda x: x * flowrate)

# define df for activity for hyd and dehy test (columns: u, v, w, x, y, z)
# define df for activity for hyd or dehy test (columns: u, v, w, x)
if activity_test == 'Dehalogenase+Hydrogenase':
    activity = data.iloc[:, 20:26].dropna()
    activity.columns = ["min", "Fractions", "Activity BV [%]", "STABWN BV", "Activity MV [%]", "STABWN MV"]
    activity['frac_start [ml]'] = frac
    # Masking zero standard deviations
    std_dev_masked_BV = [np.nan if sd == 0 else sd for sd in activity['STABWN BV']]
    std_dev_masked_MV = [np.nan if sd == 0 else sd for sd in activity['STABWN MV']]
elif activity_test == 'Dehalogenase' or activity_test == 'Hydrogenase':
    activity = data.iloc[:, 20:24].dropna()
    activity.columns = ["min", "Fractions", "Activity [%]", "STABWN"]
    activity['frac_start [ml]'] = frac
    # Masking zero standard deviations
    std_dev_masked = [np.nan if sd == 0 else sd for sd in activity['STABWN']]
print(activity)


y1 = df['UV1_280nm'].dropna().apply(lambda x: x +6)
# remove zero values from 360 (and 410nm)
# TODO: es muss verdammt nochmal iwie diese kack baseline zuverlässig von 410 und 360 abgezogen werden
df_nonull = df[df['UV3_410nm'] != 0]

# y2_min = df_nonull['UV2_360nm'].min()
y2_start = df['UV2_360nm'][5]
y2 = df['UV2_360nm'].dropna().apply(lambda x: x -y2_start) #-y2_start
print(y2)

# y3_min= df_nonull['UV3_410nm'].min()
y3_start = df['UV3_410nm'][5]
y3 = df['UV3_410nm'].dropna().apply(lambda x: x -y3_start +1) #-y3_start

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title(figure_name, fontsize=22, y=1.3)

# include second, third and furth y-axis sharing the same x-axis (twinx)
# par1 = host.twinx()
# par2 = host.twinx()
par3 = host.twinx()

# set x and y axis + labels + define tick size
host.set_xlabel("Elution volume (ml)", fontsize=18)
host.set_ylabel("Absorption (mAU)", fontsize=18, color='black')
# Force y-axis to display only integer labels
host.set_ylim(y1min, y1max) #set y min max to avoid that the axis starts automaticalyy at -2 or sth
host.yaxis.set_major_locator(ticker.MaxNLocator(integer=True)) #Ensure y-axis labels are integers
host.yaxis.set_major_locator(ticker.MultipleLocator(stepy1)) #defines the steps in the y axis
host.tick_params(axis='y', labelsize=14)
# par1.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
# par1.tick_params(axis='y', labelsize=14)
# par2.set_ylabel("Absorption at 410 nm [mAU]", fontsize=18)
if activity_test == "Dehalogenase":
    par3.set_ylabel('Relative dehalogenase activity (%)', fontsize=18, color='red')
if activity_test == "Hydrogenase":
    par3.set_ylabel('Relative hydrogenase activity (%)', fontsize=18, color='blue')
if activity_test == "Dehalogenase+Hydrogenase":
    par3.set_ylabel('Relative enzymatic activity (%)', fontsize=18, color='black')
par3.tick_params(axis='y', labelsize=14)

# add a second x-axis
host2 = host.twiny()
host2.set_xlabel("Collected fractions", fontsize=18, labelpad=10)

# add a third x-axis for showing expected molecular weight
host3 = host.twiny()
host3.set_xlabel("Molecular weight (kDa)", fontsize=18, labelpad=10)

host.set_xlim(xmin, xmax)
host.set_ylim(y1min, y1max)
host2.set_xlim(xmin, xmax)
host3.set_xlim(xmin, xmax)
# par1.set_ylim(-10, 100)
# par2.set_ylim(700, 1200)
par3.set_ylim(y2min, y2max) #set y min max to avoid that the axis starts automaticalyy at -2 or sth
par3.yaxis.set_major_locator(ticker.MaxNLocator(integer=True)) #Ensure y-axis labels are integers
par3.yaxis.set_major_locator(ticker.MultipleLocator(stepy2)) #defines the steps in the y axis

# set x1-ticks (elution volume) +define step size
host.xaxis.set_ticks(np.arange(xmin, xmax, stepx1))
host.tick_params(axis='x', labelsize=14)

# set x2 ticks (collected fractions)
x3 = [str(i + 1) for i in range(frac.size)]  # zählt einfach nur durch die anzahl der fractions durch
host2.xaxis.set_ticks(frac['min.10'].values.tolist())  # wo die ticks gesetzt werden
host2.xaxis.set_ticklabels(x3)  # was an den ticks steht
host2.tick_params(axis='x', labelsize=14) #10
host2.set_xlim(xmin, xmax)  # needs to be repeated here?! otherwise out of range

# shift the x-ticks (collected fractions) slightly more right to align labels with ticks
# create offset transform (x=12pt) -> shifting distance
dx, dy = 16, 0 #12
offset = ScaledTranslation(dx / fig.dpi, dy / fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to xticklabels
for label in host2.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# set host3 ticks (Molecular weight [kDa])
host3.xaxis.set_ticks(np.arange(xmin, xmax, stepx1))
# only integer in range() accepted
array = np.arange(xmin, xmax, stepx1)

# calculate molecular sizes according to regression from standard run (27.07.2021) y=-2.8799x + 6.986 !!! Ve=x in ml !!!!
if column_type == "small":
    # small column calibration 20210723
    #x4 = [round(10 ** ((i/1.42) * (-6.1891) + 9.4435)) for i in array]  # int function shows number without decimal places, round function defines decimal places
    x4 = [round(10 ** ((i/1.07) * (-1.6859) + 4.4949)) for i in array]
if column_type == "big":
    #big column approach 20230109 Schebnem
    ##x4 = [int(10 ** (i * (-0.1875) + 4.4974)) for i in array]
    #big column approach 20231204 Digitonin
    x4 = [int(10 ** (i * (-0.1894) + 4.6396)) for i in array]

host3.xaxis.set_ticklabels(x4)
host3.tick_params(axis='x', labelsize=14, rotation=45)

if activity_test == 'Dehalogenase':
    bars = par3.bar(activity['frac_start [ml]'], activity['Activity [%]'], fraction_size,
                    align='edge', color='red', alpha=0.3,
                    edgecolor="black", yerr=std_dev_masked,
                    error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
elif activity_test == 'Hydrogenase':
    bars = par3.bar(activity['frac_start [ml]'], activity['Activity BV'], fraction_size,
                    align='edge', color='blue', alpha=0.3,
                    edgecolor="black", yerr=std_dev_masked,
                    error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
elif activity_test == 'Dehalogenase+Hydrogenase':
    bars1 = par3.bar(activity['frac_start [ml]'], activity['Activity BV [%]'], fraction_size,
                     align='edge', color='blue', alpha=0.3,
                     edgecolor="black", yerr=std_dev_masked_BV,
                     error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
    bars2 = par3.bar(activity['frac_start [ml]'] , activity['Activity MV [%]'], fraction_size,
                     align='edge', color='red', alpha=0.4,
                     edgecolor="black", yerr=std_dev_masked_MV,
                     error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})

# set colors for the absorbance plots
color1 = 'blue'
color2 = 'orangered'
color3 = 'black'

# power_smooth einleiten um plot smooth zu machen und so
xnew = np.linspace(min(elution.tolist()), max(elution.tolist()), 912)
spl1 = make_interp_spline(elution.tolist(), y1, k=3)
power_smooth1 = spl1(xnew)
spl2 = make_interp_spline(elution.tolist(), y2, k=3)
power_smooth2 = spl2(xnew)
spl3 = make_interp_spline(elution.tolist(), y3, k=3)
power_smooth3 = spl3(xnew)

p1, = host.plot(xnew, power_smooth1, color=color1, label="280 nm")
p2, = host.plot(xnew, power_smooth2, color=color2, label="360 nm")
p3, = host.plot(xnew, power_smooth3, color=color3, label="410 nm")

# TODO: get a joined legend of host and par1
# set legend
# px = p1, + p2, + p3,
# labs = [l.get_label() for l in px]
# host.legend(px, labs, loc=10)

host.legend(loc=2)
# par1.legend(loc=0)

# set grid
host.xaxis.grid(color='gray', linestyle='-', linewidth=0.5)
host.yaxis.grid(color='gray', linestyle='-', linewidth=0.5)

# define positions of scales / right, left, top, bottom
par3.spines['right'].set_position(('outward', 0))
host3.spines['top'].set_position(('outward', 50))

# Adjust spacings w.r.t. figsize
fig.tight_layout()

plt.savefig(os.path.join(directory_path, figure_name + '.svg'), format='svg')
plt.show()
fig.savefig(figure_name)
