import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
from matplotlib.transforms import ScaledTranslation
from functions import get_flowrate_changes, time_to_elution_volume

# TODO check if delimiter ; oder tab. ändert sich dauernd was zur hölle???
# TODO es muss aufgefordrt werden Activity werte hinzuzufügen und eine extra spalte dafür eingefügt werden (column 23)
# TODO zwei spalten für aktivitätstest in das csv file einfügen
# TODO in df elution ergeben sich doppelte werte die dann nicht geplottet werden können!!! ganz schlecht!! rundungsproblem bei berechnung

activity_test = 'Dehalogenase'

figure_name = '20240111_AEX_Digitonin'

fraction_size = 1

xmin = 0
xmax = 30
step = 5  # steps on x-axis
ymin = -2

# ymin = df['UV1_280nm'].min() + 5
# ymax = df['UV1_280nm'].max() + 5
ymax = 50

# open csv file of interest but skip the first two rows
all_data = pd.read_csv(
    r'C:\Users\hellmold\Nextcloud\Experiments\Anion_exchange_chromatography\20240111_AEX_Digitonin.csv',
    skiprows=2, delimiter=';')  # falls mit tabs getrennt '\t' oder';' regex doesnt work

# change the headers (3 x mAU) to unique headers for UV absorbance
all_data.rename(columns={'min.8': 'time_point', ' ml/min': 'flow_rate'}, inplace=True)

names = all_data.columns.tolist()
names[1] = 'UV1_280nm'
names[3] = 'UV2_360nm'
names[5] = 'UV3_410nm'
names[11] = '%B'
names[22] = 'Activity [%]'
all_data.columns = names

chromatogram_df = pd.DataFrame(all_data, columns=['min', 'UV1_280nm', 'UV2_360nm', 'UV3_410nm', '%B'])
chromatogram_df = chromatogram_df.drop_duplicates(subset=['min'], keep='last').dropna()

# df_flowrate contains timepoints and flowrate per timepoint
# by calling get_flowrate_changes the change of flowrate is assigned to a timepoint
flowrate_df = all_data.iloc[:, 16:18].dropna()
flowrate_list_not_filtered = get_flowrate_changes(flowrate_df)
flowrate_list = [row for row in flowrate_list_not_filtered if row[1] != 0.0] ## hier müssen die flow rate = 0 einträge gelöscht werden, da ansonsten doppelte werte entstehen
print(flowrate_list)

# "elution" ist die von "min" in elution volumen umgerechnete spalte
# calc elution volume for every timepoint in the experiment
df_elution = pd.DataFrame(chromatogram_df.iloc[:, [0]].dropna())
series_elution = df_elution['min'].dropna()
elution = series_elution.apply(lambda x, f=flowrate_list: time_to_elution_volume(x, f))

# "frac" ist die spalte mit den gesammelten fraktionen in dem data sheet
# calc elution start volume for every fraction
df_frac = pd.DataFrame(all_data.iloc[:, [20]].dropna())
series_frac = df_frac['min.10'].dropna()
frac = series_frac.apply(lambda x, f=flowrate_list: time_to_elution_volume(x, f))
print(frac)

# define df for activity (columns: u, v, w (activity)
activity_df = all_data.iloc[:, 20:23].dropna()
activity_df['frac_start [ml]'] = frac
activity_df.drop(activity_df[(activity_df['(Fractions)'] == 'Waste')].index, inplace=True)

y1 = chromatogram_df['UV1_280nm']
# remove zero values from 360 (and 410nm)
# TODO: es muss verdammt nochmal iwie diese kack baseline zuverlässig von 410 und 360 abgezogen werden

y2_start = chromatogram_df['UV2_360nm'][3]
y2 = chromatogram_df['UV2_360nm'].dropna().apply(lambda x: x - y2_start)

y3_start = chromatogram_df['UV3_410nm'][3]
y3 = chromatogram_df['UV3_410nm'].dropna().apply(lambda x: x - y3_start)

y4 = all_data['%B'].dropna()
y4_min = 0
y4_max = 100

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title(figure_name, fontsize=25, y=1.3)

# include second, third and furth y-axis sharing the same x-axis (twinx)
# par1 = host.twinx()
# par2 = host.twinx()
par3 = host.twinx()
par4 = host.twinx()

# set x and y axis + labels + define tick size
host.set_xlabel("Elution volume [ml]", fontsize=18)
host.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
host.tick_params(axis='y', labelsize=14)
# par1.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
# par1.tick_params(axis='y', labelsize=14)
if activity_test == "Dehalogenase":
    par3.set_ylabel('Dehalogenase Activity [%]', fontsize=18, color='red')
if activity_test == "Hydrogenase":
    par3.set_ylabel('Hydrogenase Activity [%]', fontsize=18, color='blue')
par3.tick_params(axis='y', labelsize=14)
par4.set_ylabel('Concentration B [%]', fontsize=18, color='black')
par4.tick_params(axis='y', labelsize=14)

# add a second x-axis
host2 = host.twiny()
host2.set_xlabel("Collected fractions", fontsize=18, labelpad=10)
host.set_xlim(xmin, xmax)
host.set_ylim(ymin, ymax)
host2.set_xlim(xmin, xmax)
# par1.set_ylim(-10, 100)
par4.set_ylim(y4_min, y4_max)
# par2.set_ylim(700, 1200)
# par3.set_ylim(ymin, ymax)

# set x1-ticks
host.xaxis.set_ticks(np.arange(xmin, xmax, step))
host.tick_params(axis='x', labelsize=14)

# TODO: automaddisch die unnötigen fractionen aus df löschen ohne die händisch aus csv file entfernen zu müssen
# set x2-ticks (collected fractions) -> delete useless fractions in csv file manually
col_frac_to_list = activity_df['(Fractions)'].tolist()
x3 = col_frac_to_list

host2.xaxis.set_ticks(activity_df['frac_start [ml]'].values.tolist())  # wo die ticks gesetzt werden
host2.xaxis.set_ticklabels(x3)  # was an den ticks steht
host2.tick_params(axis='x', labelsize=14)
host2.set_xlim(xmin, xmax)  # needs to be repeated here?! otherwise out of range

# shift the x-ticks (collected fractions) slightly more right to align labels with ticks
# create offset transform (x=18pt) -> shifting distance
dx, dy = 18, 0
offset = ScaledTranslation(dx / fig.dpi, dy / fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to xticklabels
for label in host2.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# set colors
color1 = 'blue'
color2 = 'orange'
color3 = 'grey'
color4 = 'black'

# power_smooth einleiten um plot smooth zu machen und so
xnew = np.linspace(min(elution), max(elution), 910)
spl1 = make_interp_spline(elution, y1, k=3)
power_smooth1 = spl1(xnew)
spl2 = make_interp_spline(elution, y2, k=3)
power_smooth2 = spl2(xnew)
spl3 = make_interp_spline(elution, y3, k=3)
power_smooth3 = spl3(xnew)
spl4 = make_interp_spline(elution, y4, k=3)
power_smooth4 = spl4(xnew)

p1, = host.plot(xnew, power_smooth1, color=color1, label="280 nm")
p2, = host.plot(xnew, power_smooth2, color=color2, label="360 nm")
p3, = host.plot(xnew, power_smooth3, color=color3, label="410 nm")
if activity_test == 'Dehalogenase':
    par3.bar(activity_df['frac_start [ml]'], activity_df['Activity [%]'], fraction_size, align='edge', color='red', alpha=0.3,
             edgecolor="black")
if activity_test == 'Hydrogenase':
    par3.bar(activity_df['frac_start [ml]'], activity_df['Activity [%]'], fraction_size, align='edge', color='blue', alpha=0.3,
             edgecolor="black")
p4, = par4.plot(xnew, power_smooth4, color=color4, label="Concentration B [%]")

# TODO: get a joined legend of host and par1
# set legend
# px = p1, + p2, + p3,
# labs = [l.get_label() for l in px]
# host.legend(px, labs, loc=10)

# host.legend(loc=1)
par4.legend(loc=1)
host.legend(loc=2)

# set grid
host.xaxis.grid(color='gray', linestyle='-', linewidth=0.5)
host.yaxis.grid(color='gray', linestyle='-', linewidth=0.5)

# define positions of scales / right, left, top, bottom
# par3.spines['right'].set_position(('outward', 60))
par4.spines['right'].set_position(('outward', 80))

# Adjust spacings w.r.t. figsize
fig.tight_layout()

# host.yaxis.label.set_color(p1.get_color())
# par1.yaxis.label.set_color(p2.get_color())
# par2.yaxis.label.set_color(p3.get_color())

plt.show()
fig.savefig(figure_name)
