import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
from matplotlib.transforms import ScaledTranslation

from SEC_chromatogram import new_names
from functions import get_flowrate_changes, time_to_elution_volume, remove_duplicates_preserve_order

# TODO check if delimiter ; oder tab. ändert sich dauernd was zur hölle???
# TODO es muss aufgefordrt werden Activity werte hinzuzufügen und eine extra spalte dafür eingefügt werden (column 23)
# TODO zwei spalten für aktivitätstest in das csv file einfügen
# TODO in df elution ergeben sich doppelte werte die dann nicht geplottet werden können!!! ganz schlecht!! rundungsproblem bei berechnung

# Set plot font globally to Times New Roman and set globally font size to 20
# locally set fontsize will overwrite globally defined setting!
plt.rcParams.update({'font.family': 'Arial'})


activity_test = 'Hydrogenase'

figure_name = '20250806_AC'

fraction_size = 1

xmin = 4
xmax = 20
step = 2  # steps on x-axis
ymin = 3
# ymin = df['UV1_280nm'].min() + 5
# ymax = df['UV1_280nm'].max() + 5
ymax = 15

# open the csv file of interest but skip the first two rows
data = pd.read_csv(
    r'C:\Users\hellmold\Nextcloud\Experiments\Affinity_chromatography\20250820_Exp2\20250820_Exp2_AC_StrepTactin_XT_Flow.csv',
    skiprows=2, delimiter=',')  # falls mit tabs getrennt '\t' oder';' regex doesnt work


new_names = {
    data.columns[0]: 'min_280nm',
    data.columns[1]: 'UV1_280nm',
    data.columns[6]: 'min_410nm',
    data.columns[7]: 'UV2_410nm',
    data.columns[8]: 'min_450nm',
    data.columns[9]: 'UV3_450nm',
    data.columns[4]: 'min_%B',
    data.columns[5]: '%B',
    data.columns[20]: 'time_point',
    data.columns[21]: 'flow_rate',
    data.columns[45]: 'Activity (nkat)',
    data.columns[46]: 'Standarddeviation (nkat)'
}


# Rename the columns
data = data.rename(columns=new_names)

#chromatogram_df = pd.DataFrame(data, columns=['min_280nm', 'UV1_280nm', 'min_410nm', 'UV2_410nm','min_450nm', 'UV3_450nm', 'min_%B', '%B'])


df_280 = data[['min_280nm', 'UV1_280nm']].dropna().drop_duplicates(subset=['min_280nm'], keep='last')
df_410 = data[['min_410nm', 'UV2_410nm']].dropna().drop_duplicates(subset=['min_410nm'], keep='last')
df_450 = data[['min_450nm', 'UV3_450nm']].dropna().drop_duplicates(subset=['min_450nm'], keep='last')
df_B = data[['min_%B', '%B']].dropna().drop_duplicates(subset=['min_%B'], keep='last')

# df_flowrate contains timepoints and flowrate per timepoint
# by calling get_flowrate_changes the change of flowrate is assigned to a timepoint
df_flowrate = data.iloc[:, 20:22].dropna()
flowrate_list_not_filtered = get_flowrate_changes(df_flowrate)
flowrate_list = [row for row in flowrate_list_not_filtered if row[1] != 0.0] ## hier müssen die flow rate = 0 einträge gelöscht werden, da ansonsten doppelte werte entstehen


# "elution" ist die von "min" in elution volumen umgerechnete spalte
# calc elution volume for every timepoint in the experiment
elution_280 = df_280['min_280nm'].dropna().apply(lambda x, f=flowrate_list: time_to_elution_volume(x, f))
elution_410 = df_410['min_410nm'].dropna().apply(lambda x, f=flowrate_list: time_to_elution_volume(x, f))
elution_450 = df_450['min_450nm'].dropna().apply(lambda x, f=flowrate_list: time_to_elution_volume(x, f))
elution_B = df_B['min_%B'].dropna().apply(lambda x, f=flowrate_list: time_to_elution_volume(x, f))


# "frac" ist die spalte mit den gesammelten fraktionen in dem data sheet
# calc elution start volume for every fraction
df_frac = pd.DataFrame(data.iloc[:, [40]].dropna())
series_frac = df_frac['min.21'].dropna()
frac = series_frac.apply(lambda x, f=flowrate_list: time_to_elution_volume(x, f))
#print(frac)

# define df for activity (columns: u, v, w (activity), x (STABWN))
activity_df = data.iloc[:, 20:24].dropna()
activity_df['frac_start [ml]'] = frac
activity_df.drop(activity_df[(activity_df['(Fractions)'] == 'Waste')].index, inplace=True)
#print(activity_df['(Fractions)'])
print(activity_df)

y1 = df_280['UV1_280nm']
# remove zero values from 360 (and 410nm)
# TODO: es muss verdammt nochmal iwie diese kack baseline zuverlässig von 410 und 360 abgezogen werden

y2_start = df_410['UV2_410nm'][3]
y2 = df_410['UV2_410nm'].dropna().apply(lambda x: x - y2_start +2)

y3_start = df_450['UV3_450nm'][3]
y3 = df_450['UV3_450nm'].dropna().apply(lambda x: x - y3_start +2)

y4 = data['%B'].dropna()
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
host.set_xlabel("Elution volume (ml)", fontsize=18)
host.set_ylabel("Absorption (mAU)", fontsize=18, color='black')
host.tick_params(axis='y', labelsize=14)
# par1.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
# par1.tick_params(axis='y', labelsize=14)

if activity_test == "Dehalogenase":
    par3.set_ylabel('Dehalogenase activity (nkat)', fontsize=18, color='red')
if activity_test == "Hydrogenase":
    par3.set_ylabel('Hydrogenase activity (nkat)', fontsize=18, color='blue')
par3.tick_params(axis='y', labelsize=14)
par4.set_ylabel('Concentration B (%)', fontsize=18, color='black')
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
dx, dy = 14, 0
offset = ScaledTranslation(dx / fig.dpi, dy / fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to xticklabels
for label in host2.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# set colors
color1 = 'blue'
color2 = 'orangered'
color3 = 'grey'
color4 = 'black'

# Clean up the data before spline interpolation
elution_280_clean, y1_clean = remove_duplicates_preserve_order(elution_280, y1)
elution_410_clean, y2_clean = remove_duplicates_preserve_order(elution_410, y2)
elution_450_clean, y3_clean = remove_duplicates_preserve_order(elution_450, y3)
elution_B_clean, y4_clean = remove_duplicates_preserve_order(elution_B, y4)

# Now for the spline interpolation, use these separate series
xnew_280 = np.linspace(min(elution_280_clean), max(elution_280_clean), 912)
spl1 = make_interp_spline(elution_280_clean, y1_clean, k=3)
power_smooth1 = spl1(xnew_280)

xnew_410 = np.linspace(min(elution_410_clean), max(elution_410_clean), 912)
spl2 = make_interp_spline(elution_410_clean, y2_clean, k=3)
power_smooth2 = spl2(xnew_410)

xnew_450 = np.linspace(min(elution_450_clean), max(elution_450_clean), 912)
spl3 = make_interp_spline(elution_450_clean, y3_clean, k=3)
power_smooth3 = spl3(xnew_450)

xnew_B = np.linspace(min(elution_B_clean), max(elution_B_clean), 912)
spl4 = make_interp_spline(elution_B_clean, y4_clean, k=3)
power_smooth4 = spl4(xnew_B)

p1, = host.plot(xnew_280, power_smooth1, color=color1, label="280 nm")
p2, = host.plot(xnew_410, power_smooth2, color=color2, label="410 nm")
p3, = host.plot(xnew_450, power_smooth3, color=color3, label="450 nm")
p4, = par4.plot(xnew_B, power_smooth4, color=color4, label="Concentration B (%)")

# Masking zero standard deviations
std_dev_masked = [np.nan if sd == 0 else sd for sd in activity_df['STABWN']]

if activity_test == 'Dehalogenase':
    par3.bar(activity_df['frac_start [ml]'], activity_df['Activity [%]])'], fraction_size, align='edge', color='red', alpha=0.3,
             edgecolor="black", yerr=std_dev_masked, error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
    par3.set_ylim(bottom=0)
if activity_test == 'Hydrogenase':
    par3.bar(activity_df['frac_start [ml]'], activity_df['Activity [%]'], fraction_size, align='edge', color='blue', alpha=0.3,
             edgecolor="black", yerr=std_dev_masked, error_kw={'ecolor': 'black', 'capsize': 5, 'capthick': 1, 'elinewidth': 1})
    par3.set_ylim(bottom=0)


# Collect legend handles and labels from host and par4
handles_host, labels_host = host.get_legend_handles_labels()
handles_par4, labels_par4 = par4.get_legend_handles_labels()

# Combine the handles and labels
combined_handles = handles_host + handles_par4
combined_labels = labels_host + labels_par4

# Plot the combined legend on the host axis
host.legend(combined_handles, combined_labels, loc='upper right')  # Adjust 'loc' as needed


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

plt.savefig('Aktueller_plot', format='svg')
fig.savefig(figure_name)
plt.show()

