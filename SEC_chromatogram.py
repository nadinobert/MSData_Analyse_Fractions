import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
from matplotlib.transforms import ScaledTranslation

# TODO check if delimiter ; oder tab. ändert sich dauernd was zur hölle???
# TODO es soll aufgefordrt werden Activity werte hinzuzufügen und eine extra spalte dafür eingefügt werden (column 23)

# Set plot font globally to Times New Roman and set globally font size to 20
# locally set fontsize will overwrite globally defined setting!
plt.rcParams.update({'font.family': 'Times New Roman'})

# open csv file of interest but skip the first two rows
data = pd.read_csv(
    r'C:\Users\hellmold\Nextcloud\Experiments\Size_exclusion_chromatography\20240123_SEC_Digitonin\20240123_SEC_Digitonin.csv',
    skiprows=2, delimiter=';')  # falls mit tabs getrennt '\t' oder';'

# change the headers (3 x mAU) to unique headers for UV absorbance
names = data.columns.tolist()
names[1] = 'UV1_280nm'
names[3] = 'UV2_360nm'
names[5] = 'UV3_410nm'
data.columns = names

df = pd.DataFrame(data, columns=['min', 'UV1_280nm', 'UV2_360nm', 'UV3_410nm'])
df = df.drop_duplicates(subset=['min'], keep='last').dropna()

# parameters for plot/ figure
figure_name = "20240123_SEC_Digitonin"
flowrate = 0.4  # [ml/min]
xmin = 6
xmax = 24
step = 2  # steps on x-axis
ymin = -4
ymax = 100
fraction_size = 1

column_type = "big"             #small or big
activity_test = "Dehalogenase"  #Dehalogenase, Hydrogenase


# define the different columns as plot and x values
# "elution" ist die von "min" in elution volumen umgerechnete spalte
# "frac" ist die spalte mit den gesammelten fraktionen in dem data sheet
elution = df['min'].dropna().apply(lambda x: x * flowrate)
frac = data.iloc[:, [20]].dropna().apply(lambda x: x * flowrate)

# define df for activity (columns: u, v, w, x)
activity = data.iloc[:, 20:23].dropna()
activity.columns = ["min", "Fractions", "Activity [%]"]
activity['frac_start [ml]'] = frac

y1 = df['UV1_280nm'].dropna().apply(lambda x: x + 2)
# remove zero values from 360 (and 410nm)
# TODO: es muss verdammt nochmal iwie diese kack baseline zuverlässig von 410 und 360 abgezogen werden
df_nonull = df[df['UV3_410nm'] != 0]

# y2_min = df_nonull['UV2_360nm'].min()
y2_start = df['UV2_360nm'][5]
y2 = df['UV2_360nm'].dropna().apply(lambda x: x - y2_start)
print(y2)

# y3_min= df_nonull['UV3_410nm'].min()
y3_start = df['UV3_410nm'][5]
y3 = df['UV3_410nm'].dropna().apply(lambda x: x - y3_start)

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title(figure_name, fontsize=25, y=1.3)

# include second, third and furth y-axis sharing the same x-axis (twinx)
# par1 = host.twinx()
# par2 = host.twinx()
par3 = host.twinx()

# set x and y axis + labels + define tick size
host.set_xlabel("Elution volume [ml]", fontsize=18)
host.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
host.tick_params(axis='y', labelsize=14)
# par1.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
# par1.tick_params(axis='y', labelsize=14)
# par2.set_ylabel("Absorption at 410 nm [mAU]", fontsize=18)
if activity_test == "Dehalogenase":
    par3.set_ylabel('Dehalogenase Activity [%]', fontsize=18, color='red')
if activity_test == "Hydrogenase":
    par3.set_ylabel('Hydrogenase Activity [%]', fontsize=18, color='blue')
par3.tick_params(axis='y', labelsize=14)
par3.tick_params(axis='y', labelsize=14)

# add a second x-axis
host2 = host.twiny()
host2.set_xlabel("Collected fractions", fontsize=18, labelpad=10)

# add a third x-axis for showing expected molecular weight
host3 = host.twiny()
host3.set_xlabel("Molecular weight [kDa]", fontsize=18, labelpad=10)

host.set_xlim(xmin, xmax)
host.set_ylim(ymin, ymax)
host2.set_xlim(xmin, xmax)
host3.set_xlim(xmin, xmax)
# par1.set_ylim(-10, 100)
# par2.set_ylim(700, 1200)
# par3.set_ylim(ymin, ymax)

# set x1-ticks (elution volume)
host.xaxis.set_ticks(np.arange(xmin, xmax, step))
host.tick_params(axis='x', labelsize=14)

# set x2 ticks (collected fractions)
x3 = [str(i + 1) for i in range(frac.size)]  # zählt einfach nur durch die anzahl der fractions durch
host2.xaxis.set_ticks(frac['min.10'].values.tolist())  # wo die ticks gesetzt werden
host2.xaxis.set_ticklabels(x3)  # was an den ticks steht
host2.tick_params(axis='x', labelsize=10) #10
host2.set_xlim(xmin, xmax)  # needs to be repeated here?! otherwise out of range

# shift the x-ticks (collected fractions) slightly more right to align labels with ticks
# create offset transform (x=18pt) -> shifting distance
dx, dy = 7, 0 #12
offset = ScaledTranslation(dx / fig.dpi, dy / fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to xticklabels
for label in host2.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# set host3 ticks (Molecular weight [kDa])
host3.xaxis.set_ticks(np.arange(xmin, xmax, step))
# only integer in range() accepted
array = np.arange(xmin, xmax, step)

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

# set colors
color1 = 'blue'
color2 = 'orange'
color3 = 'grey'

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
if activity_test == 'Dehalogenase':
    par3.bar(activity['frac_start [ml]'], activity['Activity [%]'], fraction_size, align='edge', color='red', alpha=0.3,
             edgecolor="black")
if activity_test == 'Hydrogenase':
    par3.bar(activity['frac_start [ml]'], activity['Activity [%]'], fraction_size, align='edge', color='blue', alpha=0.3,
             edgecolor="black")
#par3.bar(activity['frac_start [ml]'], activity['Activity [%]'], fraction_size, align='edge', color='red', alpha=0.3,
#         edgecolor="black")
# par3.bar(list(map(lambda i: x2.tolist()[i - 1], frac)), act, 1, align='edge', color='grey', alpha=0.55,
#         edgecolor="black")

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

plt.show()
fig.savefig(figure_name)
