import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
from matplotlib.transforms import ScaledTranslation

# open csv file of interest but skip the first two rows
data = pd.read_csv(
    r'C:\Users\hellmold\Nextcloud\Experiments\Size_exclusion_chromatography\20211124_Reactor_Effl_Membrane\20211124_SEC_Reactor_Effl.csv',
    skiprows=2, delimiter='\t')  # falls mit tabs getrennt '\t' oder';'

#column 23 contains activity values
x2 = data['Activity'].dropna()
#print(x2)
print(x2.item())
x22 = x2.tolist()
#print(x22)

# change the headers (3 x mAU) to unique headers for UV absorbance
names = data.columns.tolist()
names[1] = 'UV1_280nm'
names[3] = 'UV2_360nm'
names[5] = 'UV3_410nm'
data.columns = names

df = pd.DataFrame(data, columns=['min', 'UV1_280nm', 'UV2_360nm', 'UV3_410nm'])
df = df.drop_duplicates(subset=['min'], keep='last').dropna()

# parameters for plot/ figure
flowrate = 0.5
xmin = 7
xmax = 22
step = 1  # steps on x-axis
ymin = df['UV1_280nm'].min() + 5
#ymax = df['UV1_280nm'].max() + 5
ymax = 220

# define the different columns as plot and x values
# "elution" ist die von "min" in elution volumen umgerechnete spalte
# "frac" ist die spalte mit den gesammelten fraktionen in dem data sheet
elution = df['min'].dropna().apply(lambda x: x * flowrate)
frac = data.iloc[:, [20]].dropna().apply(lambda x: x * flowrate)

y1 = df['UV1_280nm']
# remove zero values from 360 (and 410nm)
# TODO: es müssen beide fälle gleichzeitig abgedeckt werden, dass es 0 werte bei 360 und 410 gibt!
df_nonull = df[df['UV3_410nm'] != 0]
# get min value for UV2_360 and subtract. baseline at 360 nm ~55
y2_min = df_nonull['UV2_360nm'].min()
y2 = df['UV2_360nm'].dropna().apply(lambda x: x - (y2_min + 65))
# get min value for UV3_410 and subtract. basline at 410 nm ~80
y3_min = df_nonull['UV3_410nm'].min()
y3 = df['UV3_410nm'].dropna().apply(lambda x: x - (y3_min + 90))

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title("Size Exclusion Chromatography (SEC)", fontsize=25, y=1.3)

# include second, third and furth y-axis sharing the same x-axis (twinx)
par1 = host.twinx()
# par2 = host.twinx()
par3 = host.twinx()

# set x and y axis + labels + define tick size
host.set_xlabel("Elution volume [ml]", fontsize=18)
host.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
host.tick_params(axis='y', labelsize=14)
par1.set_ylabel("Absorption [mAU]", fontsize=18, color='black')
par1.tick_params(axis='y', labelsize=14)
# par2.set_ylabel("Absorption at 410 nm [mAU]", fontsize=18)
par3.set_ylabel('Activity [%]', fontsize=18, color='grey')
par3.tick_params(axis='y', labelsize=14)

# add a second x-axis
host2 = host.twiny()
host2.set_xlabel("Collected fractions", fontsize=18, labelpad=10)

# add a third x-axis
host3 = host.twiny()
host3.set_xlabel("Molecular weight [kDa]", fontsize=18, labelpad=10)

host.set_xlim(xmin, xmax)
host.set_ylim(ymin, ymax)
host2.set_xlim(xmin, xmax)
host3.set_xlim(xmin, xmax)
par1.set_ylim(0, 100)
# par2.set_ylim(700, 1200)
# par3.set_ylim(ymin, ymax)

# set x1-ticks
host.xaxis.set_ticks(np.arange(xmin, xmax, step))
host.tick_params(axis='x', labelsize=14)

# set x2 ticks (collected fractions)
x3 = [str(i + 1) for i in range(frac.size)]  # zählt einfach nur durch die anzahl der fractions durch
host2.xaxis.set_ticks(frac['min.10'].values.tolist())  # wo die ticks gesetzt werden
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

# set host3 ticks (Molecular weight [kDa])
host3.xaxis.set_ticks(np.arange(xmin, xmax, step))
# only integer in range() accepted
array = np.arange(xmin, xmax, step)

# calculate molecular sizes according to regression from standard run (27.07.2021) y=-2.8799x + 6.986 !!! Ve=x in ml !!!!
# x4 = [round(10 ** (i * (-2.8799) + 6.986), 1) for i in
#      array]  # int function shows number without decimal places, round function defines decimal places
x4 = [int(10 ** (i * (-0.1894) + 4.6396)) for i in
      array]
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
p2, = par1.plot(xnew, power_smooth2, color=color2, label="360 nm")
p3, = par1.plot(xnew, power_smooth3, color=color3, label="410 nm")
par3.bar(list(map(lambda i: x2.tolist()[i - 1], frac)), act, 1, align='edge', color='grey', alpha=0.55,
         edgecolor="black")

# TODO: get a joined legend of host and par1
# set legend
#px = p1, + p2, + p3,
#labs = [l.get_label() for l in px]
#host.legend(px, labs, loc=10)

host.legend(loc=2)
par1.legend(loc=0)

# set grid
host.xaxis.grid(color='gray', linestyle='-', linewidth=0.5)
host.yaxis.grid(color='gray', linestyle='-', linewidth=0.5)

# define positions of scales / right, left, top, bottom
par3.spines['right'].set_position(('outward', 60))
host3.spines['top'].set_position(('outward', 50))

# Adjust spacings w.r.t. figsize
fig.tight_layout()

# host.yaxis.label.set_color(p1.get_color())
# par1.yaxis.label.set_color(p2.get_color())
# par2.yaxis.label.set_color(p3.get_color())

plt.show()
fig.savefig('SEC_20211124_DBT_Reactor.png')
