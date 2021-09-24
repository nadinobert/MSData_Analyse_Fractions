import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline

# open csv file of interest but skip the first two rows
data = pd.read_csv(r'20210923_SEC.csv', skiprows=2, delimiter=';')  # falls mit tabs getrennt '\t'
df = pd.DataFrame(data, columns=['min', 'UV1_280nm', 'UV2_360nm', 'UV3_410nm'])
df = df.drop_duplicates(subset=['min'], keep='last').dropna()
# define the different columns as plot and x values
# x2 ist die spalte mit den gesammelten fraktionen in dem data sheet
#x1 = df['min']
x1 = df['min'].dropna().apply(lambda x: x * 0.02)
print(x1)
# x2 = data['fractions'].dropna().apply(lambda x: x * 0.7)
x2 = data.iloc[:, [20]].dropna()
print(x2)
print(type(x2))
y1 = df['UV1_280nm']
y2 = df['UV2_360nm']
y3 = df['UV3_410nm']
# act = data['activity'].dropna()
# print(act)
#frac = [10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]

xmin = 0.2
xmax = 2.5
step = 0.1
ymin = 0
ymax = 450

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title("Size Exclusion Chromatography (SEC)", fontsize=25, y=1.3)

# incorporate second, third and forth y-axis sharing the same x-axis (twinx)
par1 = host.twinx()
par2 = host.twinx()
par3 = host.twinx()

# set x and y axis + labels
host.set_xlabel("Elution volume [ml]", fontsize=18)
host.set_ylabel("Absorption at 280 nm [mAU]", fontsize=18, color='blue')
host.tick_params(axis='y', labelsize=14)
par1.set_ylabel("Absorption at 360 nm [mAU]", fontsize=18)
par2.set_ylabel("Absorption at 410 nm [mAU]", fontsize=18)
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
par3.set_ylim(ymin, ymax)
host3.set_xlim(xmin, xmax)
par1.set_ylim(300, 1200)
par2.set_ylim(300, 1200)

# set colors
color1 = 'blue'
color2 = 'orange'
color3 = 'grey'

# power_smooth einleiten um plot smooth zu machen und so
xnew = np.linspace(min(x1.tolist()), max(x1.tolist()), 912)

spl1 = make_interp_spline(x1.tolist(), y1, k=3)
power_smooth1 = spl1(xnew)
spl2 = make_interp_spline(x1.tolist(), y2, k=3)
power_smooth2 = spl2(xnew)
spl3 = make_interp_spline(x1.tolist(), y3, k=3)
power_smooth3 = spl3(xnew)

p1, = host.plot(xnew, power_smooth1, color=color1, label="280 nm")
p2, = host.plot(xnew, power_smooth2, color=color2, label="360 nm")
p3, = host.plot(xnew, power_smooth3, color=color3, label="410 nm")
# par3.bar(list(map(lambda i: x2.tolist()[i - 1], frac)), act, 1, align='edge', color='grey', alpha=0.55,
#         edgecolor="black")

# set legend
# lns = [p1]
# host.legend(handles=lns, loc='best')

# set x1-ticks
host.xaxis.set_ticks(np.arange(xmin, xmax, step))
host.tick_params(axis='x', labelsize=14)

# set x2 ticks (collected fractions)
x3 = [str(i + 1) for i in range(x2.size)]   # zählt einfach nur durch die anzahl der fractions durch
print(x2['min.10'].values.tolist())
host2.xaxis.set_ticks(x2['min.10'].values.tolist())  # wo die ticks gesetzt werden
#print(x2.tolist())
host2.xaxis.set_ticklabels(x3)  # was an den ticks steht
host2.tick_params(axis='x', labelsize=14)

# shift the x-ticks slightly more right to align labels with ticks
# create offset transform (x=25pt) -> shifting distance
from matplotlib.transforms import ScaledTranslation

dx, dy = 33, 0
offset = ScaledTranslation(dx / fig.dpi, dy / fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to all xticklabels
for label in host2.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

host2.set_xlim(xmin, xmax)
# set x3 ticks (Molecular weight [kDa])
# according to regression from standard run (27.07.2021) y=2.1239+1.9654*i-1.3452*i**2
host3.xaxis.set_ticks(np.arange(xmin, xmax, step))
# only integer in range() accepted
array = np.arange(xmin, xmax, step)

# x4 = [round((10**(i*(-2.8799)+6.986)), 1) for i in array]  # round to 1 digit after comma
x4 = [round(10 ** (2.1239 + 1.9654 * i / 9.6075 - 1.3452 * (i / 9.6075) ** 2), 1) for i in
      array]  # int function shows number without decimal places, round function defines decimal places
host3.xaxis.set_ticklabels(x4)
host3.tick_params(axis='x', labelsize=14)
# set grid
host.xaxis.grid(color='gray', linestyle='-', linewidth=0.5)
host.yaxis.grid(color='gray', linestyle='-', linewidth=0.5)

# define positions of scales / right, left, top, bottom
# par2.spines['right'].set_position(('outward', 60))
host3.spines['top'].set_position(('outward', 50))

# Adjust spacings w.r.t. figsize
fig.tight_layout()

# host.yaxis.label.set_color(p1.get_color())
# par1.yaxis.label.set_color(p2.get_color())
# par2.yaxis.label.set_color(p3.get_color())

plt.show()
fig.savefig('SEC_20210127_LB119.png')
