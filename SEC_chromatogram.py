import sqlite3
import pandas as pd
from sqlalchemy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# open csv file of interest but skip the first two rows
data = pd.read_csv(r'20210803_ReactorEff_DBT.csv', skiprows=2, delimiter='\t')
df = pd.DataFrame(data, columns=['ml', '280 nm', '360 nm', '410 nm'])

# define the different columns as plot and x values
# x2 ist die spalte mit den gesammelten fraktionen in dem data sheet
x1 = pd.DataFrame(data, columns=['ml'])
x2 = data.iloc[:,[20]].dropna()
y1 = pd.DataFrame(data, columns=['280 nm'])
y2 = pd.DataFrame(data, columns=['360 nm'])
y3 = pd.DataFrame(data, columns=['410 nm'])

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title("Size exclusion chromatography (SEC)", fontsize=25, y=1.3)

# incorporate second and third y-axis sharing the same x-axis (twinx)
par1 = host.twinx()
par2 = host.twinx()

# set x and y axis + labels
host.set_xlabel("Elution volume [ml]", fontsize=18)
host.set_ylabel("Absorption at 280 nm [mAU]", fontsize=18)
par1.set_ylabel("Absorption at 360 nm [mAU]", fontsize=18)
par2.set_ylabel("Absorption at 410 nm [mAU]", fontsize=18)
# add a second x-axis
host2 = host.secondary_xaxis('top')
host2.set_xlabel("Collected fractions", fontsize=18)
# add a third x-axis
host3 = host.twiny()
host3.set_xlabel("Molecular weight [kDa]", fontsize=18)

host.set_xlim(1.3, 2.5)
host.set_ylim(0, 1100)
host3.set_xlim(1.3, 2.5)
par1.set_ylim(300, 1100)
par2.set_ylim(300, 1100)

# set colors
color1 = 'blue'
color2 = 'orange'
color3 = 'grey'

#cs = interp1d(x1, y1)
p1, = host.plot(x1, y1, color=color1, label="280 nm")
p2, = par1.plot(x1, y2, color=color2, label="360 nm")
p3, = par2.plot(x1, y3, color=color3, label="410 nm")

#set legend
#lns = [p1, p2, p3]
#host.legend(handles=lns, loc='best')

# set x1-ticks
host.xaxis.set_ticks(np.arange(1.3, 2.5, step=0.1))

# set x2 ticks (collected fractions)
x3 = [str(i + 1) for i in range(x2.size)]
#print(x2["ml.10"].tolist())
host2.xaxis.set_ticks(x2["ml.10"].tolist())     # wo die ticks gesetzt werden
host2.xaxis.set_ticklabels(x3)  # was an den ticks steht

# shift the x-ticks slightly more right to align labels with ticks
# create offset transform (x=25pt) -> shifting distance
from matplotlib.transforms import ScaledTranslation
dx, dy = 25, 0
offset = ScaledTranslation(dx/fig.dpi, dy/fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to all xticklabels
for label in host2.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# set x3 ticks (Molecular weight [kDa])
# according to regression from standard run (27.07.2021)
#y = -2.8799x+6.986
host3.xaxis.set_ticks(np.arange(1.3, 2.5, step=0.1))
#x4 = [ (10**(i*(-2.8799)+6.986)) for i in range(1, 13)]
# only integer in range() accepted
array = np.arange(1.3, 2.5, 0.1)
print(array)
x4 = [round((10**(i*(-2.8799)+6.986)), 1) for i in array]  # round to 1 digit after comma
host3.xaxis.set_ticklabels(x4)
#print(x4)

# set grid
host.xaxis.grid(color='gray', linestyle='-', linewidth=0.5)
host.yaxis.grid(color='gray', linestyle='-', linewidth=0.5)

# define positions of scales / right, left, top, bottom
par2.spines['right'].set_position(('outward', 60))
host3.spines['top'].set_position(('outward', 60))

# Adjust spacings w.r.t. figsize
fig.tight_layout()

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())

plt.show()
#fig.savefig('SEC_20210803_Chromatogram.png')