import sqlite3
import pandas as pd
from sqlalchemy import *
import matplotlib.pyplot as plt
import numpy as np

# open csv file of interest but skip the first two rows
data = pd.read_csv(r'20210803_ReactorEff_DBT.csv', skiprows=2, delimiter='\t')
df = pd.DataFrame(data, columns=['ml', '280 nm', '360 nm', '410 nm'])

# define the different columns as plot and x values
x1 = pd.DataFrame(data, columns=['ml'])
x2 = data.iloc[:,[20]].dropna()
y1 = pd.DataFrame(data, columns=['280 nm'])
y2 = pd.DataFrame(data, columns=['360 nm'])
y3 = pd.DataFrame(data, columns=['410 nm'])

# Create figure and subplot manually
fig = plt.figure()
fig.set_size_inches(12, 8)
host = fig.add_subplot(111)
plt.title("SEC 20210803")

par1 = host.twinx()
par2 = host.twinx()

host.set_xlim(0.5, 2.5)
host.set_ylim(0, 1100)
par1.set_ylim(300, 1100)
par2.set_ylim(300, 1100)

# set x and y axis + labels
host.set_xlabel("Elution volume [ml]")
host.set_ylabel("Absorption at 280 nm [mAU]")
par1.set_ylabel("Absorption at 360 nm [mAU]")
par2.set_ylabel("Absorption at 410 nm [mAU]")
# add a second x-axis
host2 = host.secondary_xaxis('top')
host2.set_xlabel("Collected fractions")

color1 = plt.cm.viridis(0)
color2 = plt.cm.viridis(0.5)
color3 = plt.cm.viridis(.9)

p1, = host.plot(x1, y1, color=color1, label="280 nm")
p2, = par1.plot(x1, y2, color=color2, label="360 nm")
p3, = par2.plot(x1, y3, color=color3, label="410 nm")

#set legend
lns = [p1, p2, p3]
host.legend(handles=lns, loc='best')

# set x1-ticks
host.xaxis.set_ticks(np.arange(0.5, 2.5, step=0.1))

# set x2 ticks
x3 = ["frac" + str(i + 1) for i in range(x2.size)]
print(x2["ml.10"].tolist())
host2.xaxis.set_ticks(x2["ml.10"].tolist())
host2.xaxis.set_ticklabels(x3)

# set grid
host.xaxis.grid(color='gray', linestyle='-', linewidth=0.5)
host.yaxis.grid(color='gray', linestyle='-', linewidth=0.5)

# right, left, top, bottom
par2.spines['right'].set_position(('outward', 60))

# Adjust spacings w.r.t. figsize
fig.tight_layout()

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())

plt.show()
