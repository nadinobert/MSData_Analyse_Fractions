# for plotting the chart
import matplotlib.pyplot as plt
# for creating a proxy artist for the legend
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from natsort import index_natsorted
from sqlalchemy import *

engine = create_engine('sqlite:///C:/Users/hellmold/Nextcloud/ms_data.sqlite')
conn = engine.connect()

# data contains the joined columns and rows from db tables "proteins", "result" and "analysis" where the requested date and protein descriptions match
data = pd.read_sql_query('''SELECT * FROM
(SELECT proteins.accession, CASE WHEN proteins.abundance <> '' THEN proteins.abundance ELSE 0 END AS abundance, result.sample, proteins.description
FROM proteins
inner join result on result.id = proteins.result_id
inner join analysis on analysis.id = result.analysis_id
where analysis.id = 29 AND 
(description LIKE '%rdhA%'
   OR description LIKE '%rdhB%'
   OR description LIKE '%OmeA%'
   OR description LIKE '%OmeB%'
   OR description LIKE '%hupL%'
   OR description LIKE '%hupX%'
   OR description LIKE '%hupS%'))
WHERE abundance <> 0
;''', conn)

# name the plot
# plot should be named after "comment"- entry in "analysis" table
titel = pd.read_sql_query('''
SELECT * FROM analysis
WHERE analysis.id = 29
;''', conn)

plotname = titel['comment'].iloc[0]
print(plotname)

conn.close()

# in data the selected, joined entries get split in the field "description" to get the accessions
# new column 'colour" is created
# according to the split protein description in data, a colour for each entry is defined
subgroup = data['description'].str.split(' ', expand=True)
data['colour'] = [
    "#FA1912" if ele == 'rdhA'
    else "#FBCA0A" if ele == 'hupX'
    else "#07B0EF" if ele == 'hupS'
    else "#1A7FC4" if ele == "hupL"
    else "#08AF57" if ele == 'omeA'
    else "g" if ele == 'omeB'
    else "#F57AB1" if ele == 'rdhB'
    else " "
    for ele in subgroup[0]]

# order subplots/ experiments according to sample name
# "sample" entries in data gets naturally ordered and the corresponding sorting is stored in new column "order"
order = np.argsort(index_natsorted(data["sample"]))
data['order'] = order
data = data.set_index('order')
print(data)
data.sort_index(inplace=True)
print(data)

#pd.set_option("display.max_rows", None, "display.max_columns", None)
#print(data)

# create plot figure and axis to add on the bar charts
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot()

# iterate through the unique values of 'sample' in data
# -> .loc like index calling
# -> ax. adds a subplot different to plt
offset = 0
index = 0
x_ticks = []
labels = pd.Series()
for experiment in data['sample'].unique():
    subgroup = data.loc[data['sample'] == experiment]
    labels = labels.append(subgroup['accession'])

    # use a for loop to fill x_ticks array with positions for accession description on x-axis
    # index + offset = numeric position
    for i in range(len(subgroup['accession'])):
        x_ticks.append(index + offset)
        index += 1

    # add the selected rows as a sub bar chart to the axis
    ax.bar(subgroup.index + offset, subgroup['abundance'], log=true, color=subgroup['colour'], label=subgroup['accession'])

    # add text to print the subgroup name above the sub chart
    ax.text((((subgroup.tail(1).index + subgroup.head(1).index) // 2)[0]) - 0.5 + offset, subgroup['abundance'].max() * 1.1,
            subgroup['sample'].unique()[0])

    # increment offset by one (after a loop over one experiment -> defined by sample name) so there is a space free between the groups
    offset += 1

# add protein accessions as label to the x-axis
ax.set_xticks(x_ticks)
ax.set_xticklabels(labels, rotation=70, ha='right')

# shift the x-ticks slightly more right to align labels with ticks
# create offset transform (x=7pt) -> shifting distance
from matplotlib.transforms import ScaledTranslation
dx, dy = 7, 0
offset = ScaledTranslation(dx/fig.dpi, dy/fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to all xticklabels
for label in ax.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# legend
# create a dataframe to store all unique protein groups with short descriptor (e.g. 'hupX')
# with the corresponding color
legend = pd.DataFrame()
legend['des'] = data['description'].str.split(' ', expand=True)[0].unique()
legend['color'] = [
    "#FA1912" if ele == 'rdhA'
    else "#FBCA0A" if ele == 'hupX'
    else "#07B0EF" if ele == 'hupS'
    else "#1A7FC4" if ele == "hupL"
    else "#08AF57" if ele == 'omeA'
    else "g" if ele == 'omeB'
    else "#F57AB1" if ele == 'rdhB'
    else " "
    for ele in legend['des']]

# create a empty patch list (search for 'Proxy artists' for more information on creating legend
# https://matplotlib.org/stable/tutorials/intermediate/legend_guide.html#sphx-glr-tutorials-intermediate-legend-guide-py
patches = []
# iterate over the rows in the above created legend df
for i, row in legend.iterrows():
    # append a new patch with the color and the descriptor of the row element
    patches.append(mpatches.Patch(color=row['color'], label=row['des']))

# create a legend based on the patch list
plt.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)

# add some space below the plot and to the right
plt.subplots_adjust(bottom=0.2, right=0.9)

# name axis
ax.set_ylabel('intensity')

# name the plot
# name according the "comment" entry from "analysis" extracted previously when sql connection was established
plt.title(plotname)

# show the plot or save it as a .png
plt.show()
#plt.savefig('nadinePlot2.png')