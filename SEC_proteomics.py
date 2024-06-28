# for plotting the chart
import matplotlib.pyplot as plt
# for creating a proxy artist for the legend
from matplotlib.transforms import ScaledTranslation
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from natsort import index_natsorted
from sqlalchemy import *

engine = create_engine('sqlite:///C:/Users/hellmold/Code/MSData_Analyse_Fractions/ms_data.sqlite')
conn = engine.connect()

# Set plot font globally to Times New Roman and set globally font size to 20
# locally set fontsize will overwrite globally defined setting!
plt.rcParams.update({'font.size': 20, 'font.family': 'Times New Roman'})

##TODO obacht filter für samples die BN im namen enthalten!!!
##TODO Anpassen: OmeB und RdhB sollen auch mit nur einem detektierten peptid angezeigt werden
##TODO labels die den subgroups zugeordnet werden überragen den plot! -> jetzt einfach hard gecoded... kein ahnung Zeile 156

# data contains the joined columns and rows from db tables "proteins", "result" and "analysis" where the requested date and protein descriptions match
# no contaminants
# Filter for non-transmembrane proteins in the OHR complex:
    # only proteins under the top ten ranked proteins
    # mind 2 peptide
    # min abundance 10^5
# no mind peptide or min abundance filter for RdhB and OmeB applied
data = pd.read_sql_query('''
        SELECT * FROM (
    SELECT *, RANK() OVER (
        PARTITION BY t.result_id
        ORDER BY t.abundance DESC
    ) AS rank FROM (
        SELECT
            r.sample,
            proteins.result_id,
            proteins.accession,
            proteins.description,
            proteins.coverage,
            proteins.numPeptides,
            proteins.abundance,
            proteins.MW,
            proteins.numUniquePeptides
        FROM proteins
        INNER JOIN result r ON r.id = proteins.result_id
        INNER JOIN analysis a ON a.id = r.analysis_id
        WHERE r.analysis_id = 94 -- Adjust as per your database schema
            AND proteins.accession LIKE 'cbdbA%'
            AND proteins.abundance <> ''
            AND (
                proteins.numUniquePeptides >= 2
                OR (
                    proteins.description LIKE '%RdhB%'
                    AND proteins.numUniquePeptides >= 1
                )
                OR (
                    proteins.description LIKE '%OmeB%'
                    AND proteins.numUniquePeptides >= 1
                )
            )
            AND (
                proteins.abundance > 100000
                OR proteins.description IN ('RdhA', 'OmeA', 'HupL', 'HupS', 'HupX') -- Removed RdhB and OmeB from here
                OR proteins.description LIKE '%RdhB%' -- Include RdhB if abundance condition is met
                OR proteins.description LIKE '%OmeB%' -- Include OmeB if abundance condition is met
            )
            AND r.sample NOT LIKE '%FT%'
            AND r.sample NOT LIKE '%frac%'
            AND r.sample NOT LIKE '%B%'
    ) AS t
) AS u
WHERE u.rank < 11
    AND (
        u.description LIKE '%rdhA%' -- AND u.accession = 'cbdbA0238'
        OR u.description LIKE '%rdhB%'
        OR u.description LIKE '%OmeA%'
        OR u.description LIKE '%OmeB%'
        OR u.description LIKE '%hupL%'
        OR u.description LIKE '%hupS%'
        OR u.description LIKE '%hupX%'
    ) 
''', conn)


print(data)

conn.close()

# in data the selected, joined entries get split in the field "description" to get the accessions
# new column 'colour" is created
# according to the split protein description in data, a colour for each entry is defined

# protein comA muss vorher rausgefiltert werden cbdbA0031
data = data[data['accession'] != 'cbdbA0031']

subgroup = data['description'].str.split(' ', expand=True)
#print(subgroup)

data['colour'] = [
    "#FA1912" if ele == 'rdhA'
    else "#FBCA0A" if ele == 'hupX'
    else "#07B0EF" if ele == 'hupS'
    else "#1A7FC4" if ele == "hupL"
    else "#08AF57" if ele == 'omeA'
    else "#A0d663" if ele == 'omeB'
    else "#F57AB1" if ele == 'rdhB'
    else " "
    for ele in subgroup[0]]

# 1. sorting data according to color (Reihenfolge der proteingruppen im finalen plot!!)
# 2. order subplots/ experiments according to sample name
# "sample" entries in data gets naturally ordered and the corresponding sorting is stored in new column "order"
data = data.sort_values('colour')
order = np.argsort(index_natsorted(data["sample"]))
data['order'] = order
data = data.set_index('order')
data.sort_index(inplace=True)

# create plot figure and axis to add on the bar charts
fig = plt.figure(figsize=(16, 8))  # default figsize=12,8
ax = fig.add_subplot()
fig.tight_layout()

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

    # use a for loop to fill x_ticks array with position number for accession description on x-axis
    # index + offset = numeric position
    for i in range(len(subgroup['accession'])):
        x_ticks.append(index + offset)
        index += 1

    # add the selected rows as a sub bar chart to the axis
    ax.bar(subgroup.index + offset, subgroup['abundance'], log=true, color=subgroup['colour'],
           label=subgroup['accession'], width=0.8)

    # add text to print the subgroup name above the sub chart
    ax.text((((subgroup.tail(1).index + subgroup.head(1).index) // 2)[0]) - 0.5 + offset -0.5, # the last value (offset- x ) defines the position of the label
            subgroup['abundance'].max() * 1.1,
            subgroup['sample'].unique()[0])
    # increment offset by one (after a loop over one experiment -> defined by sample name) so there is a space free between the groups
    offset += 3

# add protein accessions as label to the x-axis
ax.set_xticks(x_ticks)
ax.set_xticklabels(labels, rotation=70, ha='right', fontsize=15)
ax.set_ylim(100000, 100000000) # hard coded y-lim --> damit die scheiß labels nicht rausragen

# shift the x-ticks slightly more right to align labels with ticks
# create offset transform (x=7pt) -> shifting distance
from matplotlib.transforms import ScaledTranslation
dx, dy = 7, 0
offset = ScaledTranslation(dx/fig.dpi, dy/fig.dpi, scale_trans=fig.dpi_scale_trans)

# apply offset transform to all xticklabels
for label in ax.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# legend
# create a dataframe to store all unique protein groups with short descripton (e.g. 'hupX')
# with the corresponding color
legend = pd.DataFrame()
legend['des'] = data['description'].str.split(' ', expand=True)[0].unique()
legend['color'] = [
    "#FA1912" if ele == 'rdhA'
    else "#FBCA0A" if ele == 'hupX'
    else "#07B0EF" if ele == 'hupS'
    else "#1A7FC4" if ele == "hupL"
    else "#08AF57" if ele == 'omeA'
    else "#A0d663" if ele == 'omeB'
    else "#F57AB1" if ele == 'rdhB'
    else " "
    for ele in legend['des']]

legend = legend.sort_values('color')

# create a empty patch list (search for 'Proxy artists' for more information on creating legend
# https://matplotlib.org/stable/tutorials/intermediate/legend_guide.html#sphx-glr-tutorials-intermediate-legend-guide-py
patches = []
# iterate over the rows in the above created legend df
for i, row in legend.iterrows():
    # append a new patch with the color and the descriptor of the row element
    patches.append(mpatches.Patch(color=row['color'], label=(row['des'][0].upper() + row['des'][1:])))

# create a legend based on the patch list
plt.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)

# add some space below the plot and to the right
plt.subplots_adjust(bottom=0.2, right=0.9, top=0.9)

# name axis
ax.set_ylabel(ylabel='Intensity')

# name the plot
# name according the "comment" entry from "analysis" extracted previously when sql connection was established
#plt.title('Proteomic Analysis of Collected Fractions', fontsize=25, y=1.05)

# Adjust spacings w.r.t. figsize
fig.tight_layout()

# show the plot or save it as a .png
plt.show()
fig.savefig('20240401_Ingel')