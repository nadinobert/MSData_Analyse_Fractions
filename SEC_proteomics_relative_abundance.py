# for plotting the chart
import matplotlib.pyplot as plt
# for creating a proxy artist for the legend
import matplotlib.patches as mpatches
import numpy as np
import output as output
import pandas as pd
from natsort import index_natsorted
from sqlalchemy import *

engine = create_engine('sqlite:///C:/Users/hellmold/Code/MSData_Analyse_Fractions/ms_data.sqlite')
conn = engine.connect()


# data contains the joined columns and rows from db tables "proteins", "result" and "analysis" where the requested date and protein descriptions match
data = pd.read_sql_query('''
    SELECT * FROM (
        SELECT 
            proteins.accession, 
            CASE WHEN proteins.abundance <> '' THEN CAST(proteins.abundance AS FLOAT) ELSE 0 END AS abundance, 
            result.sample, 
            proteins.description, 
            proteins.numPeptides
        FROM proteins
        INNER JOIN result ON result.id = proteins.result_id
        INNER JOIN analysis ON analysis.id = result.analysis_id
        WHERE 
            analysis.id = 93 AND
            (
                (description LIKE '%rdhA%' AND proteins.accession = 'cbdbA0238') OR
                description LIKE '%rdhB%' OR
                description LIKE '%OmeA%' OR
                description LIKE '%OmeB%' OR
                description LIKE '%hupL%' OR
                description LIKE '%hupS%' OR
                description LIKE '%hupX%'
            )
    ) AS subquery
    WHERE abundance <> 0
''', conn)


# name the plot
# plot should be named after "comment"- entry in "analysis" table
titel = pd.read_sql_query('''
SELECT * FROM analysis
WHERE analysis.id = 93
;''', conn)

conn.close()

print(data)

# in data the selected, joined entries get split in the field "description" to get the accessions
# new column 'colour" is created
# according to the split protein description in data, a colour for each entry is defined

# protein comA muss vorher rausgefiltert werden cbdbA0031
data = data[data['accession'] != 'cbdbA0031']

# filter für y scale range -> only proteins >10^4 intensity
data = data[data['abundance'] > 10000]
#filter for proteins with min 2 peptides detected
data = data[data['numPeptides'] > 0]

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

## For the relative abundance plot:
# Group by accession and calculate the total abundance for each accession number across all samples
accession_total_abundance = data.groupby("accession")["abundance"].sum().reset_index()

# Create a new DataFrame with distinct protein accession numbers as the first column
max_sample = data["sample"].max()
new_df = pd.DataFrame(columns=["accession"] + [str(i).zfill(2) for i in range(1, 50)])

# Iterate through each row in accession_total_abundance
for index, row in accession_total_abundance.iterrows():
    accession = row["accession"]
    total_abundance = row["abundance"]
    # iterate through the original dataframe
    protein_data = {}
    for sample in range(1, 50):
        sample_abundance = \
        data[(data["accession"] == accession) & (data["sample"].astype(str).str.zfill(2) == str(sample).zfill(2))][
            "abundance"]
        if not sample_abundance.empty:
            protein_data[str(sample).zfill(2)] = sample_abundance.values[0] / total_abundance
        else:
            protein_data[str(sample).zfill(2)] = 0
    new_row = pd.DataFrame([[accession] + list(protein_data.values())], columns=new_df.columns)
    new_df = new_df.append(new_row, ignore_index=True)

print(new_df)

# Plotting relative abundance
plt.figure(figsize=(10, 6))
plotname = titel['comment'].iloc[0]
# Iterate through each row in the new DataFrame
for index, row in new_df.iterrows():
    accession = row["accession"]
    abundances = row[1:]  # Abundances start from the second column
    colour = data.loc[data["accession"] == accession, "colour"].iloc[0]
    plt.plot(range(1, 50), abundances, label=accession, color=colour)

plt.xlabel('Sample')
plt.ylabel('Abundance Proportion')
plt.title('Abundance Proportion of Proteins across Samples')
plt.legend()
plt.grid(True)
plt.show()

## Plot for the protein per sample abundance
# 1. sorting data according to color (Reihenfolge der proteingruppen im finalen plot!!)
# 2. order subplots/ experiments according to sample name
# "sample" entries in data gets naturally ordered and the corresponding sorting is stored in new column "order"
data = data.sort_values('colour')
order = np.argsort(index_natsorted(data["sample"]))
data['order'] = order
data = data.set_index('order')
data.sort_index(inplace=True)

# create plot figure and axis to add on the bar charts
fig = plt.figure(figsize=(20, 8))   #default figsize=12,8
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

    # use a for loop to fill x_ticks array with position number for accession description on x-axis
    # index + offset = numeric position
    for i in range(len(subgroup['accession'])):
        x_ticks.append(index + offset)
        index += 1

    # add the selected rows as a sub bar chart to the axis
    ax.bar(subgroup.index + offset, subgroup['abundance'], log=true, color=subgroup['colour'], label=subgroup['accession'])

    # add text to print the subgroup name above the sub chart
    ax.text((((subgroup.tail(1).index + subgroup.head(1).index) // 2)[0]) - 0.5 + offset, subgroup['abundance'].max() * 1.1,
            subgroup['sample'].unique()[0], fontsize=14)
    # increment offset by one (after a loop over one experiment -> defined by sample name) so there is a space free between the groups
    offset += 1

# add protein accessions as label to the x-axis
ax.set_xticks(x_ticks)
ax.set_xticklabels(labels, rotation=70, ha='right', fontsize=9)     # default fontsize:14

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
plt.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., fontsize=16)

# add some space below the plot and to the right
plt.subplots_adjust(bottom=0.2, right=0.9)

# name axis
ax.set_ylabel(ylabel='Intensity', fontsize=18)

# name the plot
# name according the "comment" entry from "analysis" extracted previously when sql connection was established
#plt.title('Proteomic Analysis of Collected Fractions', fontsize=25, y=1.05)

# Adjust spacings w.r.t. figsize
fig.tight_layout()

# show the plot or save it as a .png
plt.show()
fig.savefig('20231218_Insolution')