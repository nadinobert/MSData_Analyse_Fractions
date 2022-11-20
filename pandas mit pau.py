import matplotlib.pyplot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import float64
from sqlalchemy import *

engine = create_engine('sqlite:///C:/Users/nadin/Documents/UFZ/ms_data.sqlite')
conn = engine.connect()

data = pd.read_sql_query('''SELECT proteins.accession, proteins.abundance, result.sample, proteins.description
FROM proteins
inner join result on result.id = proteins.result_id
inner join analysis on analysis.id = result.analysis_id
where date = '2021-04-23' AND
(description LIKE '%rdhA%'
   OR description LIKE '%rdhB%'
   OR description LIKE '%OmeA%'
   OR description LIKE '%OmeB%'
   OR description LIKE '%hupL%'
   OR description LIKE '%hupS%'
   OR description LIKE '%hupX%')
;''', conn)

temp = data['description'].str.split(' ', expand=True)
print(temp[0])
data['colour'] = [
    "#FA1912" if ele == 'rdhA'
    else "#FBCA0A" if ele == 'hupX'
    else "#07B0EF" if ele == 'hupS'
    else "#1A7FC4" if ele == "hupL"
    else "#08AF57" if ele == 'omeA'
    else "g" if ele == 'omeB'
    else "#F57AB1" if ele == 'rdhB'
    else " "
    for ele in temp[0]]
#data.plot('accession', 'abundance', kind='bar', log=true, color=data['colour'])
data = data.replace(to_replace='None', value=np.nan).dropna()   #remove the Non entries
print(data.info())
data['abundance'] = pd.to_numeric(data['abundance'])
data['abundance'] = data['abundance'].astype(str).astype(float64)
data.plot('accession', 'abundance', kind='bar', log=true)
plt.show()
#print(data.iloc[:, 0])  #print first column from dataframe "data"
#print(data.iloc[:, 1])  #print seconf column
#print(data.dtypes)      #error no numeric data to plot -> datatypes seems to be objects instead of numbers
#print(data.info())

#datagrouped = data.groupby('sample')
#group1 = datagrouped.get_group('DDMNaCLUN')
#print(datagrouped.get_group)
#group2 = datagrouped.get_group('DDM')
#print(group1.info())
#plt.bar(group1.index, group1['accession'], ['abundance'], log=true, color=data['colour'])
#plt.bar(group2.index + 1, group2['abundance'], log=true, color=data['colour'])
#plt.show()
