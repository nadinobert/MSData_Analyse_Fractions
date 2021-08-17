import sqlite3
import pandas as pd
from sqlalchemy import *
import matplotlib.pyplot as plt
import numpy as np

# open csv file of interest but skip the first two rows
data = pd.read_csv(r'20210803_ReactorEff_DBT.csv', skiprows=2, delimiter='\t')
df = pd.DataFrame(data, columns= ['ml', '280 nm', '360 nm', '410 nm'])
print(df)
#d = {'col1': [1, 2], 'col2': [3, 4]}
#df = pd.DataFrame(data=d)
plt.figure()
df.plot(x='ml', y="280 nm", secondary_y=['360 nm', '410 nm'], title='SEC 20210803', figsize=(15,10), fontsize=12, xticks=np.arange(0.0, 3.66, 0.1), xlim=[0.5, 2.5], xlabel='Elution volume', ylabel='Absorption [mAU]', grid=True)
#ax = df['360 nm'].plot(x='ml', y="280 nm",secondary_y=True)
plt.show()
#plt.savefig('nadinePlot.png')

#engine = create_engine('sqlite:///C:/Users/nadin/Documents/UFZ/ms_data.sqlite')
#conn = engine.connect()

#fetch all rows from created table
#print(conn.execute('''SELECT * FROM fraction_info''').fetchall())   #conn.execute() returns a sqlite3.Cursor object. Cursors can be thoght of as iterators
