import sqlite3

from sqlalchemy import create_engine

engine = create_engine('sqlite:///C:/Users/nadin/Nextcloud/ms_data.sqlite')
conn = engine.connect()
c = engine.cursor()

with open('BL21_seperated.txt', 'r') as TM:
    lines = TM.readlines()

for line in lines:
    data = line.split()
    acc = data[0]
    st = data[1]
    en = data[2]

c.execute("""INSERT INTO BL21_transmembrane(Accession, Start, End) VALUES(%s, %s, %s) """, (acc, st, en))

c.commit()
c.close()