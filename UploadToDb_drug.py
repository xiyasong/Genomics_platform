#!/usr/bin/env python
import pickle, sys, io, base64, os
import pandas as pd
import numpy as np
import sqlalchemy
from sqlalchemy.types import VARCHAR

database_username = 'szareporting'
database_password = 'SZA123reporting!'
#database_ip   = '178.157.14.88'
database_ip   = '16.171.6.238'
#variant_report_db = sqlalchemy.create_engine('mysql+mysqlconnector://{0}:{1}@{2}/{3}'.format(database_username, database_password, database_ip, 'szareporting_variant_report_genomics_v3'))
#patientdata_db = sqlalchemy.create_engine('mysql+mysqlconnector://{0}:{1}@{2}/{3}'.format(database_username, database_password, database_ip, 'szareporting_patientdata'))
drug_db = sqlalchemy.create_engine('mysql+mysqlconnector://{0}:{1}@{2}/{3}'.format(database_username, database_password, database_ip, 'szareporting_patientdata'))
select_group = "123456"
target_groups = {
	'1': 'Basic (for healthy subjects): Usually used for Health predipositions/Disease risk',
	'2': 'Basic (for healthy subjects): Usually used for Carrier-screening',
	'3': 'Basic (for healthy subjects): Usually used for Heriditary-cancer risk syndrome',
	'4': 'Basic (for healthy subjects): Usually used for Newborn-screening',
	'5': 'Extended (for potential patients): Expanded_mono_rare_diseases',
	'6': 'Extended (for potential patients): Expanded_clinvar_mono_diseases'
}
select_group_name = [target_groups[i] for i in list(select_group)]

pat = sys.argv[1].split('.')[0]+ '.drug'
print(pat)

inheritance = pd.read_csv(sys.argv[1], sep = '\t')
#inheritance['Target Class'] = inheritance['final_target_group'].str.replace(',', ': ')
#inheritance['REPORT_TAG'] = 1
#del inheritance['Target.group']
#del inheritance['Database_type']
# inheritance.drop_duplicates().reset_index(drop = True).to_sql(con=variant_report_db, name=pat,if_exists='replace')

# patientdata = pd.DataFrame()
# patientdata.loc['Name', ['Patient Info', 'TYPE']] = ['XX', 'SELF']
# patientdata.loc['DOB', ['Patient Info', 'TYPE']] = ['2005-02-25', 'SELF']
inheritance.reset_index().to_sql(con=drug_db, name=pat,if_exists='replace')
