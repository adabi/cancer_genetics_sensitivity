import pubchempy as pcp
import pandas as pd
from padelpy import padeldescriptor

drugs_df = pd.read_csv('Data/Raw/Drugs.csv')
drugs = list(drugs_df['drug_id'])
smiles_lst = []

# This list keeps track of the drugs for which a smiles was extracted successfully
drugs_with_smiles = []

# Find the isomeric smiles for all the drugs
# for drug in drugs:
#
#     # Sometimes multiple Pubchem ids are listed and separated by a comma, just go with the first
#     pubchem_id = drugs_df['pubchem'][drugs_df['drug_id'] == drug].tolist()[0].split(',')[0]
#
#     try:
#         # If the Pubchem id is an integer, use that
#         pubchem_id = int(pubchem_id)
#         smiles = pcp.Compound.from_cid(pubchem_id).isomeric_smiles
#         smiles_lst.append(smiles)
#         drugs_with_smiles.append(drug)
#
#     except ValueError:
#         # If it's something else like '-' or 'several' use the drug name to find the smiles
#         drug_name = drugs_df['drug_name'][drugs_df['drug_id'] == drug].tolist()[0]
#
#         # get_compounds returns a list of compounds, grab the first one
#         compound = pcp.get_compounds(drug_name, 'name')
#
#         if compound:
#             smiles = compound[0].isomeric_smiles
#             smiles_lst.append(smiles)
#             drugs_with_smiles.append(drug)
#
# with open('smiles.smi', 'w+') as file:
#     for item in smiles_lst:
#         file.write(f'{item}\n')
#
# with open('drugs_with_smiles.csv', 'w+') as file:
#     for drug in drugs_with_smiles:
#         file.write(f'{drug}\n')

padeldescriptor(
                mol_dir='smiles.smi',
                d_file='descriptors.csv',
                convert3d=True,
                retain3d=True,
                d_2d=True,
                d_3d=True,
                fingerprints=False,
            )