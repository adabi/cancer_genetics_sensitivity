import pubchempy as pcp
import pandas as pd
from padelpy import padeldescriptor, from_smiles
import numpy as np


def find_smiles(drugs_file):
    drugs_df = pd.read_csv(drugs_file)
    drugs = list(drugs_df['drug_id'])
    smiles_lst = []

    # This list keeps track of the drugs for which a smiles was extracted successfully
    drugs_with_smiles = []

    # Find the isomeric smiles for all the drugs
    for drug in drugs:

        # Sometimes multiple Pubchem ids are listed and separated by a comma, just go with the first
        pubchem_id = drugs_df['pubchem'][drugs_df['drug_id'] == drug].tolist()[0].split(',')[0]

        try:
            # If the Pubchem id is an integer, use that
            pubchem_id = int(pubchem_id)
            smiles = pcp.Compound.from_cid(pubchem_id).isomeric_smiles
            smiles_lst.append(smiles)
            drugs_with_smiles.append(drug)

        except ValueError:
            # If it's something else like '-' or 'several' use the drug name to find the smiles
            drug_name = drugs_df['drug_name'][drugs_df['drug_id'] == drug].tolist()[0]

            # get_compounds returns a list of compounds, grab the first one
            compound = pcp.get_compounds(drug_name, 'name')[0]

            if compound:
                smiles = compound.isomeric_smiles
                smiles_lst.append(smiles)
                drugs_with_smiles.append(drug)

    with open('smiles.smi', 'w+') as file:
        for item in smiles_lst:
            file.write(f'{item}\n')

    with open('./Data/Clean/drugs_with_smiles.csv', 'w+') as file:
        for drug in drugs_with_smiles:
            file.write(f'{drug}\n')


def calculate_descriptors(drugs_file):
    find_smiles(drugs_file)
    # Find the descriptors from the smiles and store it
    padeldescriptor(
                    mol_dir='smiles.smi',
                    d_file='./Data/Clean/descriptors.csv',
                    convert3d=True,
                    retain3d=True,
                    d_2d=True,
                    d_3d=True,
                    fingerprints=False)

    # The descriptors are sometimes out of order, sort them so they match the drugs_with_smiles order
    descriptors_df = pd.read_csv('./Data/Clean/descriptors.csv')
    descriptors_df['Index'] = [int(x[2]) for x in descriptors_df.Name.str.split("_")]
    descriptors_df.set_index('Index', drop=True, inplace=True)
    descriptors_df.sort_index(inplace=True)

    # Write the correct order back to the file
    descriptors_df.to_csv('./Data/Clean/descriptors.csv', index=False)

    # Replace all the missing features and infinities with 0s (might have to change this later if
    # the model doesn't shake out)
    cols = (list(descriptors_df.columns))[1:]
    descriptors_df[cols] = descriptors_df[cols].replace({np.nan: 0, 'infinity': 0})
    descriptors_df.to_csv('./Data/Clean/descriptors_replaced.csv')


