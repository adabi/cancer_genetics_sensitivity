import pubchempy as pcp
import pandas as pd
from padelpy import padeldescriptor, from_smiles
import numpy as np
from PIL import Image
from collections import OrderedDict
import scikit_posthocs
import os
import contextlib


def find_smiles(drugs_file):
    drugs_df = pd.read_csv(drugs_file)
    drugs = list(drugs_df['drug_id'])
    smiles_lst = []
    cid_lst = []
    # This list keeps track of the drugs for which a smiles was extracted successfully
    drugs_with_smiles = []

    # Find the isomeric smiles for all the drugs
    for drug in drugs:

        # Sometimes multiple Pubchem ids are listed and separated by a comma, just go with the first
        pubchem_id = drugs_df['pubchem'][drugs_df['drug_id'] == drug].tolist()[0].split(',')[0]

        try:
            # If the Pubchem id is an integer, use that
            pubchem_id = int(pubchem_id)
            smiles = pcp.Compound.from_cid(pubchem_id)
            smiles_lst.append(smiles)
            drugs_with_smiles.append(drug)
            cid_lst.append(pubchem_id)

        except ValueError:
            # If it's something else like '-' or 'several' use the drug name to find the smiles
            drug_name = drugs_df['drug_name'][drugs_df['drug_id'] == drug].tolist()[0]
            compounds = pcp.get_compounds(drug_name, 'name')
            if compounds:
                # get_compounds returns a list of compounds, grab the first one
                compound = compounds[0]
                smiles = compound.isomeric_smiles
                smiles_lst.append(smiles)
                drugs_with_smiles.append(drug)
                cid_lst.append(compound.cid)

    with open('smiles.smi', 'w+') as file:
        for item in smiles_lst:
            file.write(f'{item}\n')

    drugs_with_smiles_dict = {'drug_id': drugs_with_smiles, 'cid': cid_lst}
    drugs_with_smiles_df = pd.DataFrame(data=drugs_with_smiles_dict)
    drugs_with_smiles_df.to_csv(path_or_buf='./Data/Clean/drugs_with_smiles.csv')


def calculate_descriptors(drugs_file):
    print("Calculating descriptors..")
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
    # Replace nan and infinity values with 0s
    descriptors_df[cols] = descriptors_df[cols].replace({np.nan: 0, '': 0, 'Infinity': 0, '-Infinity': 0})
    descriptors_df = descriptors_df.apply(lambda x: np.array(x).astype(float))
    # Some values are huge which causes issues during numpy calculations
    # Make them not so enormous
    descriptors_df[descriptors_df > 10E100] = 10E100
    # Get rid of outliers using the Extreme Studentized Deviate test
    descriptors_df = descriptors_df.apply(scikit_posthocs.outliers_gesd)
    # Get rid of columns which have all 0s
    descriptors_df = descriptors_df.loc[:, (descriptors_df != 0).any(axis=0)]
    # Scale every column
    descriptors_df = descriptors_df.apply(scale_array)
    descriptors_df.to_csv('./Data/Clean/descriptors_scaled.csv', index=False)
    with contextlib.suppress(FileNotFoundError):
        os.remove('smiles.smi')


# Grab and store the pixel data for each drug
def calculate_drug_pixel_data(drugs_with_smiles_file):
    print("Calculating pixel data...")
    drugs_with_smiles_df = pd.read_csv(drugs_with_smiles_file)
    cid_lst = list(drugs_with_smiles_df['cid'])
    pixels_dict = OrderedDict()
    for i in range(3600):
        pixels_dict[f'pixel{i}'] = []
    for cid in cid_lst:
        # Download the picture of the compound from PubChem
        pcp.download('PNG', 'drug.png', int(cid), 'cid', overwrite=True)
        # Convert to single-channel greyscale
        img = Image.open('drug.png').convert('L')
        # Get the pixel data as a numpy array
        pixels = np.array(img)
        # The background for these images is grey and not white
        # Turn all grey pixels into white pixels
        pixels[pixels == 245] = 255
        # Make any non-grey pixel completely black
        # This ensures that all atoms and bonds have the same pixel intensity
        pixels[pixels < 245] = 0

        img = Image.fromarray(pixels)
        # Downsample using antialiasing to 60 by 60 pixels
        img = img.resize((60, 60), Image.ANTIALIAS)
        # Grab pixel data again
        pixels = np.array(img)
        # Reverse the pixels so that the darkest areas have the highest value
        pixels = 255-pixels
        # Flatten
        pixels = pixels.flatten()
        # Scale
        #pixels = scale_array(pixels)
        for i, pixel in enumerate(pixels):
            pixels_dict[f'pixel{i}'].append(pixel)

    drug_id_lst = list(drugs_with_smiles_df['drug_id'])
    pixels_dict['cid'] = cid_lst
    pixels_dict.move_to_end('cid', last=False)
    pixels_dict['drug_id'] = drug_id_lst
    pixels_dict.move_to_end('drug_id', last=False)
    drug_pixel_df = pd.DataFrame(data=pixels_dict)
    drug_pixel_df.to_csv('./Data/Clean/drug_pixels.csv', index=False)
    with contextlib.suppress(FileNotFoundError):
        os.remove('drug.png')


# This function grabs descriptor and drug data and combines it all
def calculate_drug_data(raw_drugs_file):
    calculate_descriptors(drugs_file=raw_drugs_file)
    calculate_drug_pixel_data(
       drugs_with_smiles_file='./Data/Clean/drugs_with_smiles.csv')
    print("Combining all data...")
    descriptors_df = pd.read_csv('./Data/Clean/descriptors_scaled.csv')
    pixels_df = pd.read_csv('./Data/Clean/drug_pixels.csv')
    total_df = pd.concat((descriptors_df, pixels_df), axis=1)
    cols = total_df.columns.tolist()
    cols.insert(0, cols.pop(cols.index('cid')))
    cols.insert(0, cols.pop(cols.index('drug_id')))
    total_df = total_df.reindex(columns=cols)
    total_df.to_csv('./Data/Clean/drugs_fulldata.csv', index=False)


def scale_array(array):
    scaled_array = (array - np.min(array))/(np.max(array) - np.min(array))
    return scaled_array

descriptors_df = pd.read_csv('./Data/Clean/descriptors.csv')
descriptors_df['Index'] = [int(x[2]) for x in descriptors_df.Name.str.split("_")]
descriptors_df.set_index('Index', drop=True, inplace=True)
descriptors_df.sort_index(inplace=True)
descriptors_df = descriptors_df.iloc[:, 0:1172]
print(descriptors_df.shape)
pixels_df = pd.read_csv('./Data/Clean/drug_pixels.csv')
descriptors_df['drug_id'] = pixels_df['drug_id']
# Replace all the missing features and infinities with 0s (might have to change this later if
# the model doesn't shake out)
cols = (list(descriptors_df.columns))[1:]
# Replace nan and infinity values with 0s
descriptors_df[cols] = descriptors_df[cols].replace({np.nan: 0, '': 0, 'Infinity': 0, '-Infinity': 0})
descriptors_df = descriptors_df.loc[:, (descriptors_df != 0).any(axis=0)]
descriptors_df = descriptors_df.drop(columns=['Name'])
descriptors_df.to_csv('./Data/Clean/descriptors_unscaled.csv', index=False)
