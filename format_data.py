import pandas as pd
from collections import Counter, defaultdict
import numpy as np
from drug_data import calculate_drug_data

def clean_genetic_features():
    # The csv file is strangely formatted, so we need to clean it up for use with pandas
    with open("Data/Raw/genetic_features_all.csv") as file:
        buff = []
        for line in file:
            line_splt = line.split(',')
            if len(line_splt) > 9:
                # This combines the genes_in_segment into a single dash-delimited column
                # instead of that variable column number monstrosity
                modified_line = (line_splt[0:8] + [' - '.join(line_splt[8:])])
                buff.append(','.join(modified_line))
            else:
                buff.append(','.join(line_splt))

        with open("Data/Clean/genetic_features_clean.csv", 'w+') as write_file:
            write_file.write(''.join(map(str, buff)))


def find_useful_mutations():
    """This function finds mutations which show up as either sequence variation or number variation
    from the cancer_genes.csv file and stores them in cancer_genes_relevant.csv
    I am not using these now, but they may come in handy if the model does not manage to converge"""

    genetic_features_df = pd.read_csv('Data/Clean/genetic_features_clean.csv',
                                      dtype={'recurrent_gain_loss': str, 'genes_in_segment': str})
    cancer_genes_df = pd.read_csv('Data/Raw/cancer_genes.csv')

    seqvar_count_lst = []
    seqvar_count_dict = Counter(genetic_features_df['genetic_feature'])

    for gene in cancer_genes_df['gene_symbol']:
        mutations = 0
        if f'{gene}_mut' in seqvar_count_dict:
            mutations += seqvar_count_dict[f'{gene}_mut']
        seqvar_count_lst.append(mutations)

    numb_count_lst = []
    numb_count_dict = defaultdict(int)
    for row in genetic_features_df['genes_in_segment']:
        try:
            lst = [x.strip() for x in row.split('-')]
            for item in lst:
                numb_count_dict[item] += 1
        except:
            pass

    for gene in cancer_genes_df['gene_symbol']:
        numb_count_lst.append(numb_count_dict[gene])

    cancer_genes_df.insert(2, "count_seq", seqvar_count_lst)
    cancer_genes_df.insert(3, "count_numb", numb_count_lst)

    # The relevant genes are genes that show up in either the sequence or the
    # copy number mutations, drop the genes that show up in neither
    cancer_genes_relevant_df = cancer_genes_df[(cancer_genes_df['count_seq'] > 0) |
                                               (cancer_genes_df['count_numb'] > 0)]
    cancer_genes_relevant_df.to_csv('Data/cancer_genes_relevant.csv')


# This matches every cosmic sample ID to its genetic features
def consolidate_genetic_features(genetic_features_file, microsatellite_file):
    genetic_features_df = pd.read_csv(genetic_features_file)
    # Initialize empty dictionary
    genetic_features_dict = {x:[] for x in pd.unique(genetic_features_df['genetic_feature'])}
    genetic_features_dict['cosmic_sample_id'] = []
    # Construct a list of unique cosmic sample ids
    cosmic_sample_id_lst = pd.unique(genetic_features_df['cosmic_sample_id'])
    for sample in cosmic_sample_id_lst:
        # Get a slice of the data frame where the features are mutated
        test_df = genetic_features_df[(genetic_features_df['cosmic_sample_id'] == sample) &
                                      (genetic_features_df['is_mutated'] == 1)]
        # Store the mutated features in a list
        mutations_lst = test_df['genetic_feature'].tolist()
        genetic_features_dict['cosmic_sample_id'].append(sample)
        for key in genetic_features_dict:
            if key == 'cosmic_sample_id':
                pass
            elif key in mutations_lst:
                # If the feature is in the mutated features list its value is encoded as 1
                genetic_features_dict[key].append(1)
            else:
                # If it is not mutated or is missing, encode it 0
                genetic_features_dict[key].append(0)

    # Add the microsatellite instability:
    microsatellite_df = pd.read_csv(microsatellite_file, dtype={'COSMIC identifier': 'int64'})
    # Initialize the column at 0
    microsatellite_df['microsatellite'] = 0
    # If the instability is high, change it to 1
    microsatellite_df.loc[microsatellite_df['MSI']=='MSI-H', 'microsatellite'] = 1
    # Rename the key column so we can merge the two dataframes
    microsatellite_df = microsatellite_df.rename(columns={'COSMIC identifier': 'cosmic_sample_id'})
    # Select the only two columns we need
    microsatellite_df = microsatellite_df[['cosmic_sample_id', 'microsatellite']]
    genetic_features_df = pd.DataFrame(data=genetic_features_dict)
    # Merge and save
    to_save = pd.merge(left=genetic_features_df, right=microsatellite_df, on='cosmic_sample_id')
    to_save.to_csv("./Data/Clean/sample_id_features.csv", index=False)


def create_training_file(drugs_data_file, genetic_features_file, IC50_vals_file):
    """This function combines the drug, genetic features, and IC50 value files
    to create a full file that pairs features with outputs and is ready for ML purposes"""

    # We load all the files and downcast for memory efficiency
    drugs_data_df = pd.read_csv(drugs_data_file)
    drugs_data_df = drugs_data_df.apply(func=downcast_columns)
    genetic_features_df = pd.read_csv(genetic_features_file)
    genetic_features_df = genetic_features_df.apply(func=downcast_columns)
    ic50_df = pd.read_csv(IC50_vals_file)

    # We only need three columns from the ic50 file
    ic50_df = ic50_df[['Drug Id', 'Cosmic sample Id', 'IC50']]
    ic50_df = ic50_df.rename(columns={'Drug Id': 'drug_id',
                                      'Cosmic sample Id': 'cosmic_sample_id'})
    ic50_df = ic50_df.apply(func=downcast_columns)
    print("Merging ic50 and genetic features data...")
    df_train = pd.merge(ic50_df, genetic_features_df, on='cosmic_sample_id')
    print("Adding drug descriptor data...")
    df_train = pd.merge(df_train, drugs_data_df, on='drug_id')
    features = np.array(df_train.drop(columns=['IC50', 'drug_id ', 'cosmic_sample_id']))
    print("Saving to file...")
    np.save(file='./Data/Clean/features.npy', arr=features)
    outputs = np.array(df_train['IC50'])
    np.save(file='./Data/Clean/outputs.npy', arr=outputs)


def downcast_columns(column):
    if column.dtype == 'int64' or column.dtype == 'int8':
        return column.astype('int8')
    else:
        return column.astype('float16')


clean_genetic_features()
calculate_drug_data('./Data/Raw/Drugs.csv')
consolidate_genetic_features('./Data/Clean/genetic_features_clean.csv', './Data/Raw/microsattelite_data.csv')
create_training_file(drugs_data_file='./Data/Clean/drugs_fulldata.csv',
                     genetic_features_file='./Data/Clean/sample_id_features.csv',
                     IC50_vals_file='./Data/Raw/IC50_vals.csv')


