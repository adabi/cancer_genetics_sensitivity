import pandas as pd
from collections import Counter
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


# This function finds mutations which show up as either sequence variation or number variation
# and stores them in cancer_genes_relevant
def find_useful_mutations():

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
    numb_count_dict = {}
    for row in genetic_features_df['genes_in_segment']:
        try:
            lst = [x.strip() for x in row.split('-')]
            for item in lst:
                try:
                    numb_count_dict[item] += 1
                except:
                    numb_count_dict[item] = 0
        except:
            pass

    for gene in cancer_genes_df['gene_symbol']:
        try:
            numb_count_lst.append(numb_count_dict[gene])
        except KeyError:
            numb_count_lst.append(0)

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
    # Construct a list of unqiue cosmic sample ids
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
                # If it is not mutated or is missing, it's encoded as 0
                genetic_features_dict[key].append(0)

    # Now we add the microsatellite instability:
    microsatellite_df = pd.read_csv('./Data/Raw/microsattelite_data.csv')
    microsatellite_df['microsatellite'] = 0
    microsatellite_df.loc[microsatellite_df['MIS'] ]
    to_save = pd.DataFrame(data=genetic_features_dict)
    to_save.to_csv("./Data/Clean/sample_id_features.csv", index=False)


#clean_genetic_features()
#calculate_drug_data('./Data/Raw/Drugs.csv')
#consolidate_genetic_features('./Data/Clean/genetic_features_clean.csv')


