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


#clean_genetic_features()
#find_useful_mutations()
#calculate_drug_data('./Data/Raw/Drugs.csv')
genetic_features_df = pd.read_csv('./Data/clean/genetic_features_clean.csv')
genetic_features_dict = {x:[] for x in pd.unique(genetic_features_df['genetic_feature'])}
genetic_features_dict['cosmic_sample_id'] = []
cosmic_sample_id_lst = pd.unique(genetic_features_df['cosmic_sample_id'])
for sample in cosmic_sample_id_lst:
    test_df = genetic_features_df[(genetic_features_df['cosmic_sample_id'] == sample) &
                                  (genetic_features_df['is_mutated'] == 1)]
    mutations_lst = test_df['genetic_feature'].tolist()
    genetic_features_dict['cosmic_sample_id'].append(sample)
    for key in genetic_features_dict:
        if key == 'cosmic_sample_id':
            pass
        elif key in mutations_lst:
            genetic_features_dict[key].append(1)
        else:
            genetic_features_dict[key].append(0)


to_save = pd.DataFrame(data=genetic_features_dict)
to_save.to_csv("./Data/Clean/sample_id_features.csv", index=False)


