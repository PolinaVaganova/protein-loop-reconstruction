import pandas as pd

path_to_annotation_df = '../0_collect_data/annotations/all.csv'

annotation_df = pd.read_csv(path_to_annotation_df, sep=',')

mask_with_space_group = (annotation_df['gap_num'].astype('float64') == 1) & (
        annotation_df['chain_num'].astype('float64') == 1) & (4 <= annotation_df['gap_len'].astype('float64')) & (
                                annotation_df['gap_len'].astype('float64') <= 8) & (
                                annotation_df['is_gap_in_ss'] == 'False') & (
                                annotation_df['is_gap_on_terminal_end'] == 'False') & (
                                annotation_df['hetatm'].astype('float64') == 0) & (
                                annotation_df['unk'].astype('float64') == 0) & (
                                8 <= annotation_df['a'].astype('float64')) & (
                                annotation_df['a'].astype('float64') <= 100) & (
                                8 <= annotation_df['b'].astype('float64')) & (
                                annotation_df['b'].astype('float64') <= 100) & (
                                8 <= annotation_df['c'].astype('float64')) & (
                                annotation_df['c'].astype('float64') <= 100)

filtered_with_space_groups = annotation_df[mask_with_space_group]
filtered_with_space_groups.to_csv('monomers.csv', index=False)

mask_without_space_group = (annotation_df['gap_num'].astype('float64') == 1) & (
        annotation_df['chain_num'].astype('float64') == 1) & (4 <= annotation_df['gap_len'].astype('float64')) & (
                                   annotation_df['gap_len'].astype('float64') <= 8) & (
                                   annotation_df['is_gap_in_ss'] == 'False') & (
                                   annotation_df['is_gap_on_terminal_end'] == 'False') & (
                                   annotation_df['hetatm'].astype('float64') == 0) & (
                                   annotation_df['unk'].astype('float64') == 0) & (
                                   8 <= annotation_df['a'].astype('float64')) & (
                                   annotation_df['a'].astype('float64') <= 100) & (
                                   8 <= annotation_df['b'].astype('float64')) & (
                                   annotation_df['b'].astype('float64') <= 100) & (
                                   8 <= annotation_df['c'].astype('float64')) & (
                                   annotation_df['c'].astype('float64') <= 100) & (
                                   annotation_df['Z'].astype('float64') > 7)

filtered_without_space_groups = annotation_df[mask_without_space_group]
filtered_without_space_groups.to_csv('monomers_to_ucell.csv', index=False)

mask_without_space_group = (annotation_df['gap_num'].astype('float64') == 1) & (
        annotation_df['chain_num'].astype('float64') == 1) & (4 <= annotation_df['gap_len'].astype('float64')) & (
                                   annotation_df['gap_len'].astype('float64') <= 8) & (
                                   annotation_df['is_gap_in_ss'] == 'False') & (
                                   annotation_df['is_gap_on_terminal_end'] == 'False') & (
                                   annotation_df['hetatm'].astype('float64') == 0) & (
                                   annotation_df['unk'].astype('float64') == 0) & (
                                   8 <= annotation_df['a'].astype('float64')) & (
                                   annotation_df['a'].astype('float64') <= 100) & (
                                   8 <= annotation_df['b'].astype('float64')) & (
                                   annotation_df['b'].astype('float64') <= 100) & (
                                   8 <= annotation_df['c'].astype('float64')) & (
                                   annotation_df['c'].astype('float64') <= 100) & (
                                   annotation_df['Z'].astype('float64') <= 7)

filtered_without_space_groups = annotation_df[mask_without_space_group]
filtered_without_space_groups.to_csv('monomers_to_supercell.csv', index=False)
