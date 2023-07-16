# split into multiple columns, 1 column for every mutation
def split_m1_to_m15(df):
    df[['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']] = \
        df['mutant'].str.split(':', 14, expand=True)
    return df

def reduce_df_to_DMSscore_mutations(df):
    df = df.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
    return df

