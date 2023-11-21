import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def rarefy_dataframe(df, depth):
    # Initialize an empty DataFrame to hold the rarefied counts
    rarefied_df = pd.DataFrame(index=df.index, columns=df.columns)

    # For each row in the DataFrame
    for index, row in df.iterrows():
        # If the sum of the row is less than the depth, just copy the row
        if row.sum() <= depth:
            rarefied_df.loc[index] = row
        else:
            # Otherwise, randomly subsample the row to the specified depth
            rarefied_row = subtract_from_list(row, depth)
            # Update the rarefied DataFrame with the sampled counts
            rarefied_df.loc[index] = rarefied_row
    return rarefied_df

def subtract_from_list(lst, target):
    # Calculate the difference between the sum of the list and the target sum
    diff = sum(lst) - target
    # Ensure the difference is not greater than the sum of the list
    if diff > 0:
        # Randomly distribute this difference across the list
        while diff > 0:
            for i in range(len(lst)):
                if diff > 0 and lst[i] > 0:
                    lst[i] -= 1
                    diff -= 1
    return lst

def main():

    test_matrix='/home/jbyoung/scratch/test_data/20200930-113836_count.txt'

    df = pd.read_csv(test_matrix, index_col=0, sep='\t')
    
    totals=df.sum(axis=0)
    max_depth = totals.max()
    subsampling_steps = np.linspace(1, max_depth, num=10).astype(int)
    ## transpose the df
    df_t = df.transpose()
    #df_t.set_index(df_t.columns[0], inplace=True)
    
    df_list = []
    for i in subsampling_steps[5:]:
        rarefied_df = rarefy_dataframe(df_t, i)
        rarefied_df['depth'] = i
        df_list.append(rarefied_df)
    df_list
    melted_dfs = [x.reset_index().melt(id_vars=['index', 'depth']) for x in df_list]

    plot_df = pd.concat(melted_dfs)

    summary = plot_df[plot_df['value'] >= 10].groupby(['index', 'depth']).size()

    