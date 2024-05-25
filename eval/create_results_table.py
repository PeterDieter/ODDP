import pandas as pd

grid = False
max_wait = 1200

df = pd.read_csv("eval/data_dafd.csv")
df = df[df['max_wait'] == max_wait]
df = df[df['instance'] == grid]


idx = df.groupby(['a_meth', 'r_meth', 'bundling'])['averageDelay'].idxmin()

# Create the new DataFrame using these indices
result_df = df.loc[idx].reset_index(drop=True)
result_df = result_df[result_df['r_meth'] != 'n']
result_df = result_df.sort_values(by=['bundling', 'a_meth'], ascending = False)
result_df = result_df.drop(columns=['instance', 'max_wait', 'alpha', 'beta'])

print(result_df.to_latex(index=False, float_format="%.2f"))