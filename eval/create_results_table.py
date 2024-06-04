import pandas as pd
import matplotlib.pyplot as plt
import tikzplotlib

grid = True
max_wait = 1500

df1 = pd.read_csv("eval/data_dafd.csv")
df = df1[df1['max_wait'] == max_wait]
df = df[df['instance'] == grid]

idx = df.groupby(['a_meth', 'r_meth', 'bundling'])['averageDelay'].idxmin()

# Create the new DataFrame using these indices
result_df = df.loc[idx].reset_index(drop=True)
result_df = result_df[result_df['r_meth'] == 's']
result_df = result_df.sort_values(by=['bundling', 'a_meth'], ascending = False)
result_df = result_df.drop(columns=['instance', 'max_wait', 'alpha', 'beta', 'percReassignedCouriers'])

print(result_df.to_latex(index=False, float_format="%.2f"))
df1 = df1[df1['a_meth'] == "w"]
df1 = df1[df1['r_meth'] == "s"]
df1 = df1.sort_values(by=['instance', 'max_wait', 'alpha'], ascending = False)
df1 = df1[df1['alpha'] != 1]
df1 = df1.drop(columns=['beta', 'percReassignedCouriers'])
print(df1)


df1 = df1[df1['instance'] == False]
df1 = df1[df1['bundling'] == False]

# Unique values of max_wait
max_wait_values = df1['max_wait'].unique()
bundling_values = df1['bundling'].unique()
instance_values = df1['instance'].unique()
print(df1['alpha'].unique())

# Grayscale colors
grayscale_colors = ['0.2', '0.5', '0.8']

# Plotting
fig = plt.figure(figsize=(10, 6))
i = 0
for max_wait in max_wait_values:
    for instance in instance_values:
        for bundle in bundling_values:
            subset = df1[df1['max_wait'] == max_wait]
            subset = subset[subset['instance'] == instance]
            subset = subset[subset['bundling'] == bundle]
            plt.plot(subset['alpha'], subset['averageDelay'], label=f'max_wait = {max_wait}', marker='o', color=grayscale_colors[i])
            i += 1

def tikzplotlib_fix_ncols(obj):
    """
    workaround for matplotlib 3.6 renamed legend's _ncol to _ncols, which breaks tikzplotlib
    """
    if hasattr(obj, "_ncols"):
        obj._ncol = obj._ncols
    for child in obj.get_children():
        tikzplotlib_fix_ncols(child)


# Adding labels and title
plt.xlabel('Alpha')
plt.ylabel('Average Delay')
plt.title('Average Delay vs Alpha for Different Max Wait Times')
plt.xlim(0.35,1)
plt.ylim(-1,70)
plt.legend()
plt.grid(True)
tikzplotlib_fix_ncols(fig)
tikzplotlib.save("test.tex")
plt.show()