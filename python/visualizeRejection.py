import plotly.express as px
import pandas as pd

# Read the TXT file into a Pandas DataFrame
file_path = 'clientStatistics_nearest.txt'
df = pd.read_csv(file_path, delim_whitespace=True, names=['x', 'y', 'time', 'orders', 'rejections'])
df['rejectionRate'] = df['rejections']/df['orders']

file_path = 'clientStatistics_CFA.txt'
df2 = pd.read_csv(file_path, delim_whitespace=True, names=['x', 'y', 'time', 'orders', 'rejections'])
df2['rejectionRate'] = df2['rejections']/df2['orders']

df2['rejectionDifference'] = df2['rejectionRate'] - df['rejectionRate']

df2['time'] = df2['time']*0.5 + 14
print(df2)
# Create a heatmap using Plotly Express
fig = px.density_heatmap(df2, x='x', y='y', z='rejectionDifference', histfunc='avg', nbinsx=30, nbinsy=30,range_color=[-0.15,0.15],color_continuous_scale=['green', 'yellow', 'red'], animation_frame='time')

# Show the plot
fig.show()