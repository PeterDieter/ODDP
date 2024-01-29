import plotly.express as px
import pandas as pd

# Read the TXT file into a Pandas DataFrame
file_path = 'clientStatistics_CFA.txt'
df = pd.read_csv(file_path, delim_whitespace=True, names=['x', 'y', 'time', 'orders', 'rejections'])
df['rejectionRate'] = df['rejections']/df['orders']
df['time'] = df['time']*0.5 + 14
print(df)
# Create a heatmap using Plotly Express
fig = px.density_heatmap(df, x='x', y='y', z='rejectionRate', histfunc='avg', nbinsx=30, nbinsy=30,range_color=[0,0.5],color_continuous_scale=['green', 'yellow', 'red'], animation_frame='time')

# Show the plot
fig.show()