import pandas as pd
import tikzplotlib
import matplotlib.pyplot as plt 
import plotly.graph_objects as go
import json

instance = "zip"

# Read the TXT file into a Pandas DataFrame
file_path = 'eval/clientStatistics_nearest_' + instance + '.txt'
df = pd.read_csv(file_path, delim_whitespace=True, names=['x', 'y', 'time', 'orders', 'delayTime', 'delayOccurence', 'bundled'])
df['delayRate'] = df['delayTime']/df['orders']

file_path = 'eval/clientStatistics_CFA_' + instance + '.txt'
df2 = pd.read_csv(file_path, delim_whitespace=True, names=['x', 'y', 'time', 'orders', 'delayTime', 'delayOccurence', 'bundled'])
df2['delayRate'] = df2['delayTime']/df2['orders']

df2['difference'] = (df['delayRate']-df2['delayRate'])/(df['delayRate'])# / df['delayRate']
df3 = df2.groupby(['x', 'y']).mean()
df2 = df2[df2['time'] == 11] 
df2 = df2.fillna(0)
df3 = df3.reset_index()
print(df3)

if (instance == "grid"):
    print(type(df3))
    heatmap_data = df3.pivot('y', 'x', 'delayRate')

    # Plot the heatmap
    plt.imshow(heatmap_data, cmap='viridis', origin='lower', vmin=0, vmax=20)
    plt.colorbar()

    # Add labels
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D Heatmap')

    tikzplotlib.save("test.tex")
    # Show plot
    plt.show()
else:
    with open('data/getirStores.json') as fp:
        warehouses = json.load(fp)

    # Create the density heatmap
    fig = go.Figure()

    fig.add_trace(go.Densitymapbox(
        lat=df3['x'],
        lon=df3['y'],
        
        z=df3['delayRate'],
        radius=5,
        showscale=True,
        colorscale='viridis',
        zmin = 0,
        zmax = 60
    ))

    fig.add_trace(go.Scattermapbox(
                lat=[val.get('latitude') for val in warehouses.values()],
                lon=[val.get('longitude') for val in warehouses.values()],
                mode='markers+text',
                marker=dict(
                    symbol = "circle",
                    size=23,
                    color='black',
                    opacity=0.8,
                ),
                textposition='top right',
                textfont=dict(size=16, color='black'),
                text=['Warehouse: ' + str(i) for i in range(len(warehouses))]
            ))

    # Set the layout for the map
    fig.update_layout(
        mapbox_style="carto-positron",
        mapbox_center_lat=41.95,
        mapbox_center_lon=-87.75,
        mapbox_zoom=10.8,
        margin={"r":0,"t":0,"l":0,"b":0}
    )
    # Show the figure
    fig.show()

