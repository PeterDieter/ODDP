import pandas as pd
import plotly.graph_objects as go
import json
import numpy as np



def plotData(clients, warehouses, limit):
    df = clients.to_numpy()
    matrix = df[:,3:].astype(int)
    minValues = matrix.min(axis=1)
    notInLimit = np.where(minValues > limit)[0]
    clients = clients.drop(notInLimit)
    
    fig = go.Figure(go.Scattermapbox(lat=clients.Latitude, lon=clients.Longitude, marker=go.scattermapbox.Marker(
                size=5,
                color='grey',
                opacity=1
            )))

    fig.add_trace(go.Scattermapbox(
            lat=[val.get('latitude') for val in warehouses.values()],
            lon=[val.get('longitude') for val in warehouses.values()],
            mode='markers',
            marker=dict(
                symbol = "circle",
                size=23,
                color='black',
                opacity=1
            ),
            hoverinfo='text'
        ))
    fig.update_layout(mapbox_style="carto-positron", mapbox_center_lon=-87.7493, mapbox_center_lat=41.958, mapbox_zoom=9)
    fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
    fig.show()


if __name__ == "__main__":
    with open('data/getirStores.json') as fp:
        getirStores = json.load(fp)

    Stopls = pd.read_csv("data/allDurations15.csv", header=0) 

    plotData(Stopls, getirStores, 500)