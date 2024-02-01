import geopandas as gpd
import plotly.graph_objects as go
import pandas as pd
import json
import math
import sys
import numpy as np
from shapely.geometry import Point


def eucl_dist(x1, x2, y1, y2):
    return math.sqrt((x1-x2)**2 + (y1-y2)**2)

def getClosestWarehouse(warehouses, point):
    bestWarehouseIdx, counter = -1, 0
    bestDist = float('inf')
    for key in warehouses:
        dist = eucl_dist(point.y, warehouses[key]['latitude'],point.x, warehouses[key]['longitude'])
        if dist < bestDist:
            bestDist = dist
            bestWarehouseIdx = counter
        counter += 1
    #print(bestWarehouseIdx,dist)
    return bestWarehouseIdx

def plotData(clients, warehouses, quads, polygon, colorEmpty, plotCustomers):
    df = pd.DataFrame(clients, columns =['ID', 'Longitude', 'Latitude'])
    fig = go.Figure(go.Scattermapbox())
    
    colors = { 'Red' : (255,0,0), 'White' : (255,255,255), 'Black' : (0, 0, 0), 'Lime': (0,255,0), 'Blue': (0,0,255), 'Yellow': (255,255,0), 'Cyan' :(0,255,255), 'Fuchsia': (255,0,255), 'Purple': (128,0,128), 'Teal' : (0,128,128), 'Navy' : (0,0,128), 'Maroon' : (128,0,0), 'Olive' : (128,128,0), 'Green' : (0,128,0) }
    groupedQ1 = quads.groupby(['WarehouseID'])
    if colorEmpty:
        new_quads = pd.DataFrame()
        for name, group in groupedQ1:
            for index, row in group.iterrows():
                if (row["WarehouseID"] != -1):
                    new_quads = pd.concat([new_quads,group])
                    break
                else:
                    pnt = Point((row["LonSW"]+row["LonNE"])/2, (row["LatNE"] + row["LatSW"])/2)
                    if (pnt.within(polygon["geometry"][0])):
                        group.loc[index, "WarehouseID"] = getClosestWarehouse(warehouses, pnt)
            new_quads = pd.concat([new_quads,group])
        new_quads = new_quads.drop_duplicates(subset=['ID'])
    else:
        new_quads = quads
    groupedQ2 = new_quads.groupby(['WarehouseID'])
    counter = 0
    for name, group2 in groupedQ2:
        allLat, allLon = [], []
        
        for index, row in group2.iterrows():
            allLon.extend((row["LonSW"],row["LonNE"],row["LonNE"],row["LonSW"], None)),
            allLat.extend((row["LatSW"],row["LatSW"],row["LatNE"],row["LatNE"], None)),
        
        if counter != 999:
            fig.add_trace(go.Scattermapbox(
                fill = "toself",
                lon = allLon,
                lat = allLat,
                marker = dict(size = 0, color = list(colors)[counter]),
            ))
        counter += 1

    if plotCustomers:
        fig.add_trace(go.Scattermapbox(
            lat=df.Latitude, lon=df.Longitude,
            marker=go.scattermapbox.Marker(
                color='blue',
                opacity=0.1,
            ),
        ))

    fig.add_trace(go.Scattermapbox(
            lat=[val.get('latitude') for val in warehouses.values()],
            lon=[val.get('longitude') for val in warehouses.values()],
            mode='markers',
            marker=dict(
                size=20,
                color='black',
                opacity=1,
                symbol = "circle"
            ),
            hoverinfo='text'
        ))

    fig.update_layout(mapbox_style="carto-positron", mapbox_center_lon=-87.7493, mapbox_center_lat=41.958, mapbox_zoom=9)
    fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
    fig.show()

if __name__ == "__main__":
    
    fill = False
    cust = False
    if ( (len(sys.argv) > 3) ):
      print("Usage: plotGrid.py [-f] [-c]")
      print("\t- [-f]\t fill gridf\t (default False)")
      print("\t- [-c]\t draw customers\t (default False)")
      sys.exit()
    # process arguments
    elif ( (len(sys.argv) > 1) ):
      for i in range(1,len(sys.argv)):
        if ( sys.argv[i] == "-f" or sys.argv[i] == "--fill" ):
          fill = True
        if ( sys.argv[i] == "-c" or sys.argv[i] == "--customer" ): 
          cust = True

    
    quads = pd.read_csv("data/heatmapData/wh_quad_assignment.txt", delimiter=" ")
    quads.columns =['ID', 'LatSW', 'LonSW', 'LatNE', 'LonNE', 'WarehouseID', 'w1', 'w2']

    with open('data/getirStores.json') as fp:
        getirStores = json.load(fp)

    
    df = pd.read_csv("data/allDurations15.csv", header=0) 
    df = df.to_numpy()
    clients = df[:,:3]
    matrix = df[:,3:].astype(int)

    # Remove clients where distance to warehouses is over limit
    minValues = matrix.min(axis=1)
    notInLimit = np.where(minValues > 900)[0]
    clients = np.delete(clients, notInLimit, axis=0)
    matrix = np.delete(matrix, notInLimit, axis=0)

    polygon = gpd.read_file("data/Polygon_900s/Polygon_900s.shp")

    plotData(clients, getirStores, quads, polygon, fill, cust)
