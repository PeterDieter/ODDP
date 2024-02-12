import geojson
from shapely.geometry import shape, Point
import pandas as pd
import random
import csv
import shapely.wkt


my_data = pd.read_csv('data/zipbypopulation.csv', header=None, index_col=0).squeeze('columns').to_dict()
print(my_data)
total_population = sum(my_data.values())

with open("data/ZipCodesChicago.json") as f:
    gj = geojson.load(f)
gj = gj["data"]

def polygon_random_points (poly, num_points, zipCode):
    min_x, min_y, max_x, max_y = poly.bounds
    points = []
    while len(points) < num_points:
            random_point = Point([random.uniform(min_x, max_x), random.uniform(min_y, max_y)])
            if (random_point.within(poly)):
                points.append(random_point)
    return [[point.xy[0][0], point.xy[1][0], zipCode] for point in points ]

customerList = []
total_customers = 50000
for zipCode in range(len(gj)):
    pol = shapely.wkt.loads(gj[zipCode][8])
    # if (int(features[zipCode]["properties"]["ZIP"]) in my_data):
    numCustomersInZIP = int(total_customers* my_data[int(gj[zipCode][10])]/total_population)
    points = polygon_random_points(pol, numCustomersInZIP, int(int(gj[zipCode][10])))
    customerList += points

fields = ['Longitude', 'Latitude', 'Zip']
with open('data/generated_customers.csv', 'w') as f:
     
    # using csv.writer method from CSV package
    write = csv.writer(f)
    write.writerow(fields)
    write.writerows(customerList)

