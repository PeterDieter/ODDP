import pandas as pd
import numpy as np
import json
import random
import geopandas as gpd
import utm
import math
import itertools
import collections


def create_instance(fileName: str, limit: int=900, totalCouriers: int=80, pickersPerWarehouse: int=3, meanComissionTime: int=120, meanServiceTimeAtClient: int=60, gridStepSize: int=300):
    """Create a .txt file of a problem instance
    Args:
        fileName (str): File in which we save the instance parameters
        limit (int): Clients whose minimum distance (in seconds) to a warehouse is above this limit are excluded 
        couriersPerWarehouse (int): number of couriers that a warehouse has at the start
        pickersPerWarehouse (int): number of pickers that a warehouse has
        meanComissionTime (int): Time a picker needs on average to comission an order (expoentially distributed) (in seconds)
        meanServiceTimeAtClient (int): Mean time a courier needs at the clients door to deliver the order (expoentially distributed) (in seconds)
    Returns:
        None 
    """
    
    df = pd.read_csv("data/allDurationsGenerated.csv", header=0) 
    df = df.to_numpy()
    random.seed(422)
    rndIdxs = random.sample(range(len(df)), round(len(df)*0.75))
    
    # clients = df[rndIdxs,:3]
    # matrix = df[rndIdxs,3:].astype(int)
    
    # Comment the following block out in case of creating test instance
    clients = df[:,:3]
    matrix = df[:,3:].astype(int)
    # clients = np.delete(clients, rndIdxs, axis=0)
    # matrix = np.delete(matrix, rndIdxs, axis=0)

    # Remove clients where distance to warehouses is over limit
    minValues = matrix.min(axis=1)
    notInLimit = np.where(minValues > limit)[0]
    clients = np.delete(clients, notInLimit, axis=0)
    matrix = np.delete(matrix, notInLimit, axis=0)

    with open('data/getirStores.json') as fp:
        getirStores = json.load(fp)

    warehouses = np.array([[val.get('longitude'),val.get('latitude')] for val in getirStores.values()])
    warehouses = np.c_[warehouses, np.ones(len(warehouses)), np.ones(len(warehouses))*pickersPerWarehouse]
    argminVec = matrix.argmin(axis=1)
    closestWarehouseCounter = collections.Counter(argminVec)
    distribution = np.ones(len(warehouses))
    for key, value in closestWarehouseCounter.items():
        distribution[key] = matrix[np.where(argminVec== key)[0],key].sum()
    
    summ = 0
    shareVector = []
    for key, value in closestWarehouseCounter.items():
        warehouses[key,2] = int(math.floor(distribution[key]/sum(distribution)*totalCouriers))
        shareVector.append(distribution[key]/sum(distribution)*totalCouriers% 1)
        summ += int(round(distribution[key]/sum(distribution)*totalCouriers))
    
    print(summ, totalCouriers, shareVector)

    if(totalCouriers-summ > 0):
        for i in range(len(shareVector)):
            #if shareVector[i] > 0.5:            
            #    shareVector[i] = shareVector[i] - 1
            sorted_list = sorted(shareVector, reverse=True)
            res = [i for i,x in enumerate(shareVector) if x in itertools.islice(sorted_list, totalCouriers-summ)]
        for j in res:
            summ += 1
            warehouses[j,2] += 1
    elif(totalCouriers-summ < 0):
        for i in range(len(shareVector)):
            if shareVector[i] > 0.5:            
                shareVector[i] = shareVector[i] - 1
            sorted_list = sorted(shareVector, reverse=False)
            res = [i for i,x in enumerate(shareVector) if x in itertools.islice(sorted_list, summ-totalCouriers)]
        for j in res:
            summ -= 1
            warehouses[j,2] -= 1
    
    print(summ)
    print(warehouses)
    # Grid Stuff
    points = gpd.read_file('data/Polygon_1200s/Polygon_1200s.shp')
    xmin, ymin, xmax, ymax = points.total_bounds

    # Create corners of rectangle to be transformed to a grid
    sw = utm.from_latlon(ymin, xmin)
    ne = utm.from_latlon(ymax, xmax)

    # Iterate over 2D area
    gridpointsAll = []
    x = sw[1]
    while x < ne[1]:
        y = sw[0]
        gridpoints = []
        while y < ne[0]:
            p = utm.to_latlon(y, x, sw[2], sw[3])
            gridpoints.append(p)
            y += gridStepSize
        gridpointsAll.append(gridpoints)
        x += gridStepSize

    counter = 0
    gridDict = {}
    list1, list2 = [], []
    for rowIndex, row in enumerate(gridpointsAll):
        if rowIndex == len(gridpointsAll)-1:
            continue
        for colIndex, col in enumerate(row):
            if colIndex == len(row)-1:
                continue
            counter += 1
            list1.append(gridpointsAll[rowIndex+1][colIndex+1][1]-gridpointsAll[rowIndex][colIndex][1])
            list2.append(gridpointsAll[rowIndex+1][colIndex+1][0]-gridpointsAll[rowIndex][colIndex][0])

    stepLat = round(sum(list1) / len(list1), 8)
    stepLon = round(sum(list2) / len(list2), 8)


    with open("instances/"+fileName+".txt", 'w') as f:
        f.write("\n".join([
            "{} : {}".format(k, v)
            for k, v in [
                ("NAME", fileName),
                ("NUMBER_CLIENTS", len(clients)),
                ("NUMBER_WAREHOUSES", len(warehouses)),
                ("MEAN_COMMISSION_TIME", meanComissionTime),
                ("MEAN_SERVICE_AT_CLIENT_TIME", meanServiceTimeAtClient)]
        ]))
        f.write("\n")

        f.write("\n".join([
            "{} : {} {}".format(k, v, u)
            for k, v, u in [
                ("GRID_SW_COORDINATES", xmin, ymin),
                ("GRID_NE_COORDINATES", xmax, ymax),
                ("STEPSIZE_LAT", stepLat, ""),
                ("STEPSIZE_LON", stepLon, "")]
        ]))
        f.write("\n")
        
        f.write("WAREHOUSE_SECTION\n")
        f.write("\n".join([
            "{}\t{}\t{}\t{}\t{}".format(i, y, z, int(q), int(p))
            for i, (y, z, q, p) in enumerate(warehouses)
        ]))
        f.write("\n")


        f.write("CLIENT_SECTION\n")
        f.write("\n".join([
            "{}\t{}\t{}".format(i, y, z)
            for i, (x, y, z) in enumerate(clients)
        ]))
        f.write("\n")

        f.write("EDGE_WEIGHT_SECTION\n")
        for row in matrix:
            f.write("\t".join(map(str, row)))
            f.write("\n")
        
        
        f.write("EOF\n")


def create_grid_instance(fileName: str, totalCouriers: int=20, pickersPerWarehouse: int=3, meanComissionTime: int=120, meanServiceTimeAtClient: int=60):
    """Create a .txt file of a problem instance
    Args:
        fileName (str): File in which we save the instance parameters
        limit (int): Clients whose minimum distance (in seconds) to a warehouse is above this limit are excluded 
        couriersPerWarehouse (int): number of couriers that a warehouse has at the start
        pickersPerWarehouse (int): number of pickers that a warehouse has
        meanComissionTime (int): Time a picker needs on average to comission an order (expoentially distributed) (in seconds)
        meanServiceTimeAtClient (int): Mean time a courier needs at the clients door to deliver the order (expoentially distributed) (in seconds)
    Returns:
        None 
    """
    warehouses = np.array([[24,24],
                            [24,74],
                            [74,24],
                            [74,74]])
    clients = np.array([[i, j] for i in range(100) for j in range(100) if not ([i,j] == warehouses).all(1).any()])

    matrix = np.zeros((len(clients),len(warehouses)))
    clientCounter, warehouseCounter = 0, 0
    for client in clients:
        for warehouse in warehouses:
            dist = abs(client[0] - warehouse[0]) + abs(client[1] - warehouse[1])
            matrix[clientCounter, warehouseCounter] = int(dist)*10
            warehouseCounter += 1
        warehouseCounter = 0
        clientCounter += 1
    matrix = matrix.astype(int)

    warehouses = np.c_[warehouses, np.ones(len(warehouses)), np.ones(len(warehouses))*pickersPerWarehouse]
    closestWarehouseCounter = collections.Counter(matrix.argmin(axis=1))
    for key, value in closestWarehouseCounter.items():
        warehouses[key,2] = int(round(value/len(matrix)*totalCouriers))
    

    with open("instances/"+fileName+".txt", 'w') as f:
        f.write("\n".join([
            "{} : {}".format(k, v)
            for k, v in [
                ("NAME", fileName),
                ("NUMBER_CLIENTS", len(clients)),
                ("NUMBER_WAREHOUSES", len(warehouses)),
                ("MEAN_COMMISSION_TIME", meanComissionTime),
                ("MEAN_SERVICE_AT_CLIENT_TIME", meanServiceTimeAtClient)]
        ]))
        f.write("\n")
        
        f.write("WAREHOUSE_SECTION\n")
        f.write("\n".join([
            "{}\t{}\t{}\t{}\t{}".format(i, y, z, int(q), int(p))
            for i, (y, z, q, p) in enumerate(warehouses)
        ]))
        f.write("\n")


        f.write("CLIENT_SECTION\n")
        f.write("\n".join([
            "{}\t{}\t{}".format(i, y, z)
            for i, (y, z) in enumerate(clients)
        ]))
        f.write("\n")

        f.write("EDGE_WEIGHT_SECTION\n")
        for row in matrix:
            f.write("\t".join(map(str, row)))
            f.write("\n")
        
        
        f.write("EOF\n")



if __name__ == "__main__":
    create_instance(fileName = "zip", limit=900, totalCouriers=80, pickersPerWarehouse=3, meanComissionTime=180, meanServiceTimeAtClient=120, gridStepSize=800)
    create_grid_instance(fileName = "grid", totalCouriers=40, pickersPerWarehouse=5, meanComissionTime=180, meanServiceTimeAtClient=120)