# The On-Demand Delivery Problem: Online Assignment of Orders to Warehouses and Couriers

In the problem at hand, we operate multiple warehouse (depots) in a service region from which couriers start their trip to serve an order. When an order arrives, we need to assign it to a warehouse and a courier. The goal is to minimize the tardiness of customers.

Visualization of the problem (made in [visualizeSimulation.py](python/visualizeSimulation.py)):

<p align="center">
<img src="animation.gif" width="300" height="400" align="center">
</p>


## Data
For the warehouses, we use Getir stores in Chicago. For the order data, we use publicly available anonymized customer data of Amazon customers in Chicago or self generated customers based on the population in a zip code area. We then draw random orders from these customers. Interarrival times, service times and comission times are all assumed to be exponentially distributed. An instance can be created in [createInstance.py](python/createInstance.py):


## C++ compiling 
Data is prepared in Python and an instance is then passed to C++. The raw and processed data is contained in folder [data](data). All code related to preprocessing data (including Isochrone API and DistanceMatrix API) and creating code to create instances is contained in [python](python).

```
make clean
make
```

## Running the program

You can then execute the code with:

```
./onlineAssignment instanceName maxWaiting methodName
```

where **instanceName** gives the path to the .txt file containing the instance information. The third parameter is the **max waiting time** in seconds (int). The **methodName** is a string that determines the method which will be applied for the assignment problem. For example:

```
./onlineAssignment instances/grid.txt 1800 nearestWarehouse
```

Currently, the following assigning strategies are available:
1. nearestWarehouse: In this policy, we assign each order to the nearest warehouse. If order cannot be served on time, we reject.
2. reassignment: We check if order can be assigned to any warehouse. We choose the warehouse with lowest waiting for the order.
3. staticPartitioning: Use Gurobi to partition the service region based on quadrants. 