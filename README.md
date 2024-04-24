# The On-Demand Delivery Problem: Online Assignment of Orders to Warehouses and Couriers

In the problem at hand, we operate multiple warehouse (depots) in a service region from which couriers start their trip to serve an order. When an order arrives, we need to assign it to a warehouse and a courier. When a courier serves the last order on his tour, he needs to be assigned to a warehouse to return to. The goal is to minimize the tardiness of customers.

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
./onlineAssignment --instance=pathToInstance --maxWaiting=int --AMethod=AMethodName --RMethod=RMethodName --b --a=float
```

where **instanceName** gives the path to the .txt file containing the instance information. The maxWaiting parameter give the maximal waiting time in seconds (int). To enable bundling of customers to one courier trip, add --b. The **AMethod** is a string that determines the method which will be applied for the assignment problem. The following assigning strategies are available:

1. n: In this policy, we assign each order to the nearest warehouse and the earliest courier available.
2. r: We check if order can be assigned to any warehouse. We choose the courier with lowest waiting for the order.
2. w: We check if order can be assigned to any warehouse. We weigh waiting time and travel time, to anticipate future demand.
3. s: Use Gurobi to partition the service region based on quadrants. 

If the weighted policy is chosen, a weighting value **a** needs to be chosen, that is a float between 0 and 1.
Concerning rebalancing, **RMethod** can take the following values:

1. s: Static. No rebalancing occurs.
2. n: Nearest: Always assign the courier to the nearest warehouse.
3. l: Level: Assign the courier to the nearest warehouse but check that each warehouse remains at a certain level.

For example:

```
./onlineAssignment --instance=instances/zip.txt --maxWaiting=1200 --AMethod=w --RMethod=s --b --a=0.6
```