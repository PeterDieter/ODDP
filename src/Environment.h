#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <assert.h>
#include <string>
#include <vector>
#include <limits.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <random>
#include <deque>

#include "Matrix.h"
#include "Data.h"
#include "Environment.h"



class Environment
{
public:
	// Constructor 
	Environment(Data* data);

	// Function to perform a simulation
	void simulate(char * argv[]);

	const double EulerConstant = std::exp(1.0);

private:
	Data* data;													// Problem parameters
	std::vector<Order*> orders;									// Vector of pointers to orders. containing information on each order
	std::vector<Order*> ordersAssignedToCourierButNotServed;	// Vector of orders that have not been served yet
	std::vector<Order*> ordersPending;							// vector of orders that have not been preocessed yet
	std::vector<Warehouse*> warehouses;							// Vector of pointers containing information on each warehouse
	std::vector<Courier*> couriers;								// Vector of pointers containing  information on each courier
	std::vector<Picker*> pickers;								// Vector of pointers  containing information on each picker
	std::vector<Route*> routes;									// Vector of pointers  containing information on each route
	Order* nextOrderBeingServed;								// Order that will be served next. Needed as we have two types of decision epoch: Order arriving and order being served (courier needs to be reassigned)
	std::vector<int> orderTimes;								// Vector of times at which clients arrive. Will be created upon initialization
	std::vector<int> clientsVector;								// Vector of clients that arrive. Same length as orderTimes vector. Will be created upon initialization
	std::vector<int> timesToComission;							// Vector of times to comission. Same length as orderTimes vector. Will be created upon initialization
	std::vector<int> timesToServe;								// Vector of times how long it takes to serve a client at his house. Same length as orderTimes vector. Will be created upon initialization
	int currentTime;											// The current time
	int nbOrdersServed;											// Number of orders served
	int timeCustomerArrives;									// Time the last customer arrived
	int timeNextCourierArrivesAtOrder;							// Time the next order is served
	int totalWaitingTime;										// Tracking the total waiting time	
	bool gridInstance;											// Bool if instance is a grid or a real city (distance measure depends on this)
	bool bundle;												// Bool if bundling orders (one courier serving multiple orders in one trip) is allowed
	bool postpone;												// Bool if we postpone the assignment decision
	int bundledOrders;											// States the total number of bundled orders
	int timeStepSize;											// As demand rates are given as a vector, this states how long one time step (1 element of the vector) is used (in seconds)

	// In this method, we reassign orders to other warehouses
	void simulation(int policy);

	// In this method we initialize the rest of the Data, such as warehouses, couriers, etc.
	void initialize();

	// Function to initialize the values of an order
	void initOrder(int currentTime, int id, Order* o);

	// Functions that assigns order to a warehouse, picker and courier, respectively
	void chooseClosestWarehouseForOrder(Order* newOrder);
	void chooseWarehouseBasedOnQuadrant(Order* newOrder);;
	void updateInformation(Order* newOrder, bool bundling);

	// Choose warehouse for an order, based on the Lower bound policy
	void chooseWarehouseForOrderReassignment(Order* newOrder, bool bundle);

	// Function that assigns a courier to the closest warehouse
	void chooseWarehouseForCourier(Courier* courier);

	// Function that adds order to a vector of orders based on the (expected) arrival time
	void addOrderToVectorArrivalTime(std::vector<Order*> & V, Order* orderToAdd);
	void addOrderToVectorDecisionTime(std::vector<Order*> & V, Order* orderToAdd);

	double getTotalWaitingTime();
	double getTotalDelays();

	// Function that returns the fastest available picker at a warehouse
	Picker* getFastestAvailablePicker(Warehouse* warehouse);

	// Function that returns the fastest available courier assigned to a warehouse
	Courier* getFastestAvailableCourier(Warehouse* warehouse);
	double insertOrderToCourierCosts(Order* newOrder, Courier* courier, bool bundle);
	std::tuple<int, Courier*>  costsToWarehouse(Order* newOrder, Warehouse* war, bool bundle);

	// Function that updates the order that will be served next
	void updateOrderBeingServedNext();

	// Function that saves a route to the list of routes
	void writeOrderStatsToClients();
	void writeClientsStatsToFile(std::string filename);

	// Functions that writes routes/orders and costs to file
	void writeCostsToFile(std::vector<float> costs, std::vector<float> averageDelayRateVector, float lambdaTemporal, float lambdaSpatial, bool is_training);
	void writeStatsToFile(std::vector<float> costs, std::vector<float> averageDelayRateVector, std::vector<float> averageWaitingTime, std::vector<float> maxWaitingTime);
	void writeMatrixToFile(std::vector<std::vector<double>> matrix, std::string filename);
	void writeCourierRoutesToFile(std::string fileNameRoutes, std::string fileNameOrders);
	
	// Function to draw an inter arrival time based on rate specified in data
	int drawFromExponentialDistribution(double lambda);
	
	// new member functions used for postponement
	int calcTimeAndEvent(int,int&);
	int calcNewdecisionTime(Order* newOrder);
	
};


#endif
