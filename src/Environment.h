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

#include <torch/torch.h>
#include <torch/script.h>
#include "Matrix.h"
#include "Data.h"
#include "Environment.h"

struct policyNetwork;

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
	std::vector<Warehouse*> warehouses;							// Vector of pointers containing information on each warehouse
	std::vector<Courier*> couriers;								// Vector of pointers containing  information on each courier
	std::vector<Picker*> pickers;								// Vector of pointers  containing information on each picker
	std::vector<Route*> routes;									// Vector of pointers  containing information on each route
	Order* nextOrderBeingServed;								// Order that will be served next. Needed as we have two types of decision epoch: Order arriving and order being served (courier needs to be reassigned)
	std::vector<int> orderTimes;								// Vector of times at which clients arrive. Will be created upon initialization
	std::vector<int> clientsVector;								// Vector of clients that arrive. Same length as orderTimes vector. Will be created upon initialization
	std::vector<int> timesToComission;							// Vector of times to comission. Same length as orderTimes vector. Will be created upon initialization
	std::vector<int> timesToServe;								// Vector of times how long it takes to serve a client at his house. Same length as orderTimes vector. Will be created upon initialization
	int currentTime;
	int nbOrdersServed;
	int rejectCount;
	int timeCustomerArrives;
	int timeNextCourierArrivesAtOrder;
	int totalWaitingTime;
	int highestWaitingTimeOfAnOrder;
	int latestArrivalTime;
	torch::Tensor assingmentProblemStates;
	torch::Tensor assingmentProblemActions;

	// In this method we apply the nearest warehouse policy.
	void tuneK(int timelimit);

	// In this method, we reassign orders to other warehouses
	void reassignmentPolicyLB(int timelimit);

	// In this method we initialize the rest of the Data, such as warehouses, couriers, etc.
	void initialize(int timeLimit, std::vector<int> vectorOfKs);

	// Function to initialize the values of an order
	void initOrder(int currentTime, Order* o);

	// Functions that assigns order to a warehouse, picker and courier, respectively
	void chooseClosestWarehouseForOrder(Order* newOrder);
	void choosePickerForOrder(Order* newOrder);
	void chooseCourierForOrder(Order* newOrder);

	// Choose warehouse for an order, based on the Lower bound policy
	void chooseWarehouseForOrderLB(Order* newOrder);
	double getOpportunityCostsLB(Warehouse* w, Order* o);
	std::vector<int> getDonePickingTimes(Warehouse* w, int endTime);

	// Function that assigns a courier to the closest warehouse
	void chooseClosestWarehouseForCourier(Courier* courier);

	// Function that deletes order from ordersNotServed vector
	void RemoveOrderFromVector(std::vector<Order*> & V, Order* orderToDelete);

	// Function that adds order to a vector of orders based on the (expected) arrival time
	void AddOrderToVector(std::vector<Order*> & V, Order* orderToAdd);

	// Function that returns the euclidean distance between two locations
	double euclideanDistance(double latFrom, double latTo, double lonFrom, double lonTo);

	// Function that deletes order from ordersNotServed vector
	void RemoveCourierFromVector(std::vector<Courier*> & V, Courier* courierToDelete);

	// Function that returns the fastest available picker at a warehouse
	Picker* getFastestAvailablePicker(Warehouse* warehouse);

	// Function that returns the fastest available courier assigned to a warehouse
	Courier* getFastestAvailableCourier(Warehouse* warehouse);

	int getNumberOfAvailablePickers(Warehouse* warehouse);
	int getNumberOfOrdersNotPickedYet(Warehouse* warehouse);

	// Function that updates the order that will be served next
	void updateOrderBeingServedNext();

	// Function that saves a route to the list of routes
	void saveRoute(int startTime, int arrivalTime, double fromLat, double fromLon, double toLat, double toLon);

	// Functions that writes routes/orders and costs to file
	void writeRoutesAndOrdersToFile(std::string fileNameRoutes, std::string fileNameOrders);
	void writeCostsToFile(std::vector<float> costs, std::vector<float> averageRejectionRateVector, float lambdaTemporal, float lambdaSpatial, bool is_training);
	void writeStatsToFile(std::vector<float> costs, std::vector<float> averageRejectionRateVector, std::vector<float> averageWaitingTime, std::vector<float> maxWaitingTime);
	void writeMatrixToFile(std::vector<std::vector<double>> matrix, std::string filename);
	// Function to draw an inter arrival time based on rate specified in data
	int drawFromExponentialDistribution(double lambda);

	// Function that returns the objective value (waiting time + penalty)
	int getObjValue();


	// Function that returns the state as a tensor
	torch::Tensor getStateAssignmentProblem(Order* order);
};


#endif