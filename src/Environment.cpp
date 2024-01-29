#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <random>

#include "Data.h"
#include "Matrix.h"
#include "Environment.h"



Environment::Environment(Data* data) : data(data)
{   
    std::cout<<"----- Create Environment -----"<<std::endl;
}

void Environment::initialize()
{
    
    // CONSTRUCTOR: First we initialize the environment by assigning 
    int courierCounter = 0;
    int pickerCounter = 0;
    totalWaitingTime = 0;
    highestWaitingTimeOfAnOrder = 0;
    latestArrivalTime = 0;
    nbOrdersServed = 0;
    rejectCount = 0;
    nextOrderBeingServed = nullptr;
    
    for (size_t ord=0; ord<ordersAssignedToCourierButNotServed.size(); ord++) {
			delete ordersAssignedToCourierButNotServed[ord];
	}
    ordersAssignedToCourierButNotServed = std::vector<Order*>(0);
    
    for (size_t ord=0; ord<couriers.size(); ord++) {
		delete couriers[ord];
	}
    couriers = std::vector<Courier*>(0);
	
    for (size_t ord=0; ord<pickers.size(); ord++) {
		delete pickers[ord];
	}
    pickers = std::vector<Picker*>(0);

	for (size_t ord=0; ord<orders.size(); ord++) {
		delete orders[ord];
	}
    orders = std::vector<Order*>(0);
	
    for (size_t ord=0; ord<routes.size(); ord++) {
		delete routes[ord];
	}
    routes = std::vector<Route*>(0);

    warehouses = std::vector<Warehouse*>(0);

    for (int wID = 0; wID < data->nbWarehouses; wID++)
    {
        Warehouse* newWarehouse = new Warehouse;
        warehouses.push_back(newWarehouse);
        newWarehouse->wareID = data->paramWarehouses[wID].wareID;
        newWarehouse->lat = data->paramWarehouses[wID].lat;
        newWarehouse->lon = data->paramWarehouses[wID].lon;
        newWarehouse->initialNbCouriers = data->paramWarehouses[wID].initialNbCouriers;
        newWarehouse->initialNbPickers = data->paramWarehouses[wID].initialNbPickers;
        newWarehouse->currentNbCustomers = 0;
        newWarehouse->ordersNotAssignedToCourier = std::vector<Order*>(0);
        newWarehouse->ordersAssigned = std::vector<Order*>(0);

        for (int cID = 0; cID < newWarehouse->initialNbCouriers; cID++)
        {
            Courier* newCourier = new Courier;
            newCourier->courierID = courierCounter;
            newCourier->assignedToWarehouse = warehouses[wID];
            newCourier->assignedToOrder = nullptr;
            newCourier->timeWhenAvailable = 0;
            couriers.push_back(newCourier);
            newWarehouse->couriersAssigned.push_back(newCourier);
            courierCounter ++;    
        }
        for (int pID = 0; pID < newWarehouse->initialNbPickers; pID++)
        {
            Picker* newPicker = new Picker;
            newPicker->pickerID = pickerCounter;
            newPicker->assignedToWarehouse = warehouses[wID];
            newPicker->timeWhenAvailable = 0;
            pickers.push_back(newPicker);
            newWarehouse->pickersAssigned.push_back(newPicker);
            pickerCounter ++;    
        }
    }

    // Now we draw the random numbers
    data->simulationTime = data->hourlyArrivalRates.size();
    orderTimes = std::vector<int>(0);
    clientsVector = std::vector<int>(0);
    timesToComission = std::vector<int>(0);
    timesToServe = std::vector<int>(0);
    int currTime = 0;
    int nextTime;
    while (currTime < data->simulationTime*1800){
        nextTime = drawFromExponentialDistribution(data->hourlyArrivalRates[currTime/1800]);
        currTime += nextTime;
        orderTimes.push_back(nextTime);
        clientsVector.push_back(data->rng() % data->nbClients);
        timesToComission.push_back(drawFromExponentialDistribution(data->meanCommissionTime));
        timesToServe.push_back(drawFromExponentialDistribution(data->meanServiceTimeAtClient));
    }
}
 

void Environment::initOrder(int currentTime, int id, Order* o)
{
    o->orderID = id;
    o->timeToComission = 180;//timesToComission[o->orderID]; // Follows expoential distribution
    o->assignedCourier = nullptr;
    o->assignedPicker = nullptr;
    o->assignedWarehouse = nullptr;
    o->client = &data->paramClients[clientsVector[o->orderID]];
    o->orderTime = currentTime;
    o->arrivalTime = 0;
    o->serviceTimeAtClient = 120;//timesToServe[o->orderID]; // Follows expoential distribution
}

int Environment::drawFromExponentialDistribution(double lambda)
{
    // Random number generator based on poisson process
    double lambdaInv = 1/lambda;
    std::exponential_distribution<double> exp (lambdaInv);
    return round(exp.operator() (data->rng));
}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void Environment::saveRoute(int startTime, int arrivalTime, double fromLat, double fromLon, double toLat, double toLon){
    Route* route = new Route;
    route->fromLat = fromLat; route->fromLon = fromLon; route->toLat = toLat; route->tolon = toLon;
    route->startTime = startTime;
    route->arrivalTime = arrivalTime;
    routes.push_back(route);
}

void Environment::writeRoutesAndOrdersToFile(std::string fileNameRoutes, std::string fileNameOrders){
	std::cout << "----- WRITING Routes IN : " << fileNameRoutes << " and Orders IN : " << fileNameOrders << std::endl;
	std::ofstream myfile(fileNameRoutes);
	if (myfile.is_open())
	{
		for (auto route : routes)
		{
            // Here we print the order of customers that we visit 
            myfile << route->startTime << " " << route->arrivalTime << " " << route->fromLat << " " << route->fromLon << " " << route->toLat << " " << route->tolon;
            myfile << std::endl;
		}
	}
	else std::cout << "----- IMPOSSIBLE TO OPEN: " << fileNameRoutes << std::endl;
    // Now the orders
    std::ofstream myfile2(fileNameOrders);
	if (myfile2.is_open())
	{
		for (auto order : orders)
		{
            if (order->accepted){
                if (order->arrivalTime == -1){
                    // Here we print the order of customers that we visit 
                    myfile2 << order->orderTime << " " << latestArrivalTime << " " << order->client->lat << " " << order->client->lon << " " << 1;
                    myfile2 << std::endl;
                }else{
                    myfile2 << order->orderTime << " " << order->arrivalTime << " " << order->client->lat << " " << order->client->lon << " " << 1;
                    myfile2 << std::endl;
                }
            }else{
                myfile2 << order->orderTime << " " << order->orderTime + 180<< " " << order->client->lat << " " << order->client->lon << " " << 0;
                myfile2 << std::endl;
            }
		}
	}
	else std::cout << "----- IMPOSSIBLE TO OPEN: " << fileNameOrders << std::endl;
}


void Environment::writeCostsToFile(std::vector<float> costs, std::vector<float> averageRejectionRateVector, float lambdaTemporal, float lambdaSpatial, bool is_training){
    std::string fileName;
    if (is_training){
        fileName = "data/experimentData/trainingData/averageCosts_" + std::to_string(data->maxWaiting) + "_" + std::to_string(lambdaTemporal) + "_" + std::to_string(lambdaSpatial) +".txt";
    }
    else{
        fileName = "data/experimentData/testData/averageCosts_" + std::to_string(data->maxWaiting) + "_"  + std::to_string(lambdaTemporal) + "_" + std::to_string(lambdaSpatial) +".txt";
    }
 
	std::cout << "----- WRITING COST VECTOR IN : " << fileName << std::endl;
	std::ofstream myfile(fileName);
	if (myfile.is_open())
	{
        int _i = 0;
        myfile << "TotalCosts " << "RejectionRate ";
        myfile << std::endl;
		for (auto cost : costs)
		{
            // Here we print the order of customers that we visit 
            myfile << cost << " " << averageRejectionRateVector[_i];
            myfile << std::endl;
            _i += 1;
		}
	}
	else std::cout << "----- IMPOSSIBLE TO OPEN: " << fileName << std::endl;
}


void Environment::writeStatsToFile(std::vector<float> costs, std::vector<float> averageRejectionRateVector, std::vector<float> averageWaitingTime, std::vector<float> maxWaitingTime){
    std::string fileName;
    fileName = "data/experimentData/trainingData/statsData_" + std::to_string(data->maxWaiting) +".txt";
   
 
	std::cout << "----- WRITING COST VECTOR IN : " << fileName << std::endl;
	std::ofstream myfile(fileName);
	if (myfile.is_open())
	{
        int _i = 0;
        myfile << "TotalCosts " << "RejectionRate " <<"MeanWaitingTime " << "MaxWaitingTime ";
        myfile << std::endl;
		for (auto cost : costs)
		{
            // Here we print the order of customers that we visit 
            myfile << cost << " " << averageRejectionRateVector[_i]<< " " << averageWaitingTime[_i]<< " " << maxWaitingTime[_i];
            myfile << std::endl;
            _i += 1;
		}
        myfile.close();
	}
	else std::cout << "----- IMPOSSIBLE TO OPEN: " << fileName << std::endl;
}

void Environment::writeMatrixToFile(std::vector<std::vector<double>> matrix, std::string filename) {
    std::ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (const auto& row : matrix) {
        for (auto it = row.begin(); it!=row.end(); ++it) {
            outputFile << *it;
            if (std::next(it) != row.end()){
                outputFile << " ";
            }
        }
        outputFile << std::endl;
    }

    outputFile.close();
}


void Environment::RemoveOrderFromVector(std::vector<Order*> & V, Order* orderToDelete) {
    V.erase(
        std::remove_if(V.begin(), V.end(), [&](Order* const & o) {
            return o->orderID == orderToDelete->orderID;
        }),
        V.end());
}

void Environment::AddOrderToVector(std::vector<Order*> & V, Order* orderToAdd) {
    bool inserted = false;
    for (int i = V.size() - 1; i >= 0; --i) {
        Order* obj = V[i]; // Access the object using the index
        if(obj->arrivalTime < orderToAdd->arrivalTime){
            V.insert(V.begin() + i +1, orderToAdd);
            inserted = true;
            break;
        }
    }
    if(!inserted){
       V.insert(V.begin(), orderToAdd); 
    }
}

Picker* Environment::getFastestAvailablePicker(Warehouse* war){
    int timeAvailable = INT_MAX;
    Picker* fastestAvailablePicker = war->pickersAssigned[0];
    for (auto picker : war->pickersAssigned){
        if(picker->timeWhenAvailable < timeAvailable){
            timeAvailable = picker->timeWhenAvailable;
            fastestAvailablePicker = picker;
        }
    }
    return fastestAvailablePicker;
}

Courier* Environment::getFastestAvailableCourier(Warehouse* war){
    int timeAvailable = INT_MAX;
    Courier* fastestAvailableCourier = war->couriersAssigned[0];
    for (auto courier : war->couriersAssigned){
        if(courier->timeWhenAvailable < timeAvailable){
            timeAvailable = courier->timeWhenAvailable;
            fastestAvailableCourier = courier;
        }
    }
    return fastestAvailableCourier;
}

void printVector(std::vector<Order*> & O) {
    std::cout << "[ ";
    for (Order* element : O) {
        std::cout << element->arrivalTime << " ";
    }
    std::cout << "]" << std::endl;
}

void printMatrix(const std::vector<std::vector<float>>& matrix) {
    for (const auto& row : matrix) {
        for (const float& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

void Environment::updateOrderBeingServedNext(){
    int highestArrivalTime = INT_MAX;
    if (ordersAssignedToCourierButNotServed.size() < 1){
        nextOrderBeingServed = nullptr;
        timeNextCourierArrivesAtOrder = INT_MAX;
    }else{
        nextOrderBeingServed = ordersAssignedToCourierButNotServed[0];
        timeNextCourierArrivesAtOrder = nextOrderBeingServed->arrivalTime;
    }
}

void Environment::writeOrderStatsToClients(){
    for (Order* o: orders){
        int hour = o->orderTime/1800;
        o->client->nbOrders[hour] += 1;
        if (!o->accepted){
            o->client->nbRejected[hour] += 1;
        }
    }
}

void Environment::writeClientsStatsToFile(std::string filename){
    std::cout << "----- WRITING Clients IN : " << filename<<std::endl;
	std::ofstream myfile(filename);
	if (myfile.is_open())
	{
		for (auto client : data->paramClients)
		{
            for (int hour = 0; hour<data->simulationTime; hour++){
                // Here we print the order of customers that we visit 
                myfile << client.lat << " " << client.lon << " " << hour << " " << client.nbOrders[hour] << " " << client.nbRejected[hour];
                myfile << std::endl;
            }

		}
	}
	else std::cout << "----- IMPOSSIBLE TO OPEN: " << filename << std::endl;    
    
}


bool Environment::isFeasible(Order* newOrder, Warehouse* warehouse){
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
    int fastestPickerTimeAvailable = getFastestAvailablePicker(warehouses[indexClosestWarehouse])->timeWhenAvailable;
    int waitingForCourierTime = getFastestAvailableCourier(warehouses[indexClosestWarehouse])->timeWhenAvailable-(currentTime + newOrder->timeToComission + std::max(0,fastestPickerTimeAvailable-currentTime));
    int waitingForPickerTime = fastestPickerTimeAvailable-currentTime;

    if (distancesToWarehouses[indexClosestWarehouse] + newOrder->timeToComission + std::max(0,waitingForPickerTime) + std::max(0,waitingForCourierTime) <= data->maxWaiting){
        return true;
    }
    return false;
}

void Environment::choosePickerForOrder(Order* newOrder) 
{
    // We choose the picker who is available fastest
    newOrder->assignedPicker = getFastestAvailablePicker(newOrder->assignedWarehouse);
    // We set the time the picker is available again to the maximum of either the previous availability time or the current time, plus the time needed to comission the order
    newOrder->assignedPicker->timeWhenAvailable = std::max(newOrder->assignedPicker->timeWhenAvailable, currentTime) + newOrder->timeToComission;
    newOrder->donePickingTime = newOrder->assignedPicker->timeWhenAvailable;
    newOrder->assignedWarehouse->currentNbCustomers += 1;
}

void Environment::chooseCourierForOrder(Order* newOrder)
{
    // We choose the courier who is available fastest
    newOrder->assignedCourier = getFastestAvailableCourier(newOrder->assignedWarehouse);
    // We set the time the courier is arriving at the order to the maximum of either the current time, or the time the picker or couriers are available (comission time for picker has already been accounted for before). We then add the distance to the warehouse
    newOrder->arrivalTime = std::max(currentTime + newOrder->timeToComission, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable + newOrder->timeToComission)) + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
    
    if(newOrder->arrivalTime<timeNextCourierArrivesAtOrder){
        timeNextCourierArrivesAtOrder = newOrder->arrivalTime;
        nextOrderBeingServed = newOrder;
    }

    if (latestArrivalTime < newOrder->arrivalTime){
        latestArrivalTime = newOrder->arrivalTime;
    }

    //saveRoute(std::max(currentTime, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable)), newOrder->arrivalTime, newOrder->assignedWarehouse->lat, newOrder->assignedWarehouse->lon, newOrder->client->lat, newOrder->client->lon);

    // Remove order from vector of orders that have not been assigned to a courier yet (If applicable)   
    RemoveOrderFromVector(newOrder->assignedWarehouse->ordersNotAssignedToCourier, newOrder);
    newOrder->assignedCourier->timeWhenAvailable = newOrder->arrivalTime + newOrder->serviceTimeAtClient + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
    //std::cout<<"Hello: "<<newOrder->orderID<<" "<<newOrder->assignedCourier->courierID<<" "<<newOrder->assignedCourier->timeWhenAvailable<<" "<<newOrder->arrivalTime<<" "<<newOrder->serviceTimeAtClient<<" "<<data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID)<<std::endl;
}


void Environment::chooseClosestWarehouseForCourier(Courier* courier)
{
    // Increment the number of order that have been served
    nbOrdersServed ++;
    totalWaitingTime += nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime;
    if (highestWaitingTimeOfAnOrder < nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime)
    {
        highestWaitingTimeOfAnOrder = nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime;
    }
    if (latestArrivalTime < courier->timeWhenAvailable){
        latestArrivalTime = courier->timeWhenAvailable;
    }
    
    saveRoute(nextOrderBeingServed->arrivalTime, courier->timeWhenAvailable, nextOrderBeingServed->client->lat, nextOrderBeingServed->client->lon, courier->assignedToWarehouse->lat, courier->assignedToWarehouse->lon);
    
    // Remove the order from the order that have not been served
    RemoveOrderFromVector(ordersAssignedToCourierButNotServed, nextOrderBeingServed);
    // Update the order that will be served next
    updateOrderBeingServedNext();
    courier->assignedToOrder = nullptr;
    courier->assignedToWarehouse->currentNbCustomers -= 1;
}

void Environment::chooseClosestWarehouseForOrder(Order* newOrder)
{
    // We just assign the order to the closest warehouse
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();

    if (isFeasible(newOrder, warehouses[indexClosestWarehouse])){
        newOrder->assignedWarehouse = warehouses[indexClosestWarehouse];
        newOrder->accepted = true;
    }else{
        newOrder->accepted = false;
        rejectCount++;
    }
}


void Environment::chooseWarehouseForOrderReassignment(Order* newOrder, float penaltyParameter)
{
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
    
    if (isFeasible(newOrder, warehouses[indexClosestWarehouse])){
        newOrder->assignedWarehouse = warehouses[indexClosestWarehouse];
        newOrder->accepted = true;
    }else{
        int warehouseCounter = 0;
        std::vector<int> indices(data->nbWarehouses);
        std::vector< int> costs(data->nbWarehouses, INT16_MAX);
        for(Warehouse* w: warehouses){
            int waitingForPickerTime = std::max(0,getFastestAvailablePicker(w)->timeWhenAvailable -currentTime);
            int waitingForCourierTime = std::max(0,getFastestAvailableCourier(w)->timeWhenAvailable-(currentTime+ newOrder->timeToComission + waitingForPickerTime));
            costs[warehouseCounter] = distancesToWarehouses[warehouseCounter]*penaltyParameter + newOrder->timeToComission + waitingForPickerTime + waitingForCourierTime;
            indices[warehouseCounter] = warehouseCounter;
            warehouseCounter += 1;
        }
        std::sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return costs[A] < costs[B];});

        if (costs[indices[0]] <= data->maxWaiting)
        {
            newOrder->assignedWarehouse = warehouses[indices[0]];
            newOrder->accepted = true;
        }else{
            newOrder->accepted = false;
            rejectCount++;
        }
  }

}

void Environment::basePolicy(int policy)
{
    std::cout<<"----- Simulation starts -----"<<std::endl;
   
    double running_costs = 0.0;
    double runningCounter = 0.0;
    double runnining_rejections = 0.0;
    std::vector< float> averageRejectionRateVector;
    std::vector< float> meanWaitingTimeVector;
    std::vector< float> maxWaitingTimeVector;

    for (int epoch = 1; epoch <= 5000; epoch++) {
        // Initialize data structures
        initialize();
        
        // Start with simulation
        int counter = 0;
        currentTime = 0;
        timeCustomerArrives = 0;
        timeNextCourierArrivesAtOrder = INT_MAX;
        while (currentTime < data->simulationTime*1800 || ordersAssignedToCourierButNotServed.size() > 0){
            // Keep track of current time
            if (counter == orderTimes.size()-1){
                currentTime = timeNextCourierArrivesAtOrder;
            }else{
                currentTime = std::min(timeCustomerArrives + orderTimes[counter] , timeNextCourierArrivesAtOrder);
            }
            if (timeCustomerArrives + orderTimes[counter] < timeNextCourierArrivesAtOrder && currentTime <= data->simulationTime*1800  && counter<orderTimes.size()-1){
                timeCustomerArrives += orderTimes[counter];
                currentTime = timeCustomerArrives;
                counter += 1;
                // Draw new order and assign it to warehouse, picker and courier. MUST BE IN THAT ORDER!!!
                Order* newOrder = new Order;
                initOrder(timeCustomerArrives, counter-1, newOrder);
                orders.push_back(newOrder);
                // We immediately assign the order to a warehouse and a picker
                if (policy==0){
                    chooseClosestWarehouseForOrder(newOrder);
                }else if(policy == 1){
                    chooseWarehouseForOrderReassignment(newOrder, 2);
                }
                
                if (newOrder->accepted){
                    newOrder->assignedWarehouse->ordersAssigned.push_back(newOrder);
                    choosePickerForOrder(newOrder);
                    chooseCourierForOrder(newOrder);
                    AddOrderToVector(ordersAssignedToCourierButNotServed, newOrder);
                }
            }else { // when a courier arrives at an order
                if (nextOrderBeingServed){
                    Courier* c = nextOrderBeingServed->assignedCourier;
                    // We choose a warehouse for the courier
                    chooseClosestWarehouseForCourier(c);
                }
            }

        }
        writeOrderStatsToClients();
        runningCounter += 1;
        runnining_rejections += (float)rejectCount/ orders.size();
        averageRejectionRateVector.push_back((float)rejectCount/(float)orderTimes.size());
        if (nbOrdersServed > 0){
            //std::cout<<"----- Iteration: " << epoch << " Number of orders that arrived: " << orders.size() << " and served: " << nbOrdersServed << " Obj. value: " << getObjValue() << ". Mean wt: " << totalWaitingTime/nbOrdersServed <<" seconds. Highest wt: " << highestWaitingTimeOfAnOrder <<" seconds. -----" <<std::endl;
            meanWaitingTimeVector.push_back(totalWaitingTime/nbOrdersServed);
            maxWaitingTimeVector.push_back(highestWaitingTimeOfAnOrder);
        }else{
            meanWaitingTimeVector.push_back(0);
            maxWaitingTimeVector.push_back(0);
        }

        if (epoch % 500 == 0) {

            std::cout << "[Iteration: " << epoch << "] Rejected requests:" << runnining_rejections / runningCounter << std::endl;
            runningCounter = 0.0;
            runnining_rejections = 0.0;
        }

        //std::cout<<"----- Simulation finished -----"<<std::endl;
        //std::cout<<"----- Number of orders that arrived: " << orders.size() << " and served: " << nbOrdersServed << " Obj. value: " << getObjValue() << ". Mean wt: " << totalWaitingTime/nbOrdersServed <<" seconds. Highest wt: " << highestWaitingTimeOfAnOrder <<" seconds. -----" <<std::endl;
        //writeRoutesAndOrdersToFile("data/animationData/routes.txt", "data/animationData/orders.txt");
    }
    writeClientsStatsToFile("clientStatistics_nearest.txt");
    //writeStatsToFile(averageCostVector, averageRejectionRateVector, meanWaitingTimeVector, maxWaitingTimeVector);
}

std::vector<float> softmax(const std::vector<float>& input, float tau) {
    std::vector<float> result;
    double sum_exp = 0.0;

    // Calculate the sum of exponentials of input elements
    for (float val : input) {
        sum_exp += 1/std::exp(val*tau);
    }

    // Calculate softmax probabilities
    for (float val : input) {
        result.push_back(1/(std::exp(val*tau) / sum_exp));
    }

    return result;
}

// Function to sample an element from a probability distribution
int Environment::sampleFromProbabilities(const std::vector<float>& probabilities) {
    // Use discrete_distribution for sampling
    std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
    return distribution(data->rng);
}

void Environment::tuneParameters()
{
    std::cout<<"----- Parameter tuning starts -----"<<std::endl;
   
    std::vector<float> ActionSpace{1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
    std::vector<std::vector<float>> CostToGoMatrix(data->hourlyArrivalRates.size(), std::vector<float>(ActionSpace.size(), 0.0));
    float learningRate = 0.0025;
    double runningCounter = 0.0;
    double runnining_rejections = 0.0;
    double exponentialParameter = -0.0002;

    for (int epoch = 1; epoch <= 90000; epoch++) {
        // Initialize data structures
        initialize();
        std::vector<int> rejectionsPerHour(data->hourlyArrivalRates.size(), 0);
        std::vector<int> actionsPerHour(data->hourlyArrivalRates.size(), 0);
        // Start with simulation
        int counter = 0;
        currentTime = 0;
        int currentHour = INT_MAX;
        float currentParameter = 1.0;
        timeCustomerArrives = 0;
        timeNextCourierArrivesAtOrder = INT_MAX;
        bool selectedRandom  = false;
        float tau = 10;
        while (currentTime < data->simulationTime*1800 || ordersAssignedToCourierButNotServed.size() > 0){
            if (currentTime/1800 != currentHour && currentTime < 1800*data->simulationTime){
                currentHour = currentTime/1800;
                if (data->rng()/ (UINT32_MAX + 1.0) > exp(exponentialParameter * epoch) || selectedRandom == true){
                    auto minColumnIter = std::min_element(CostToGoMatrix[currentHour].begin(), CostToGoMatrix[currentHour].end());
                    int minColumn = std::distance(CostToGoMatrix[currentHour].begin(), minColumnIter);
                    currentParameter = ActionSpace[minColumn];
                    actionsPerHour[currentHour] = minColumn;
                }else{
                    selectedRandom = true;
                    int picked_number = static_cast<int>(data->rng()/ (UINT32_MAX + 1.0) * ActionSpace.size()); 
                    currentParameter = ActionSpace[picked_number];
                    actionsPerHour[currentHour] = picked_number;
                }
            }
            // if (currentTime/1800 != currentHour && currentTime < 1800*data->simulationTime){
            //     currentHour = currentTime/1800; 
            //     std::vector<float> boltzmannVector = softmax(CostToGoMatrix[currentHour], tau);
            //     int picked_number = sampleFromProbabilities(boltzmannVector);
            //     currentParameter = ActionSpace[picked_number];
            //     actionsPerHour[currentHour] = picked_number;
            // }


            // Keep track of current time
            if (counter == orderTimes.size()-1){
                currentTime = timeNextCourierArrivesAtOrder;
            }else{
                currentTime = std::min(timeCustomerArrives + orderTimes[counter], timeNextCourierArrivesAtOrder);
            }
            if (timeCustomerArrives + orderTimes[counter] < timeNextCourierArrivesAtOrder && currentTime <= data->simulationTime*1800  && counter<orderTimes.size()-1){
                timeCustomerArrives += orderTimes[counter];
                currentTime = timeCustomerArrives;
                counter += 1;
                // Draw new order and assign it to warehouse, picker and courier. MUST BE IN THAT ORDER!!!
                Order* newOrder = new Order;
                initOrder(timeCustomerArrives, counter-1, newOrder);
                orders.push_back(newOrder);
                
                // We immediately assign the order to a warehouse and a picker
                chooseWarehouseForOrderReassignment(newOrder, currentParameter);
                
                if (newOrder->accepted){
                    newOrder->assignedWarehouse->ordersAssigned.push_back(newOrder);
                    choosePickerForOrder(newOrder);
                    chooseCourierForOrder(newOrder);
                    AddOrderToVector(ordersAssignedToCourierButNotServed, newOrder);
                }else{
                    rejectionsPerHour[currentHour] += 1;
                }
            }else { // when a courier arrives at an order
                if (nextOrderBeingServed){
                    Courier* c = nextOrderBeingServed->assignedCourier;
                    // We choose a warehouse for the courier
                    chooseClosestWarehouseForCourier(c);
                }
            }
        }

        if (epoch > 80000){
            writeOrderStatsToClients();
        }else{
            int sum;
            for (int hour = 0; hour<rejectionsPerHour.size(); hour++){
                sum = std::accumulate(rejectionsPerHour.begin() + hour, rejectionsPerHour.end(), 0);
                //sum = rejectionsPerHour[hour];
                CostToGoMatrix[hour][actionsPerHour[hour]] = (1-learningRate)*CostToGoMatrix[hour][actionsPerHour[hour]] + learningRate*sum;
            }
            tau = std::max((float)0, tau - (float)0.0001);
        }

        runningCounter += 1;
        runnining_rejections += (float)rejectCount/ orders.size();
        if (epoch % 500 == 0) {
            std::cout << "[Iteration: " << epoch << "] Rejected requests:" << runnining_rejections / runningCounter << std::endl;
            runningCounter = 0.0;
            runnining_rejections = 0.0;
            //printMatrix(CostToGoMatrix);
        }
    }
    writeClientsStatsToFile("clientStatistics_CFA.txt");
}

void Environment::simulate(char *argv[])
{   
    if(std::string(argv[3])== "nearestWarehouse"){
        basePolicy(0);
    }else if(std::string(argv[3])== "reassignment"){
        basePolicy(1);
    }else if(std::string(argv[3])== "tuneParameters"){
        tuneParameters();
    }else{
        std::cerr<<"Method: " << argv[3] << " not found."<<std::endl;
    }

}