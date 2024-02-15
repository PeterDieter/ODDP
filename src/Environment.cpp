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

#define pi 3.14159265358979323846


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
        newWarehouse->ordersAssigned = std::vector<Order*>(0);

        for (int cID = 0; cID < newWarehouse->initialNbCouriers; cID++)
        {
            Courier* newCourier = new Courier;
            newCourier->courierID = courierCounter;
            newCourier->assignedToWarehouse = warehouses[wID];
            newCourier->assignedToOrders.clear();
            newCourier->assignedToOrders.resize(0);
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
    o->commissionTime = timesToComission[o->orderID]; // Follows expoential distribution
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

double toRad(double degree) {
    return degree/180 * pi;
}

void printOrders(std::vector<Order*> & O) {
    std::cout << "[ ";
    for (Order* element : O) {
        std::cout <<element->orderTime<<" "<<element->donePickingTime<<" "<<element->timeCourierLeavesToOrder<<" "<<element->arrivalTime << " ";
    }
    std::cout << "]" << std::endl;
}

void Environment::writeCourierRoutesToFile(std::string fileNameRoutes, std::string fileNameOrders){
    std::cout << "----- WRITING Routes IN : " << fileNameRoutes << " and Orders IN : " << fileNameOrders << std::endl;
	std::ofstream myfile(fileNameRoutes);
    std::ofstream myfile2(fileNameOrders);
	if (myfile.is_open() && myfile2.is_open())
	{
        for(Courier* c: couriers){
            for(int route = 0; route<c->assignedToOrders.size(); route++){
                int counter = 0;
                for(Order* o : c->assignedToOrders[route]){
                    myfile2 << o->orderTime << " " << o->arrivalTime << " " << o->client->lat << " " << o->client->lon << " " << 1;
                    myfile2 << std::endl;
                    if(counter == 0){
                        myfile << c->assignedToOrders[route][0]->timeCourierLeavesToOrder << " " << c->assignedToOrders[route][0]->arrivalTime << " " << c->assignedToOrders[route][0]->arrivalTime + c->assignedToOrders[route][0]->serviceTimeAtClient  << " " << c->assignedToWarehouse->lat << " " << c->assignedToWarehouse->lon << " " << c->assignedToOrders[route][0]->client->lat  << " " << c->assignedToOrders[route][0]->client->lon;
                        myfile << std::endl;
                    }else{
                        myfile << c->assignedToOrders[route][counter]->timeCourierLeavesToOrder << " " << c->assignedToOrders[route][counter]->arrivalTime << " " << c->assignedToOrders[route][counter]->arrivalTime + c->assignedToOrders[route][counter]->serviceTimeAtClient << " " << c->assignedToOrders[route][counter-1]->client->lat << " " << c->assignedToOrders[route][counter-1]->client->lon << " " << c->assignedToOrders[route][counter]->client->lat  << " " << c->assignedToOrders[route][counter]->client->lon;
                        myfile << std::endl;
                    }
                    counter += 1;
                }

                if (c->assignedToOrders[route].size() > 0){
                    int backAtDepotTime = c->assignedToOrders[route].back()->arrivalTime + c->assignedToOrders[route].back()->serviceTimeAtClient + data->travelTime.get(c->assignedToOrders[route].back()->client->clientID, c->assignedToWarehouse->wareID);
                    myfile << c->assignedToOrders[route].back()->arrivalTime + c->assignedToOrders[route].back()->serviceTimeAtClient << " " << backAtDepotTime << " " << backAtDepotTime << " " << c->assignedToOrders[route].back()->client->lat  << " " << c->assignedToOrders[route].back()->client->lon << " " << c->assignedToWarehouse->lat << " " << c->assignedToWarehouse->lon;
                    myfile << std::endl;
                }
            }
        }
    }else std::cout << "----- Cannot Open one of the files -----" << std::endl;
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

void printVector(std::vector<int> vec) {
    std::cout << "[ ";
    for (int element : vec) {
        std::cout << element << " ";
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
        if (o->arrivalTime-o->orderTime > data->maxWaiting){
            o->client->waitingTimes[hour] += 1; // o->arrivalTime-o->orderTime - data->maxWaiting;
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
                myfile << client.lat << " " << client.lon << " " << hour << " " << client.nbOrders[hour] << " " << client.waitingTimes[hour];
                myfile << std::endl;
            }

		}
	}
	else std::cout << "----- IMPOSSIBLE TO OPEN: " << filename << std::endl;    
    
}

double Environment::getTotalWaitingTime(){
    double sum = 0;
    for (Order* o: orders){
        if (o->arrivalTime-o->orderTime > data->maxWaiting){
            sum += o->arrivalTime-o->orderTime-data->maxWaiting;
        }
    }
    return sum/ orders.size();
}

double Environment::getTotalDelays(){
    double sum = 0;
    for (Order* o: orders){
        if (o->arrivalTime-o->orderTime > data->maxWaiting){
            sum += 1;
        }
    }
    return sum/ orders.size();
}

void Environment::choosePickerForOrder(Order* newOrder) 
{
    // We choose the picker who is available fastest
    newOrder->assignedPicker = getFastestAvailablePicker(newOrder->assignedWarehouse);
    // We set the time the picker is available again to the maximum of either the previous availability time or the current time, plus the time needed to comission the order
    newOrder->assignedPicker->timeWhenAvailable = std::max(newOrder->assignedPicker->timeWhenAvailable, currentTime) + newOrder->commissionTime;
    newOrder->donePickingTime = newOrder->assignedPicker->timeWhenAvailable;
    newOrder->assignedWarehouse->currentNbCustomers += 1;
}

int calculateDistance(double lat1, double long1, double lat2, double long2) {
    if(lat1 == lat2 && long1 == long2){
        return 0;
    }else{
        double dist;
        dist = sin(toRad(lat1)) * sin(toRad(lat2)) + cos(toRad(lat1)) * cos(toRad(lat2)) * cos(toRad(long1 - long2));
        dist = acos(dist);
        dist = 6371 * dist * 1000 / 5;
        return int(dist);
    }
}

double Environment::insertOrderToCourierCosts(Order* newOrder, Courier* courier, bool bundle){
    if (courier->assignedToOrders.size() > 0 && bundle){
        int numberOfRoutes = courier->assignedToOrders.size()-1;
        if (courier->assignedToOrders[numberOfRoutes].front()->timeCourierLeavesToOrder > currentTime + std::max(0,getFastestAvailablePicker(courier->assignedToWarehouse)->timeWhenAvailable -currentTime) + newOrder->commissionTime){
            Order* lastOrder = courier->assignedToOrders[courier->assignedToOrders.size()-1].back();
            double distance = calculateDistance(lastOrder->client->lat, lastOrder->client->lon, newOrder->client->lat, newOrder->client->lon);
            return lastOrder->arrivalTime + lastOrder->serviceTimeAtClient + distance - currentTime;
        };
    }
    return std::max(courier->timeWhenAvailable, std::max(currentTime, getFastestAvailablePicker(courier->assignedToWarehouse)->timeWhenAvailable) + newOrder->commissionTime) - currentTime + 1*data->travelTime.get(newOrder->client->clientID, courier->assignedToWarehouse->wareID);
} 

Courier* Environment::getCourierBundling(Order* newOrder, Warehouse* war){
    double costs = __DBL_MAX__;
    double costsNew;
    Courier* fastestAvailableCourier = war->couriersAssigned[0];
    for (auto courier : war->couriersAssigned){
        costsNew = insertOrderToCourierCosts(newOrder, courier, true);
        if(costsNew < costs){
            costs = costsNew;
            fastestAvailableCourier = courier;
        }
    }
    return fastestAvailableCourier;
}

int Environment::costsToWarehouse(Order* newOrder, Warehouse* war, bool bundle){
    double costs = __DBL_MAX__;
    double costsNew;
    for (auto courier : war->couriersAssigned){
        costsNew = insertOrderToCourierCosts(newOrder, courier, bundle);
        if(costsNew < costs){
            costs = costsNew;
        }
    }
    return costs;
}

double Environment::getLeavingTimeCourier(Courier* courier){
    return std::max(courier->timeWhenAvailable, courier->assignedToOrders[courier->assignedToOrders.size()-1][0]->orderTime + data->maxWaiting - data->travelTime.get(courier->assignedToOrders[courier->assignedToOrders.size()-1][0]->client->clientID, courier->assignedToWarehouse->wareID));
}

void Environment::chooseCourierForOrder(Order* newOrder, bool bundling)
{
    // We choose the courier who is available fastest
    if (bundling){
        newOrder->assignedCourier = getCourierBundling(newOrder, newOrder->assignedWarehouse);
        if (newOrder->assignedCourier->assignedToOrders.size() == 0){// First we check if the courier already has a route planned
            int earliestTimeToLeaveDepot = std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, getFastestAvailablePicker(newOrder->assignedWarehouse)->timeWhenAvailable + newOrder->commissionTime));
            newOrder->timeCourierLeavesToOrder = earliestTimeToLeaveDepot + 0;
            newOrder->arrivalTime = newOrder->timeCourierLeavesToOrder + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
            newOrder->assignedCourier->assignedToOrders.push_back({newOrder});
        }else{ // if a route is planned, we check if we can add the order to the current route
            int numberOfRoutes = newOrder->assignedCourier->assignedToOrders.size()-1;
            if (newOrder->assignedCourier->assignedToOrders[numberOfRoutes][0]->timeCourierLeavesToOrder <= currentTime + std::max(0,getFastestAvailablePicker(newOrder->assignedWarehouse)->timeWhenAvailable -currentTime) + newOrder->commissionTime){
                int earliestTimeToLeaveDepot = std::max(newOrder->orderTime+data->maxWaiting-data->travelTime.get(newOrder->client->clientID, newOrder->assignedCourier->assignedToWarehouse->wareID)-data->maxWaiting, std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, getFastestAvailablePicker(newOrder->assignedWarehouse)->timeWhenAvailable + newOrder->commissionTime)));
                newOrder->timeCourierLeavesToOrder = earliestTimeToLeaveDepot + 0;
                newOrder->arrivalTime = newOrder->timeCourierLeavesToOrder + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
                newOrder->assignedCourier->assignedToOrders.push_back({newOrder});
            }else if(newOrder->assignedCourier->assignedToOrders[numberOfRoutes].front()->timeCourierLeavesToOrder > currentTime + std::max(0,getFastestAvailablePicker(newOrder->assignedWarehouse)->timeWhenAvailable -currentTime) + newOrder->commissionTime){
                Order* lastOrder = newOrder->assignedCourier->assignedToOrders[numberOfRoutes].back();
                newOrder->arrivalTime = lastOrder->arrivalTime + lastOrder->serviceTimeAtClient + calculateDistance(lastOrder->client->lat, lastOrder->client->lon, newOrder->client->lat, newOrder->client->lon);
                newOrder->timeCourierLeavesToOrder = lastOrder->arrivalTime + lastOrder->serviceTimeAtClient;
                newOrder->assignedCourier->assignedToOrders[numberOfRoutes].push_back(newOrder);
            }
        }
    }else{
        newOrder->assignedCourier = getFastestAvailableCourier(newOrder->assignedWarehouse); 
        newOrder->assignedCourier->assignedToOrders.push_back({newOrder});
        // We set the time the courier is arriving at the order to the maximum of either the current time, or the time the picker or couriers are available (comission time for picker has already been accounted for before). We then add the distance to the warehouse
        newOrder->timeCourierLeavesToOrder = std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, getFastestAvailablePicker(newOrder->assignedWarehouse)->timeWhenAvailable + newOrder->commissionTime));
        newOrder->arrivalTime = std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, getFastestAvailablePicker(newOrder->assignedWarehouse)->timeWhenAvailable + newOrder->commissionTime)) + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
    }
    if(newOrder->arrivalTime<timeNextCourierArrivesAtOrder){
        timeNextCourierArrivesAtOrder = newOrder->arrivalTime;
        nextOrderBeingServed = newOrder;
    }

    // The time the courier is available is the arrival time at the current order, plus the service time at client plus travelling back to the depot.  
    newOrder->assignedCourier->timeWhenAvailable = newOrder->arrivalTime + newOrder->serviceTimeAtClient + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);

}


void Environment::chooseWarehouseForCourier(Courier* courier)
{
    // Increment the number of order that have been served
    nbOrdersServed ++;
    totalWaitingTime += nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime;
    
    // Remove the order from the order that have not been served
    RemoveOrderFromVector(ordersAssignedToCourierButNotServed, nextOrderBeingServed);
    // Update the order that will be served next
    updateOrderBeingServedNext();
    //courier->assignedToOrders.erase(courier->assignedToOrders.begin());
    courier->assignedToWarehouse->currentNbCustomers -= 1;
}

void Environment::chooseClosestWarehouseForOrder(Order* newOrder)
{
    // We just assign the order to the closest warehouse
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();

    newOrder->assignedWarehouse = warehouses[indexClosestWarehouse];
}

void Environment::chooseWarehouseBasedOnQuadrant(Order* newOrder)
{
    newOrder->assignedWarehouse = warehouses[newOrder->client->inQuadrant->assignedToWarehouse->wareID];
}


void Environment::chooseWarehouseForOrderReassignment(Order* newOrder, float penaltyParameter, bool bundle)
{
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int warehouseCounter = 0;
    std::vector<int> indices(data->nbWarehouses);
    std::vector< int> costs(data->nbWarehouses, INT16_MAX);
    for(Warehouse* w: warehouses){
        int waitingForPickerTime = std::max(0,getFastestAvailablePicker(w)->timeWhenAvailable -currentTime);
        int waitingForCourierTime = std::max(0,getFastestAvailableCourier(w)->timeWhenAvailable-(currentTime+ newOrder->commissionTime + waitingForPickerTime));
        costs[warehouseCounter] = costsToWarehouse(newOrder, w, bundle)*(1-penaltyParameter) + distancesToWarehouses[warehouseCounter]*(penaltyParameter); // distancesToWarehouses[warehouseCounter]*(1-penaltyParameter) + penaltyParameter*(waitingForPickerTime + waitingForCourierTime);
        indices[warehouseCounter] = warehouseCounter;
        warehouseCounter += 1;
    }
    std::sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return costs[A] < costs[B];});
    newOrder->assignedWarehouse = warehouses[indices[0]];
}

void Environment::simulation(int policy)
{
    std::cout<<"----- Simulation starts -----"<<std::endl;
   
    double running_costs = 0.0;
    double runningCounter = 0.0;
    double runnning_waiting = 0.0;
    double running_delays = 0.0;
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
        bool bundle = true;
        while (currentTime < data->simulationTime*1800 || ordersAssignedToCourierButNotServed.size() > 0){
            // Keep track of current time
            if ((size_t) counter == orderTimes.size()-1){
                currentTime = timeNextCourierArrivesAtOrder;
            }else{
                currentTime = std::min(timeCustomerArrives + orderTimes[counter] , timeNextCourierArrivesAtOrder);
            }
            if (timeCustomerArrives + orderTimes[counter] < timeNextCourierArrivesAtOrder && currentTime <= data->simulationTime*1800  && (size_t) counter<orderTimes.size()-1){
                timeCustomerArrives += orderTimes[counter];
                currentTime = timeCustomerArrives;
                counter += 1;
                // Draw new order and assign it to warehouse, courier, picker. MUST BE IN THAT ORDER!!!
                Order* newOrder = new Order;
                initOrder(timeCustomerArrives, counter-1, newOrder);
                orders.push_back(newOrder);
                // We immediately assign the order to a warehouse and a picker
                if (policy==0){
                    chooseClosestWarehouseForOrder(newOrder);
                }else if(policy == 1){
                    chooseWarehouseForOrderReassignment(newOrder, 0, bundle);
                }else if(policy == 2){
                    chooseWarehouseBasedOnQuadrant(newOrder);
                }
                
                newOrder->assignedWarehouse->ordersAssigned.push_back(newOrder);
                chooseCourierForOrder(newOrder, bundle);
                choosePickerForOrder(newOrder);

                AddOrderToVector(ordersAssignedToCourierButNotServed, newOrder);
                
            }else { // when a courier arrives at an order
                if (nextOrderBeingServed){
                    Courier* c = nextOrderBeingServed->assignedCourier;
                    // We choose a warehouse for the courier
                    chooseWarehouseForCourier(c);
                }
            }

        }
        writeOrderStatsToClients();
        runningCounter += 1;
        runnning_waiting += getTotalWaitingTime(); // (float)rejectCount/ orders.size();
        running_delays += getTotalDelays();
        averageRejectionRateVector.push_back(getTotalWaitingTime());
        if (nbOrdersServed > 0){
            //std::cout<<"----- Iteration: " << epoch << " Number of orders that arrived: " << orders.size() << " and served: " << nbOrdersServed << " Obj. value: " << getObjValue() << ". Mean wt: " << totalWaitingTime/nbOrdersServed <<" seconds. Highest wt: " << highestWaitingTimeOfAnOrder <<" seconds. -----" <<std::endl;
            meanWaitingTimeVector.push_back(totalWaitingTime/nbOrdersServed);
        }else{
            meanWaitingTimeVector.push_back(0);
        }

        if (epoch % 1000 == 0) {

            std::cout << "[Iteration: " << epoch << "] Delayed orders:" << running_delays / runningCounter << " Average delay: "<< runnning_waiting / runningCounter  << std::endl;
            runningCounter = 0.0;
            runnning_waiting = 0.0;
            running_delays = 0.0;
        }
        //writeCourierRoutesToFile("data/animationData/routes.txt", "data/animationData/orders.txt");
    }
    writeClientsStatsToFile("clientStatistics_CFA.txt");
    std::cout<<"----- Simulations finished -----"<<std::endl;
}

void Environment::simulate(char *argv[])
{   
    if(std::string(argv[3])== "nearestWarehouse"){
        simulation(0);
    }else if(std::string(argv[3])== "reassignment"){
        simulation(1);
    }else if (std::string(argv[3])=="staticAssignment"){
        simulation(2);
    }else{
        std::cerr<<"Method: " << argv[3] << " not found."<<std::endl;
    }

}