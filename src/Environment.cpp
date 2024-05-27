#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <random>
#include <filesystem>
#include <unordered_map>

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
    bundledOrders = 0;
    nbOrdersNotAssignedToNearest = 0;
    nbRebalanced = 0;
    nextOrderBeingServed = nullptr;
    timeStepSize = 3600;
    currentTime = 0;
    timeCustomerArrives = 0;
    timeNextCourierArrivesAtOrder = INT_MAX;

    for (size_t ord=0; ord<ordersAssignedToCourierButNotServed.size(); ord++) {
			delete ordersAssignedToCourierButNotServed[ord];
	}
    ordersAssignedToCourierButNotServed = std::vector<Order*>(0);
		
		if ( ordersPending.size() > 0 ) {
			ordersPending.clear();
		}
    
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
        newWarehouse->location.lat = data->paramWarehouses[wID].location.lat;
        newWarehouse->location.lon = data->paramWarehouses[wID].location.lon;
        newWarehouse->initialNbCouriers = data->paramWarehouses[wID].initialNbCouriers;
        newWarehouse->initialNbPickers = data->paramWarehouses[wID].initialNbPickers;
        newWarehouse->servedOrders = 0;
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
    while (currTime < data->simulationTime*timeStepSize){
        //std::cout<<data->hourlyArrivalRates[currTime/timeStepSize]*2.5<<std::endl;
        nextTime = drawFromExponentialDistribution(data->hourlyArrivalRates[currTime/timeStepSize]);
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
    o->commissionTime = drawFromNormalDistribution(180,30);//180;//timesToComission[o->orderID]; // Follows expoential distribution
    o->assignedCourier = nullptr;
    o->assignedPicker = nullptr;
    o->assignedWarehouse = nullptr;
    o->client = &data->paramClients[clientsVector[o->orderID]];
    o->orderTime = currentTime;
    o->arrivalTime = 0;
    o->serviceTimeAtClient = drawFromNormalDistribution(120,20);//120;//timesToServe[o->orderID]; // Follows expoential distribution
}

int Environment::drawFromExponentialDistribution(double lambda)
{
    // Random number generator based on poisson process
    double lambdaInv = 1/lambda;
    std::exponential_distribution<double> exp (lambdaInv);
    return round(exp.operator() (data->rng));
}


int Environment::drawFromNormalDistribution(float mean, float stdv){
    std::normal_distribution d{mean, stdv};


    // draw a sample from the normal distribution and round it to an integer
    return int(d(data->rng));
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
                    myfile2 << o->orderTime << " " << o->arrivalTime << " " << o->client->location.lat << " " << o->client->location.lon << " " << 1;
                    myfile2 << std::endl;
                    if(counter == 0){
                        myfile << c->assignedToOrders[route][0]->timeCourierLeavesToOrder << " " << c->assignedToOrders[route][0]->arrivalTime << " " << c->assignedToOrders[route][0]->arrivalTime + c->assignedToOrders[route][0]->serviceTimeAtClient  << " " << c->assignedToWarehouse->location.lat << " " << c->assignedToWarehouse->location.lon << " " << c->assignedToOrders[route][0]->client->location.lat  << " " << c->assignedToOrders[route][0]->client->location.lon;
                        myfile << std::endl;
                    }else{
                        myfile << c->assignedToOrders[route][counter]->timeCourierLeavesToOrder << " " << c->assignedToOrders[route][counter]->arrivalTime << " " << c->assignedToOrders[route][counter]->arrivalTime + c->assignedToOrders[route][counter]->serviceTimeAtClient << " " << c->assignedToOrders[route][counter-1]->client->location.lat << " " << c->assignedToOrders[route][counter-1]->client->location.lon << " " << c->assignedToOrders[route][counter]->client->location.lat  << " " << c->assignedToOrders[route][counter]->client->location.lon;
                        myfile << std::endl;
                    }
                    counter += 1;
                }

                if (c->assignedToOrders[route].size() > 0){
                    int backAtDepotTime = c->assignedToOrders[route].back()->arrivalTime + c->assignedToOrders[route].back()->serviceTimeAtClient + data->travelTime.get(c->assignedToOrders[route].back()->client->clientID, c->assignedToWarehouse->wareID);
                    myfile << c->assignedToOrders[route].back()->arrivalTime + c->assignedToOrders[route].back()->serviceTimeAtClient << " " << backAtDepotTime << " " << backAtDepotTime << " " << c->assignedToOrders[route].back()->client->location.lat  << " " << c->assignedToOrders[route].back()->client->location.lon << " " << c->assignedToWarehouse->location.lat << " " << c->assignedToWarehouse->location.lon;
                    myfile << std::endl;
                }
            }
        }
    }else std::cout << "----- Cannot Open one of the files -----" << std::endl;
}

void Environment::writeCostsToFile(std::vector<float> costs, std::vector<float> averageDelayRateVector, float lambdaTemporal, float lambdaSpatial, bool is_training){
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
        myfile << "TotalCosts " << "DelayRate ";
        myfile << std::endl;
		for (auto cost : costs)
		{
            // Here we print the order of customers that we visit 
            myfile << cost << " " << averageDelayRateVector[_i];
            myfile << std::endl;
            _i += 1;
		}
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

void Environment::addOrderToVectorArrivalTime(std::vector<Order*> & V, Order* orderToAdd) {
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

void Environment::addOrderToVectorDecisionTime(std::vector<Order*> & V, Order* orderToAdd) {
    bool inserted = false;
    for (int i = V.size() - 1; i >= 0; --i) {
        Order* obj = V[i]; // Access the object using the index
        if(obj->decisionTime < orderToAdd->decisionTime){
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

std::tuple<int, Picker*>  Environment::getFastestAvailablePickerExcludePickers(Warehouse* war, std::vector<int> excludePickers){
    int earliestAvailable = INT_MAX;
    Picker* fastestAvailablePicker = war->pickersAssigned[0];
    for (auto picker : war->pickersAssigned){
        int timeWhenAvailable = picker->timeWhenAvailable;
        if(std::find(excludePickers.begin(), excludePickers.end(), picker->pickerID) != excludePickers.end()) {
            timeWhenAvailable += data->meanCommissionTime;
        }
        
        if(timeWhenAvailable < earliestAvailable){
            earliestAvailable = picker->timeWhenAvailable;
            fastestAvailablePicker = picker;
        }
    }
    return std::make_tuple(earliestAvailable, fastestAvailablePicker);
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

std::vector<int> Environment::printCouriersAssigned() {
    std::vector<int> nbCouriersVector(data->nbWarehouses, 0);
    for (Warehouse* w : warehouses) {
        nbCouriersVector[w->wareID] += w->couriersAssigned.size();
    }
    return nbCouriersVector;
}


void printVectorFloat(std::vector<float> vec) {
    std::cout << "[ ";
    for (float element : vec) {
        std::cout << element << " ";
    }
    std::cout << "]" << std::endl;
}

void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (const double& element : row) {
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
        int hour = o->orderTime/timeStepSize;
        o->client->nbOrders[hour] += 1;
        if (o->arrivalTime-o->orderTime > data->maxWaiting){
            o->client->timesDelayed[hour] += o->arrivalTime-o->orderTime - data->maxWaiting;
            o->client->waitingTimes[hour] += 1;
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
                myfile << client.location.lat << " " << client.location.lon << " " << hour << " " << client.nbOrders[hour] << " " << client.timesDelayed[hour]<< " " << client.waitingTimes[hour]<< " " << client.timesBundled[hour];
                myfile << std::endl;
            }

		}
	}
	else std::cout << "----- IMPOSSIBLE TO OPEN: " << filename << std::endl;    
}

void Environment::writeStatsToFile(double percDelayed, double averageDelay, double percBundled, double percReassignedOrders, double percReassignedCouriers){
    std::string bundleString = bundle ? "true" : "false";
    std::string gridInstanceString = gridInstance ? "true" : "false";
    std::string maxWaitingString = std::to_string(data->maxWaiting);
    std::string filename = "./results/" + gridInstanceString + "_" + maxWaitingString + "_" + assignmentMethod + "_" + rebalancingMethod + "_" + bundleString + "_" + alphaString + "_" + betaString + ".txt";
    bool exist = std::filesystem::exists(filename);
    //appendFileToWorkWith.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
    if (!exist) { // Write header only if file doesn't exist or is empty
        std::ofstream myfile(filename);
        myfile << "percDelayed" << ";" << "averageDelay" << ";" << "percBundled" << ";" << "percReassignedOrders" << ";" << "percReassignedCouriers" << std::endl;
        myfile << percDelayed << ";" << averageDelay << ";" << percBundled << ";" << percReassignedOrders << ";" << percReassignedCouriers << std::endl;
        myfile.close();
    }else{
        std::ofstream myfile(filename, std::ios_base::app); // Open file in append mode
        // Write statistics
        myfile << percDelayed << ";" << averageDelay << ";" << percBundled << ";" << percReassignedOrders << ";" << percReassignedCouriers << std::endl;
        myfile.close();
    }

  
}

double Environment::getAverageDelayAllCustomers(){
    double sum = 0;
    for (Order* o: orders){
        if (o->arrivalTime-o->orderTime > data->maxWaiting){
            sum += o->arrivalTime-o->orderTime-data->maxWaiting;
        }
    }
    return sum/ orders.size();
}

double Environment::getAverageDelayDelayedCustomers(){
    double sum = 0;
    int counter = 0;
    for (Order* o: orders){
        if (o->arrivalTime-o->orderTime > data->maxWaiting){
            sum += o->arrivalTime-o->orderTime-data->maxWaiting;
            counter += 1;
        }
    }
    return sum/ counter;
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

int Environment::calculateDistance(Location loc1, Location loc2){
    if(loc1.lat == loc2.lat && loc1.lon == loc2.lon){
        return 0;
    }else{
        if(!gridInstance){
            double dist;
            dist = sin(toRad(loc1.lat)) * sin(toRad(loc2.lat)) + cos(toRad(loc1.lat)) * cos(toRad(loc2.lat)) * cos(toRad(loc1.lon-loc2.lon));
            dist = acos(dist);
            dist = 6371 * dist * 1000 / 5;
            return int(dist);
        }else{
            return int(abs(loc1.lat - loc2.lat) + abs(loc1.lon - loc2.lon))*10;
        }

    }
}

void Environment::updateInformation(Order* newOrder)
{
    // First courier stuff
    if (bundle){
        if (newOrder->assignedCourier->assignedToOrders.size() == 0){// First we check if the courier already has a route planned
            int earliestTimeToLeaveDepot = std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable + newOrder->commissionTime));
            newOrder->timeCourierLeavesToOrder = earliestTimeToLeaveDepot;
            newOrder->arrivalTime = newOrder->timeCourierLeavesToOrder + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
            newOrder->assignedCourier->assignedToOrders.push_back({newOrder});
        }else{ // if a route is planned, we check if we can add the order to the current route
            int numberOfRoutes = newOrder->assignedCourier->assignedToOrders.size()-1;
            if (newOrder->assignedCourier->assignedToOrders[numberOfRoutes][0]->timeCourierLeavesToOrder <= currentTime + std::max(0,newOrder->assignedPicker->timeWhenAvailable -currentTime) + newOrder->commissionTime){
                int earliestTimeToLeaveDepot = std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable + newOrder->commissionTime));
                newOrder->timeCourierLeavesToOrder = earliestTimeToLeaveDepot;
                newOrder->arrivalTime = newOrder->timeCourierLeavesToOrder + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
                //std::cout<<currentTime<<" "<<earliestTimeToLeaveDepot<<" "<<currentTime+data->maxWaiting-data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID)<<" "<<newOrder->arrivalTime-newOrder->orderTime-data->maxWaiting<<std::endl;
                newOrder->assignedCourier->assignedToOrders.push_back({newOrder});
            }else{
                Order* lastOrder = newOrder->assignedCourier->assignedToOrders[numberOfRoutes].back();
                newOrder->arrivalTime = lastOrder->arrivalTime + lastOrder->serviceTimeAtClient + calculateDistance(lastOrder->client->location, newOrder->client->location);
                newOrder->timeCourierLeavesToOrder = lastOrder->arrivalTime + lastOrder->serviceTimeAtClient;
                newOrder->assignedCourier->assignedToOrders[numberOfRoutes].push_back(newOrder);
                int currentHour = newOrder->orderTime/timeStepSize;
                if (newOrder->assignedCourier->assignedToOrders[numberOfRoutes].size() == 1){
                    bundledOrders += 2;
                    newOrder->client->timesBundled[currentHour] += 1;
                    lastOrder->client->timesBundled[currentHour] += 1;
                }else{
                    bundledOrders += 1;
                    newOrder->client->timesBundled[currentHour] += 1;
                    
                }
 
            }
        }
    }else{
        newOrder->assignedCourier->assignedToOrders.push_back({newOrder});
        // We set the time the courier is arriving at the order to the maximum of either the current time, or the time the picker or couriers are available (comission time for picker has already been accounted for before). We then add the distance to the warehouse
        newOrder->timeCourierLeavesToOrder = std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable + newOrder->commissionTime));
        newOrder->arrivalTime = std::max(currentTime + newOrder->commissionTime, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable + newOrder->commissionTime)) + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
    }
    if(newOrder->arrivalTime<timeNextCourierArrivesAtOrder){
        timeNextCourierArrivesAtOrder = newOrder->arrivalTime;
        nextOrderBeingServed = newOrder;
    }

    // The time the courier is available is the arrival time at the current order, plus the service time at client plus travelling back to the depot.  
    newOrder->assignedCourier->timeWhenAvailable = newOrder->arrivalTime + newOrder->serviceTimeAtClient + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);

    // Now Picker stuff
    // We set the time the picker is available again to the maximum of either the previous availability time or the current time, plus the time needed to comission the order
    newOrder->assignedPicker->timeWhenAvailable = std::max(newOrder->assignedPicker->timeWhenAvailable, currentTime) + newOrder->commissionTime;
    newOrder->donePickingTime = newOrder->assignedPicker->timeWhenAvailable;
    newOrder->assignedWarehouse->servedOrders += 1;
}

double Environment::insertOrderToCourierCosts(Order* newOrder, Courier* courier, double alpha){
    if (courier->assignedToOrders.size() > 0 && bundle){
        int numberOfRoutes = courier->assignedToOrders.size()-1;
        if (courier->assignedToOrders[numberOfRoutes].front()->timeCourierLeavesToOrder > currentTime + std::max(0,getFastestAvailablePicker(courier->assignedToWarehouse)->timeWhenAvailable - currentTime) + newOrder->commissionTime){
            Order* lastOrder = courier->assignedToOrders[courier->assignedToOrders.size()-1].back();
            double distance = calculateDistance(lastOrder->client->location, newOrder->client->location);
            //if (lastOrder->arrivalTime + lastOrder->serviceTimeAtClient - currentTime + distance > data->maxWaiting){
            if (alpha == -1){
                return (lastOrder->arrivalTime + lastOrder->serviceTimeAtClient - currentTime) + distance;
            }else{
                return (lastOrder->arrivalTime + lastOrder->serviceTimeAtClient - currentTime) *(1-alpha) + distance*alpha;
            }
        };
    }
    
    if (alpha == -1){
        return (std::max(courier->timeWhenAvailable, std::max(currentTime, getFastestAvailablePicker(courier->assignedToWarehouse)->timeWhenAvailable) + newOrder->commissionTime) - currentTime) + data->travelTime.get(newOrder->client->clientID, courier->assignedToWarehouse->wareID);
    }else{
        return (std::max(courier->timeWhenAvailable, std::max(currentTime, getFastestAvailablePicker(courier->assignedToWarehouse)->timeWhenAvailable) + newOrder->commissionTime) - currentTime)*(1-alpha) + alpha* data->travelTime.get(newOrder->client->clientID, courier->assignedToWarehouse->wareID);
    }
} 


std::tuple<int, Courier*>  Environment::costsToWarehouse(Order* newOrder, Warehouse* war, double alpha){
    double costs = __DBL_MAX__;
    double costsNew;
    Courier* fastestAvailableCourier = war->couriersAssigned[0];
    for (auto courier : war->couriersAssigned){
        costsNew = insertOrderToCourierCosts(newOrder, courier, alpha);
        if(costsNew < costs){
            costs = costsNew;
            fastestAvailableCourier = courier;
        }
    }
    return std::make_tuple(costs, fastestAvailableCourier);
}


double Environment::insertOrderToCourierCostsExcludePickers(Order* newOrder, Courier* courier, int timePickerAvailable){
    if (courier->assignedToOrders.size() > 0 && bundle){
        int numberOfRoutes = courier->assignedToOrders.size()-1;
        if (courier->assignedToOrders[numberOfRoutes].front()->timeCourierLeavesToOrder > currentTime + std::max(0,timePickerAvailable - currentTime) + newOrder->commissionTime){
            Order* lastOrder = courier->assignedToOrders[courier->assignedToOrders.size()-1].back();
            double distance = calculateDistance(lastOrder->client->location, newOrder->client->location);
            return (lastOrder->arrivalTime + lastOrder->serviceTimeAtClient) *(1-0.65)+ 0.65*distance + 0*data->maxWaiting;
        };
    }
    int waitingTime = std::max(courier->timeWhenAvailable, std::max(currentTime, timePickerAvailable) + newOrder->commissionTime) *(1-0.65)+ 0.65* data->travelTime.get(newOrder->client->clientID, courier->assignedToWarehouse->wareID);
    if( waitingTime > data->maxWaiting){
        return waitingTime + 0*data->maxWaiting;
    }else{
        return waitingTime;
    }
}

std::tuple<int, Courier*>  Environment::costsToWarehouseExclude(Order* newOrder, Warehouse* war, int timePickerAvailable, std::vector<int> couriersToExclude){
    double costs = __DBL_MAX__;
    double costsNew;
    Courier* fastestAvailableCourier = war->couriersAssigned[0];
    for (auto courier : war->couriersAssigned){
        if(std::find(couriersToExclude.begin(), couriersToExclude.end(), courier->courierID) != couriersToExclude.end()) {
            continue;
        }
        costsNew = insertOrderToCourierCostsExcludePickers(newOrder, courier, timePickerAvailable);
        if(costsNew < costs){
            costs = costsNew;
            fastestAvailableCourier = courier;
        }
    }
    return std::make_tuple(costs, fastestAvailableCourier);
}

std::tuple<int, Picker*, Courier*>  Environment::postponeAssignment(Warehouse* war, std::vector<int> relatedOrders){
    int bestCosts = INT_MAX;
    Courier* bestCourier = war->couriersAssigned[0];
    Picker* bestPicker = war->pickersAssigned[0];
    do {
        int accumulatedCosts = 0;
        std::vector<int> excludedCouriers = {};
        std::vector<int> pickersToExclude = {};
        Courier* courier0;
        Picker* picker0;
        for(int el: relatedOrders){
            auto [timeWhenPickerAvailable, picker] = getFastestAvailablePickerExcludePickers(war, pickersToExclude);
            auto [cost, courier] = costsToWarehouseExclude(ordersPending[el], war, timeWhenPickerAvailable, excludedCouriers); 
            excludedCouriers.push_back(courier->courierID);
            pickersToExclude.push_back(picker->pickerID);
            accumulatedCosts += cost;
            if(el == 0){
                courier0 = courier;
                picker0 = picker;
            }
        }

        if (accumulatedCosts < bestCosts){
            bestCosts = accumulatedCosts;
            bestCourier = courier0;
            bestPicker = picker0;
        }
    }while (std::next_permutation(relatedOrders.begin(), relatedOrders.end()));
    
    return std::make_tuple(bestCosts, bestPicker, bestCourier);
}

std::tuple<int, Picker*, Courier*, Warehouse*>  Environment::postponeReAssignment(std::vector<int> relatedOrders){
    int bestCosts = INT_MAX;
    Courier* bestCourier;
    Picker* bestPicker;
    Warehouse* bestWarehouse;
    do {
        int accumulatedCosts = 0;
        std::vector<int> excludedCouriers = {};
        std::vector<int> pickersToExclude = {};
        Courier* courier0;
        Picker* picker0;
        Warehouse* warehouse0;
        for(int o: relatedOrders){
            std::vector<int> indices(data->nbWarehouses);
            std::vector< int> costs(data->nbWarehouses, INT16_MAX);
            std::vector< Courier* > bestCouriersPerWarehouse;
            std::vector< Picker* > bestPickersPerWarehouse;
            int warehouseCounter = 0;
            // For each warehouse, check the costs of assigning the order to the warehouse
            for(Warehouse* w: warehouses){
                auto [timeWhenPickerAvailable, picker] = getFastestAvailablePickerExcludePickers(w, pickersToExclude);
                auto [cost, courier] = costsToWarehouseExclude(ordersPending[o], w, timeWhenPickerAvailable, excludedCouriers);
                costs[warehouseCounter] = cost; 
                bestCouriersPerWarehouse.push_back(courier);
                bestPickersPerWarehouse.push_back(picker);
                indices[warehouseCounter] = warehouseCounter;
                warehouseCounter += 1;
            }
            std::sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return costs[A] < costs[B];});
            excludedCouriers.push_back(bestCouriersPerWarehouse[indices[0]]->courierID);
            pickersToExclude.push_back(bestPickersPerWarehouse[indices[0]]->pickerID);
            accumulatedCosts += costs[indices[0]];
            if(o == 0){
                courier0 = bestCouriersPerWarehouse[indices[0]];
                picker0 = bestPickersPerWarehouse[indices[0]];
                warehouse0 = warehouses[indices[0]];
            }
        }

        if (accumulatedCosts < bestCosts){
            bestCosts = accumulatedCosts;
            bestCourier = courier0;
            bestPicker = picker0;
            bestWarehouse = warehouse0;
        }
    }while (std::next_permutation(relatedOrders.begin(), relatedOrders.end()));
    
    return std::make_tuple(bestCosts, bestPicker, bestCourier, bestWarehouse);
}


void Environment::chooseWarehouseForCourierStatic(Courier* courier)
{
    // Increment the number of order that have been served
    nbOrdersServed ++;
    totalWaitingTime += nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime;
    int numberOfRoutes = courier->assignedToOrders.size()-1;
    int numberOfCustomers = courier->assignedToOrders[numberOfRoutes].size()-1;
    int arrivalTime = courier->assignedToOrders[numberOfRoutes][numberOfCustomers]->arrivalTime;
    // Remove the order from the order that have not been served
	ordersAssignedToCourierButNotServed.erase(ordersAssignedToCourierButNotServed.begin());
    // Update the order that will be served next
    updateOrderBeingServedNext();
}

void Environment::chooseWarehouseForCourierLevel(Courier* courier, float beta)
{
    // Increment the number of order that have been served
    int wareIDInitial = courier->assignedToWarehouse->wareID;
    nbOrdersServed ++;
    totalWaitingTime += nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime;
    int numberOfRoutes = courier->assignedToOrders.size()-1;
    int numberOfCustomers = courier->assignedToOrders[numberOfRoutes].size()-1;
    int arrivalTime = courier->assignedToOrders[numberOfRoutes][numberOfCustomers]->arrivalTime;
    float ratio = (float) courier->assignedToWarehouse->couriersAssigned.size()/courier->assignedToWarehouse->initialNbCouriers;
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(nextOrderBeingServed->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
    if(currentTime == arrivalTime && ratio > beta){
        auto it = std::find_if(courier->assignedToWarehouse->couriersAssigned.begin(), courier->assignedToWarehouse->couriersAssigned.end(), [&](const Courier* obj) {
            return obj->courierID == courier->courierID;
        });
        courier->assignedToWarehouse->couriersAssigned.erase(it);
        courier->assignedToWarehouse = warehouses[indexClosestWarehouse];
        warehouses[indexClosestWarehouse]->couriersAssigned.push_back(courier);
        courier->timeWhenAvailable = currentTime + distancesToWarehouses[indexClosestWarehouse];
    }else{
        courier->timeWhenAvailable = std::max(courier->timeWhenAvailable, currentTime + distancesToWarehouses[indexClosestWarehouse]);
    }

    // Remove the order from the order that have not been served
	ordersAssignedToCourierButNotServed.erase(ordersAssignedToCourierButNotServed.begin());
    // Update the order that will be served next
    updateOrderBeingServedNext();

    if (courier->assignedToWarehouse->wareID != wareIDInitial){
        nbRebalanced += 1;
    }
}

void Environment::chooseClosestWarehouseForOrder(Order* newOrder, std::vector<int> relatedOrders)
{
    // We just assign the order to the closest warehouse
    if (!postpone){
        std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
        int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
        newOrder->assignedWarehouse = warehouses[indexClosestWarehouse];
        auto [cost, courier] = costsToWarehouse(newOrder, newOrder->assignedWarehouse, -1);
        newOrder->assignedCourier = courier;
        newOrder->assignedPicker = getFastestAvailablePicker(newOrder->assignedWarehouse);
    }else{
        std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
        int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
        newOrder->assignedWarehouse = warehouses[indexClosestWarehouse];
        auto [cost, picker, courier] = postponeAssignment(newOrder->assignedWarehouse, relatedOrders);
        newOrder->assignedCourier = courier;
        newOrder->assignedPicker = picker;
    }

}

void Environment::chooseWarehouseBasedOnQuadrant(Order* newOrder, std::vector<int> relatedOrders)
{
    if (!postpone){
        newOrder->assignedWarehouse = warehouses[newOrder->client->inQuadrant->assignedToWarehouse->wareID];
        auto [cost, courier] = costsToWarehouse(newOrder, newOrder->assignedWarehouse, -1);
        newOrder->assignedCourier = courier;
        newOrder->assignedPicker = getFastestAvailablePicker(newOrder->assignedWarehouse);
    }else{
        newOrder->assignedWarehouse = warehouses[newOrder->client->inQuadrant->assignedToWarehouse->wareID];
        auto [cost, picker, courier] = postponeAssignment(newOrder->assignedWarehouse, relatedOrders);
        newOrder->assignedCourier = courier;
        newOrder->assignedPicker = picker;
    }
}


void Environment::chooseWarehouseForOrderReassignment(Order* newOrder, std::vector<int> relatedOrders, double alpha)
{
    if(!postpone){
        std::vector<int> indices(data->nbWarehouses);
        std::vector< int> costs(data->nbWarehouses, INT16_MAX);
        std::vector< Courier* > bestCouriersPerWarehouse;
        int warehouseCounter = 0;
        // For each warehouse, check the costs of assigning the order to the warehouse
        for(Warehouse* w: warehouses){
            auto [cost, courier] = costsToWarehouse(newOrder, w, alpha);
            costs[warehouseCounter] = cost; 
            bestCouriersPerWarehouse.push_back(courier);
            indices[warehouseCounter] = warehouseCounter;
            warehouseCounter += 1;
        }
        std::sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return costs[A] < costs[B];});
        newOrder->assignedWarehouse = warehouses[indices[0]];
        newOrder->assignedCourier = bestCouriersPerWarehouse[indices[0]];
        newOrder->assignedPicker = getFastestAvailablePicker(newOrder->assignedWarehouse);
    }else{
        auto [cost, picker, courier, warehouse] = postponeReAssignment(relatedOrders);
        newOrder->assignedWarehouse = warehouse;
        newOrder->assignedCourier = courier;
        newOrder->assignedPicker = picker;
    }
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
    if (newOrder->assignedWarehouse->wareID != warehouses[indexClosestWarehouse]->wareID){
        nbOrdersNotAssignedToNearest += 1;
    }
    
}


std::vector<int> Environment::getRelatedOrders(Order* order){
    int orderCounter = 0;
    std::vector<int> relatedOrders;
    for (Order* o : ordersPending){
        if (calculateDistance(o->client->location, order->client->location) < 1200){
            if (relatedOrders.size() < 7){
                relatedOrders.push_back(orderCounter);
            }else{
                break;
            }
        }
        orderCounter += 1;
    }
    return relatedOrders;
}

int Environment::calcNewdecisionTime(Order* newOrder){
    // For now we just take the decision immediately
    if (postpone){
        std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
        std::vector<int> indices(data->nbWarehouses) ; // vector with nbWarehouses ints.
        std::iota (std::begin(indices), std::end(indices), 0); // Fill with 0, 1, ..., nbWarehouses.
        std::sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return distancesToWarehouses[A] < distancesToWarehouses[B];});
        if(currentTime < 12*1800 && distancesToWarehouses[indices[0]]< 300){
            return currentTime + 20;
        }else{
            return currentTime;
        }
    }else{
        return currentTime;
    }
}

int Environment::calcTimeAndEvent(int count,int& event) {
	int curr_time = 0;
	// compute current time
	if ( (size_t) count == orderTimes.size()-1 ) {
		curr_time = timeNextCourierArrivesAtOrder;
		event = 2;
	}
	else{
        if (ordersPending.size() > 0){
            if(timeCustomerArrives + orderTimes[count] < ordersPending[0]->decisionTime && timeCustomerArrives + orderTimes[count] < timeNextCourierArrivesAtOrder){
                event = 0;
                curr_time = timeCustomerArrives + orderTimes[count];
            }else if (ordersPending[0]->decisionTime < timeNextCourierArrivesAtOrder){
                event = 1;
                curr_time = ordersPending[0]->decisionTime;
            }else{
                event = 2;
                curr_time = timeNextCourierArrivesAtOrder;
            }
        }else{
            if(timeCustomerArrives + orderTimes[count] < timeNextCourierArrivesAtOrder){
                event = 0;
                curr_time = timeCustomerArrives + orderTimes[count];
            }else{
                event = 2;
                curr_time = timeNextCourierArrivesAtOrder;
            }
        }
	}
    return curr_time;
}

std::vector<float> softmax(const std::vector<double>& input, float tau) {
    std::vector<float> result;
    float sum_exp = 0.0;

    // Calculate the sum of exponentials of input elements
    for (float val : input) {
        sum_exp += std::exp(-val/tau);
    }

    // Calculate softmax probabilities
    for (float val : input) {
        result.push_back(std::exp(-val/tau) / sum_exp);
    }

    return result;
}

int Environment::sampleFromProbabilities(const std::vector<float>& probabilities) {
    // Use discrete_distribution for sampling
    std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
    return distribution(data->rng);
}

void Environment::simulation(int AssignmentPolicy, int RebalancePolicy, float alpha, float beta)
{
    std::cout<<"----- Simulation starts -----"<<std::endl;
   
    double runningCounter = 0.0;
    double runnning_waiting = 0.0;
    double running_delays = 0.0;
    double running_bundling = 0.0;
    double running_rebalancingOperations = 0.0;
    double running_notNearestAssignments = 0.0;
    int totalOrdersServed = 0;
    int nbEpochs = 500;
    int print_every_x_epochs = 100;
    std::vector<int> servedVector(data->nbWarehouses, 0);
    std::vector<int> currentCourierDistribution(data->nbWarehouses, 0);
    std::vector<int> courierDistribution(data->nbWarehouses, 0);
		
	std::printf("[Iteration]\t Delayed orders\t Average delay\t Average percentage bundled\t Average percentage of reassigned orders\t Average percentage of reassigned couriers\n");
    for (int epoch = 1; epoch <= nbEpochs; epoch++) {
        // Initialize data structures
        initialize();
        // Start with simulation
        int counter = 0;
        while (currentTime < data->simulationTime*timeStepSize || ordersAssignedToCourierButNotServed.size() > 0){
			int event = -1;
			currentTime = calcTimeAndEvent(counter,event);
            if ( event == 0 ) { // Either a new order comes into the system
                timeCustomerArrives += orderTimes[counter];
                currentTime = timeCustomerArrives;
                counter += 1;
                Order* newOrder = new Order;
                initOrder(timeCustomerArrives, counter-1, newOrder);
                newOrder->decisionTime = calcNewdecisionTime(newOrder);
                orders.push_back(newOrder);
                addOrderToVectorDecisionTime(ordersPending, newOrder); // add Order to ordersPending based on the decision time determined before
			}else if (event == 1 ){ // Or we make a decision
                Order* order = ordersPending.front();
                if (AssignmentPolicy==0){
                    chooseClosestWarehouseForOrder(order, getRelatedOrders(order));
                }else if(AssignmentPolicy == 1){
                    chooseWarehouseForOrderReassignment(order, getRelatedOrders(order), -1);
                }else if(AssignmentPolicy == 2){
                    chooseWarehouseForOrderReassignment(order, getRelatedOrders(order), alpha);
                }else if(AssignmentPolicy == 3){
                    chooseWarehouseBasedOnQuadrant(order, getRelatedOrders(order));
                }
                updateInformation(order); 
                addOrderToVectorArrivalTime(ordersAssignedToCourierButNotServed, order);
                ordersPending.erase(ordersPending.begin());
            }else if ( event == 2 ) { // Or a courier arrives at an order
                if (nextOrderBeingServed){
                    Courier* c = nextOrderBeingServed->assignedCourier;
                    if (RebalancePolicy==0){
                        chooseWarehouseForCourierStatic(c);
                    }else if(RebalancePolicy == 1){
                        chooseWarehouseForCourierLevel(c, beta);
                    }
                }
            }else { // Exception Catcher
                std::cout << "Error: Unknown event" << std::endl;
                std::exit(-1);
				}

        }

        writeOrderStatsToClients();
        runningCounter += 1;
        runnning_waiting += getAverageDelayAllCustomers();
        running_delays += getTotalDelays();
        running_bundling += (float)bundledOrders/orders.size();
        running_notNearestAssignments += (float)nbOrdersNotAssignedToNearest/orders.size();
        running_rebalancingOperations += (float)nbRebalanced/orders.size();
        for(Warehouse* w: warehouses){
            servedVector[w->wareID] += w->servedOrders;
        }
        totalOrdersServed += orders.size();
        currentCourierDistribution = printCouriersAssigned();
        for (int i = 0; i < currentCourierDistribution.size(); ++i) {
            courierDistribution[i] += currentCourierDistribution[i];
        }
        
        writeStatsToFile(getTotalDelays(), getAverageDelayAllCustomers(), (float)bundledOrders/orders.size(), (float)nbOrdersNotAssignedToNearest/orders.size(), (float)nbRebalanced/orders.size());
        
        if (epoch % print_every_x_epochs == 0) {
			std::printf("[%9.0d]\t %.4f\t\t %.4f\t\t %.4f\t\t\t\t %.4f\t\t\t\t\t\t %.4f\n",epoch,running_delays / runningCounter,runnning_waiting / runningCounter,running_bundling / runningCounter, running_notNearestAssignments/runningCounter, running_rebalancingOperations/runningCounter);
            runningCounter = 0.0;
            runnning_waiting = 0.0;
            running_delays = 0.0;
            running_bundling = 0.0;
            running_rebalancingOperations = 0.0;
            running_notNearestAssignments = 0.0;
        }
        //writeCourierRoutesToFile("data/animationData/routes.txt", "data/animationData/orders.txt");
    }
   
    //writeClientsStatsToFile("clientStatistics_rebal_zip.txt");
    std::cout<<"----- Simulations finished -----"<<std::endl;
}

void Environment::simulate(std::unordered_map<std::string, std::string> arguments)
{   
    gridInstance = false;
    if (arguments["instance"].find("grid") != std::string::npos) {
        gridInstance = true;
    }

    bundle = false;
    postpone = false;
    int assignmentPolicy = 0;
    int rebalancingPolicy = 0;
    float alpha = -1;
    float beta = 1;

    if(arguments.find("b") != arguments.end()){
        std::cout<<"----- Customer Bundling: True -----"<<std::endl;
        bundle = true;
    }else{
        std::cout<<"----- Customer Bundling: False -----"<<std::endl;
    }

    if(arguments.find("a") != arguments.end()){
        alpha = std::stod(arguments["a"]);
    }

    if(arguments.find("beta") != arguments.end()){
        beta = std::stod(arguments["beta"]);
    }


    if( arguments["RMethod"] == "s"){
        std::cout<<"----- Rebalancing Method: Static -----"<<std::endl;
        rebalancingPolicy = 0;
        rebalancingMethod = "s";
    }else if (arguments["RMethod"]=="l"){
        std::cout<<"----- Rebalancing Method: Level with beta of " << beta << " -----"<<std::endl;
        rebalancingPolicy = 1;
        rebalancingMethod = "l";
    }else{
        std::cerr<<"Rebalancing Method: " << arguments["RMethod"] << " not found."<<std::endl;
        std::exit(-1);
    }

    if( arguments["AMethod"] == "n"){
        std::cout<<"----- Assignment Method: Nearest warehouse -----"<<std::endl;
        assignmentPolicy = 0;
        assignmentMethod = "n";
    }else if (arguments["AMethod"]=="r"){
        std::cout<<"----- Assignment Method: Reassignment -----"<<std::endl;
        assignmentPolicy = 1;
        assignmentMethod = "r";
    }
    else if(arguments["AMethod"]== "w"){
        std::cout<<"----- Assignment Method: Weighted Reassignment with alpha of " << alpha << " -----"<<std::endl;
        assignmentPolicy = 2;
        assignmentMethod = "w";
    }
    else if (arguments["AMethod"]=="s"){
        std::cout<<"----- Assignment Method: StaticPartitioning -----"<<std::endl;
        assignmentPolicy = 3;
        assignmentMethod = "s";
    }else{
        std::cerr<<"Assignment Method:: " << arguments["AMethod"] << " not found."<<std::endl;
        std::exit(-1);
    }

    alphaString = std::to_string(alpha);
    alphaString = alphaString.substr(0, alphaString.find('.') + 3);
    betaString = std::to_string(beta);
    betaString = betaString.substr(0, betaString.find('.') + 3);
    simulation(assignmentPolicy, rebalancingPolicy, alpha, beta);

}
