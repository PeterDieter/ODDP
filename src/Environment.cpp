#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <random>

#include <torch/torch.h>
#include <torch/script.h>
#include "Data.h"
#include "Matrix.h"
#include "Environment.h"



Environment::Environment(Data* data) : data(data)
{   
    std::cout<<"----- Create Environment -----"<<std::endl;
}

void Environment::initalizeForCostEstimation()
{
    int courierCounter = 0;
    int pickerCounter = 0;
    totalWaitingTime = 0;
    highestWaitingTimeOfAnOrder = 0;
    latestArrivalTime = 0;
    nbOrdersServed = 0;
    rejectCount = 0;
    nextOrderBeingServed = nullptr;
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
        newWarehouse->ordersNotAssignedToCourier = std::vector<Order*>(0);
        newWarehouse->ordersAssigned = std::vector<Order*>(0);
        newWarehouse->costsIncurred = 0;

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
        newWarehouse->ordersNotAssignedToCourier = std::vector<Order*>(0);
        newWarehouse->ordersAssigned = std::vector<Order*>(0);
        newWarehouse->costsIncurred = 0;

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
    orderTimes = std::vector<int>(0);
    clientsVector = std::vector<int>(0);
    timesToComission = std::vector<int>(0);
    timesToServe = std::vector<int>(0);
    int currTime = 0;
    int nextTime;
    while (currTime < data->simulationTime*3600){
        nextTime = drawFromExponentialDistribution(data->interArrivalTime);
        currTime += nextTime;
        orderTimes.push_back(nextTime);
        clientsVector.push_back(data->rng() % data->nbClients);
        timesToComission.push_back(drawFromExponentialDistribution(data->meanCommissionTime));
        timesToServe.push_back(drawFromExponentialDistribution(data->meanServiceTimeAtClient));
    }
}
 

void Environment::initOrder(int currentTime, Order* o)
{
    o->orderID = orders.size();
    o->timeToComission = timesToComission[o->orderID]; // Follows expoential distribution
    o->assignedCourier = nullptr;
    o->assignedPicker = nullptr;
    o->assignedWarehouse = nullptr;
    o->client = &data->paramClients[clientsVector[o->orderID]];
    o->orderTime = currentTime;
    o->arrivalTime = 0;
    o->serviceTimeAtClient = timesToServe[o->orderID]; // Follows expoential distribution
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

void Environment::choosePickerForOrder(Order* newOrder) 
{
    // We choose the picker who is available fastest
    newOrder->assignedPicker = getFastestAvailablePicker(newOrder->assignedWarehouse);
    // We set the time the picker is available again to the maximum of either the previous availability time or the current time, plus the time needed to comission the order
    newOrder->assignedPicker->timeWhenAvailable = std::max(newOrder->assignedPicker->timeWhenAvailable, currentTime) + newOrder->timeToComission;
    newOrder->donePickingTime = newOrder->assignedPicker->timeWhenAvailable;
}

void Environment::chooseCourierForOrder(Order* newOrder)
{
    // We choose the courier who is available fastest
    newOrder->assignedCourier = getFastestAvailableCourier(newOrder->assignedWarehouse);
    // We set the time the courier is arriving at the order to the maximum of either the current time, or the time the picker or couriers are available (comission time for picker has already been accounted for before). We then add the distance to the warehouse
    newOrder->arrivalTime = std::max(currentTime + newOrder->timeToComission, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable)) + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID);
    //newOrder->assignedCourier->assignedToOrder = newOrder;
    
    if(newOrder->arrivalTime<timeNextCourierArrivesAtOrder){
        timeNextCourierArrivesAtOrder = newOrder->arrivalTime;
        nextOrderBeingServed = newOrder;
    }

    if (latestArrivalTime < newOrder->arrivalTime){
        latestArrivalTime = newOrder->arrivalTime;
    }

    saveRoute(std::max(currentTime, std::max(newOrder->assignedCourier->timeWhenAvailable, newOrder->assignedPicker->timeWhenAvailable)), newOrder->arrivalTime, newOrder->assignedWarehouse->lat, newOrder->assignedWarehouse->lon, newOrder->client->lat, newOrder->client->lon);

    // Remove order from vector of orders that have not been assigned to a courier yet (If applicable)   
    RemoveOrderFromVector(newOrder->assignedWarehouse->ordersNotAssignedToCourier, newOrder);
    newOrder->assignedCourier->timeWhenAvailable = newOrder->arrivalTime + newOrder->serviceTimeAtClient + data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID)*2;
    //std::cout<<"Hello: "<<newOrder->orderID<<" "<<newOrder->assignedCourier->courierID<<" "<<newOrder->assignedCourier->timeWhenAvailable<<" "<<newOrder->arrivalTime<<" "<<newOrder->serviceTimeAtClient<<" "<<data->travelTime.get(newOrder->client->clientID, newOrder->assignedWarehouse->wareID)<<std::endl;
}


void Environment::chooseClosestWarehouseForCourier(Courier* courier)
{
    // Increment the number of order that have been served
    nbOrdersServed ++;
    totalWaitingTime += nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime;
    courier->assignedToWarehouse->costsIncurred += nextOrderBeingServed->arrivalTime - nextOrderBeingServed->orderTime;
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

int Environment::getObjValue(){
    int objectiveValue = 0;
    for (Order* order: orders){
        if (order->accepted){
            if (order->arrivalTime == -1){
                objectiveValue += data->maxWaiting;
            }else{
                objectiveValue += (order->arrivalTime-order->orderTime); 
                //std::cout<<order->arrivalTime<<" "<<order->orderTime<<" "<<order->arrivalTime-order->orderTime<<std::endl;
            }
        }else{
            objectiveValue += data->maxWaiting;
        }
    } 
    return objectiveValue;
}

void Environment::writeCostsToFile(std::vector<float> costs, std::vector<float> averageRejectionRateVector, float lambdaTemporal, float lambdaSpatial, bool is_training){
    std::string fileName;
    if (is_training){
        fileName = "data/experimentData/trainingData/averageCosts_" + std::to_string(data->maxWaiting) + "_" + std::to_string(data->interArrivalTime) + "_" + std::to_string(lambdaTemporal) + "_" + std::to_string(lambdaSpatial) +".txt";
    }
    else{
        fileName = "data/experimentData/testData/averageCosts_" + std::to_string(data->maxWaiting) + "_" + std::to_string(data->interArrivalTime) + "_"  + std::to_string(lambdaTemporal) + "_" + std::to_string(lambdaSpatial) +".txt";
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
    fileName = "data/experimentData/trainingData/statsData_" + std::to_string(data->maxWaiting) + "_" + std::to_string(data->interArrivalTime) +".txt";
   
 
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

void Environment::RemoveCourierFromVector(std::vector<Courier*> & V, Courier* courierToDelete) {
    V.erase(
        std::remove_if(V.begin(), V.end(), [&](Courier* const & o) {
            return o->courierID == courierToDelete->courierID;
        }),
        V.end());
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

int Environment::getNumberOfAvailablePickers(Warehouse* war){
    int availablePickers = 0;
    for (auto picker : war->pickersAssigned){
        if(picker->timeWhenAvailable < currentTime){
            availablePickers += 1;
        }
    }
    return availablePickers;
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

int Environment::getNumberOfOrdersNotPickedYet(Warehouse* w)
{
    int counter = 0;
    for (Order* o : w->ordersAssigned){
        if (o->donePickingTime>currentTime){
            counter += 1;
        }
    }
    return counter;
}

torch::Tensor Environment::getState(Order* order){
    // For now, the state is only the distances to the warehouses
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(order->client->clientID);
    std::vector<float> state;
    for (int i = 0; i < distancesToWarehouses.size(); i++) {
        state.push_back(distancesToWarehouses[i]);
    }
    
    std::vector<int> wareHouseLoad;
    for (Warehouse* w : warehouses){
        state.push_back(w->couriersAssigned.size());
        state.push_back(getNumberOfAvailablePickers(w));
        state.push_back(std::max(0, getFastestAvailablePicker(w)->timeWhenAvailable - currentTime));
        state.push_back(std::max(0, getFastestAvailableCourier(w)->timeWhenAvailable - currentTime));
    }

    //state.push_back(currentTime);

    // vector to tensor
    auto options = torch::TensorOptions().dtype(at::kFloat);
    torch::Tensor stateTensor = torch::from_blob(state.data(), {1, data->nbWarehouses*5}, options).clone().to(torch::kFloat);
    return stateTensor;
}

double Environment::euclideanDistance(double latFrom, double latTo, double lonFrom, double lonTo){
    double dx = latFrom - latTo;
    double dy = lonFrom - lonTo;
    return std::sqrt(dx * dx + dy * dy)*100;
}

std::vector<int> Environment::getDonePickingTimes(Warehouse* w)
{
    std::vector<int> readyTimes;
    for(Picker* p : w->pickersAssigned)
    {
        readyTimes.push_back(p->timeWhenAvailable);
    }
    return readyTimes; 
}

std::vector<int> Environment::getDoneCourierTimes(Warehouse* w)
{
    std::vector<int> readyTimes;
    for(Courier* c : w->couriersAssigned)
    {
        readyTimes.push_back(c->timeWhenAvailable);
    }
    return readyTimes; 
}

bool Environment::isFeasible(Order* newOrder, Warehouse* warehouse){
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
    int fastestPickerTimeAvailable = getFastestAvailablePicker(warehouses[indexClosestWarehouse])->timeWhenAvailable;
    int waitingForCourierTime = getFastestAvailableCourier(warehouses[indexClosestWarehouse])->timeWhenAvailable-(currentTime + newOrder->timeToComission + std::max(0,fastestPickerTimeAvailable-currentTime));
    int waitingForPickerTime = fastestPickerTimeAvailable-currentTime;

    if (distancesToWarehouses[indexClosestWarehouse] + newOrder->timeToComission + std::max(0,waitingForPickerTime) + std::max(0,waitingForCourierTime) <= data->maxWaiting){
        return true;
        newOrder->assignedWarehouse = warehouses[indexClosestWarehouse];
        newOrder->accepted = true;
    }
    return false;
}

void Environment::chooseClosestWarehouseForOrder(Order* newOrder)
{
    // For now we just assign the order to the closest warehouse
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

double Environment::getOpportunityCostsLB(Warehouse* w, Order* o)
{
    double oppCosts = 0;
    double lambda = 1/(data->interArrivalTime*data->nbWarehouses);
    std::vector<int> donePickingTimes = getDonePickingTimes(w);
    std::vector<int> doneCourierTimes = getDoneCourierTimes(w);
    sort(donePickingTimes.begin(), donePickingTimes.end()); 
    sort(doneCourierTimes.begin(), doneCourierTimes.end()); 
    
    return oppCosts;
}


void Environment::chooseWarehouseForOrderReassignment(Order* newOrder)
{
    // For now we just assign the order to the closest warehouse
    std::vector<int> distancesToWarehouses = data->travelTime.getRow(newOrder->client->clientID);
    int indexClosestWarehouse = std::min_element(distancesToWarehouses.begin(), distancesToWarehouses.end())-distancesToWarehouses.begin();
    int warehouseCounter = 0;
    std::vector<int> indices(data->nbWarehouses);
    std::vector< int> costs(data->nbWarehouses, INT16_MAX);
    std::vector< double> expectedRejections(data->nbWarehouses, 0);
    for(Warehouse* w: warehouses){
        int numberOfOrdersNotPickedYet = getNumberOfOrdersNotPickedYet(warehouses[warehouseCounter]);
        if (numberOfOrdersNotPickedYet < warehouses[warehouseCounter]->K){
            double oppCosts = getOpportunityCostsLB(w, newOrder);
            expectedRejections[warehouseCounter] = oppCosts;
            int fastestPickerTimeAvailable = getFastestAvailablePicker(w)->timeWhenAvailable;
            int waitingForCourierTime = getFastestAvailableCourier(w)->timeWhenAvailable-(currentTime+ newOrder->timeToComission + std::max(0,fastestPickerTimeAvailable-currentTime));
            int waitingForPickerTime = fastestPickerTimeAvailable-currentTime;
            costs[warehouseCounter] = distancesToWarehouses[warehouseCounter] + newOrder->timeToComission + std::max(0,waitingForPickerTime) + std::max(0,waitingForCourierTime);
            indices[warehouseCounter] = warehouseCounter;
        }
        warehouseCounter += 1;
    }
    std::sort(indices.begin(), indices.end(),[&](int A, int B) -> bool {return costs[A] < costs[B];});

    bool assigned = false;
    for (int w : indices){
        if (costs[w] <= data->maxWaiting)
        {
            newOrder->assignedWarehouse = warehouses[w];
            newOrder->accepted = true;
            assigned = true;
            break;
        }
    }

    if (assigned == false){
        newOrder->accepted = false;
        rejectCount++;
    }

}

void Environment::chooseWarehouseForOrderNN(Order* newOrder, neuralNetwork& n){
    torch::Tensor state = getState(newOrder);
    torch::Tensor prediction = n.forward(state);
    // Prediction tensor to vector
    std::vector<float> predVector(prediction.data_ptr<float>(), prediction.data_ptr<float>() + prediction.numel());
    int indexWarehouse = std::min_element(predVector.begin(), predVector.end())-predVector.begin();
    // If the index is nb.warehouses, we reject the order
    if (indexWarehouse >= data->nbWarehouses){
        newOrder->accepted = false;
        rejectCount++;
    }else{
        if (isFeasible(newOrder, warehouses[indexWarehouse])){
            newOrder->assignedWarehouse = warehouses[indexWarehouse];
            newOrder->accepted = true;
        }else{
            newOrder->accepted = false;
            rejectCount++;
        }

    }

    if (newOrder->orderID == 0){
        assingmentProblemStates = state;
        assingmentProblemActions = torch::tensor({indexWarehouse});
    }else{
        assingmentProblemStates = torch::cat({assingmentProblemStates, state});
        assingmentProblemActions = torch::cat({assingmentProblemActions, torch::tensor({indexWarehouse})});
    }
}

void Environment::basePolicy(int policy)
{
    std::cout<<"----- Simulation starts -----"<<std::endl;
   
    double running_costs = 0.0;
    double runningCounter = 0.0;
    double runnining_rejections = 0.0;
    std::vector< float> averageCostVector;
    std::vector< float> averageRejectionRateVector;
    std::vector< float> meanWaitingTimeVector;
    std::vector< float> maxWaitingTimeVector;

    for (int epoch = 1; epoch <= 500; epoch++) {
        // Initialize data structures
        initialize();
        
        // Start with simulation
        int counter = 0;
        currentTime = 0;
        timeCustomerArrives = 0;
        timeNextCourierArrivesAtOrder = INT_MAX;
        while (currentTime < data->simulationTime*3600 || ordersAssignedToCourierButNotServed.size() > 0){
            // Keep track of current time
            if (counter == orderTimes.size()-1){
                currentTime = timeNextCourierArrivesAtOrder;
            }else{
                currentTime = std::min(timeCustomerArrives, timeNextCourierArrivesAtOrder);
            }
            if (timeCustomerArrives < timeNextCourierArrivesAtOrder && currentTime <= data->simulationTime*3600  && counter<orderTimes.size()-1){
                timeCustomerArrives += orderTimes[counter];
                currentTime = timeCustomerArrives;
                counter += 1;
                // Draw new order and assign it to warehouse, picker and courier. MUST BE IN THAT ORDER!!!
                Order* newOrder = new Order;
                initOrder(timeCustomerArrives, newOrder);
                orders.push_back(newOrder);
                // We immediately assign the order to a warehouse and a picker
                if (policy==0){
                    chooseClosestWarehouseForOrder(newOrder);
                }else if(policy == 1){
                    chooseWarehouseForOrderReassignment(newOrder);
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
        double occuredCosts = getObjValue();
        running_costs += occuredCosts;
        runningCounter += 1;
        runnining_rejections += (float)rejectCount/ orders.size();
        averageCostVector.push_back(occuredCosts);
        averageRejectionRateVector.push_back((float)rejectCount/(float)orderTimes.size());
        if (nbOrdersServed > 0){
            //std::cout<<"----- Iteration: " << epoch << " Number of orders that arrived: " << orders.size() << " and served: " << nbOrdersServed << " Obj. value: " << getObjValue() << ". Mean wt: " << totalWaitingTime/nbOrdersServed <<" seconds. Highest wt: " << highestWaitingTimeOfAnOrder <<" seconds. -----" <<std::endl;
            meanWaitingTimeVector.push_back(totalWaitingTime/nbOrdersServed);
            maxWaitingTimeVector.push_back(highestWaitingTimeOfAnOrder);
        }else{
            meanWaitingTimeVector.push_back(0);
            maxWaitingTimeVector.push_back(0);
        }

        //std::cout<<"----- Simulation finished -----"<<std::endl;
        //std::cout<<"----- Number of orders that arrived: " << orders.size() << " and served: " << nbOrdersServed << " Obj. value: " << getObjValue() << ". Mean wt: " << totalWaitingTime/nbOrdersServed <<" seconds. Highest wt: " << highestWaitingTimeOfAnOrder <<" seconds. -----" <<std::endl;
        //writeRoutesAndOrdersToFile("data/animationData/routes.txt", "data/animationData/orders.txt");
    }
    //writeStatsToFile(averageCostVector, averageRejectionRateVector, meanWaitingTimeVector, maxWaitingTimeVector);
    std::cout<< "Iterations: " << runningCounter <<" Average costs: " << running_costs / runningCounter << " Average rejection rate: " <<  runnining_rejections / runningCounter <<std::endl;

}

void Environment::trainPolicy()
{
    std::cout<<"----- Train policy starts -----"<<std::endl;
    
    auto assignmentNet = std::make_shared<neuralNetwork>(data->nbWarehouses*5, data->nbWarehouses+1);
    torch::Tensor lossAssignmentNet;
    double running_costs = 0.0;
    double runningCounter = 0.0;
    double runnining_rejections = 0.0;
    std::vector< float> averageCostVector;
    std::vector< float> averageRejectionRateVector;
    std::vector< float> meanWaitingTimeVector;
    std::vector< float> maxWaitingTimeVector;

    for (int epoch = 1; epoch <= 500; epoch++) {
        // Initialize data structures
        initialize();
        
        // Start with simulation
        int counter = 0;
        currentTime = 0;
        timeCustomerArrives = 0;
        timeNextCourierArrivesAtOrder = INT_MAX;
        while (currentTime < data->simulationTime*3600  || ordersAssignedToCourierButNotServed.size() > 0){
            // Keep track of current time
            if (counter == orderTimes.size()-1){
                currentTime = timeNextCourierArrivesAtOrder;
            }else{
                currentTime = std::min(timeCustomerArrives, timeNextCourierArrivesAtOrder);
            }
            if (timeCustomerArrives < timeNextCourierArrivesAtOrder && currentTime <= data->simulationTime*3600  && counter<orderTimes.size()-1){
                timeCustomerArrives += orderTimes[counter];
                currentTime = timeCustomerArrives;
                counter += 1;
                // Draw new order and assign it to warehouse, picker and courier. MUST BE IN THAT ORDER!!!
                Order* newOrder = new Order;
                initOrder(timeCustomerArrives, newOrder);
                orders.push_back(newOrder);
                // We immediately assign the order to a warehouse and a picker
                chooseWarehouseForOrderNN(newOrder, *assignmentNet);
                
                
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
        double occuredCosts = getObjValue();
        running_costs += occuredCosts;
        runningCounter += 1;
        runnining_rejections += (float)rejectCount/ orders.size();
        averageCostVector.push_back(occuredCosts);
        averageRejectionRateVector.push_back((float)rejectCount/(float)orderTimes.size());
        if (nbOrdersServed > 0){
            //std::cout<<"----- Iteration: " << epoch << " Number of orders that arrived: " << orders.size() << " and served: " << nbOrdersServed << " Obj. value: " << getObjValue() << ". Mean wt: " << totalWaitingTime/nbOrdersServed <<" seconds. Highest wt: " << highestWaitingTimeOfAnOrder <<" seconds. -----" <<std::endl;
            meanWaitingTimeVector.push_back(totalWaitingTime/nbOrdersServed);
            maxWaitingTimeVector.push_back(highestWaitingTimeOfAnOrder);
        }else{
            meanWaitingTimeVector.push_back(0);
            maxWaitingTimeVector.push_back(0);
        }

        //std::cout<<"----- Simulation finished -----"<<std::endl;
        //std::cout<<"----- Number of orders that arrived: " << orders.size() << " and served: " << nbOrdersServed << " Obj. value: " << getObjValue() << ". Mean wt: " << totalWaitingTime/nbOrdersServed <<" seconds. Highest wt: " << highestWaitingTimeOfAnOrder <<" seconds. -----" <<std::endl;
        //writeRoutesAndOrdersToFile("data/animationData/routes.txt", "data/animationData/orders.txt");
    }
    //writeStatsToFile(averageCostVector, averageRejectionRateVector, meanWaitingTimeVector, maxWaitingTimeVector);
    std::cout<< "Iterations: " << runningCounter <<" Average costs: " << running_costs / runningCounter << " Average rejection rate: " <<  runnining_rejections / runningCounter <<std::endl;

}


void Environment::simulateFromKOn(int k, neuralNetwork& n)
{
    std::cout<<"----- Simulate from k -----"<<std::endl;
    initalizeForCostEstimation();
    // Start with simulation
    int counter = 0;
    currentTime = 0;;
    timeCustomerArrives = 0;
    timeNextCourierArrivesAtOrder = INT_MAX;
    while (currentTime < data->simulationTime*3600  || ordersAssignedToCourierButNotServed.size() > 0){
        // Keep track of current time
        if (counter == orderTimes.size()-1){
            currentTime = timeNextCourierArrivesAtOrder;
        }else{
            currentTime = std::min(timeCustomerArrives, timeNextCourierArrivesAtOrder);
        }
        if (timeCustomerArrives < timeNextCourierArrivesAtOrder && currentTime <= data->simulationTime*3600  && counter<orderTimes.size()-1){
            timeCustomerArrives += orderTimes[counter];
            currentTime = timeCustomerArrives;
            counter += 1;
            // Draw new order and assign it to warehouse, picker and courier. MUST BE IN THAT ORDER!!!
            Order* newOrder = new Order;
            initOrder(timeCustomerArrives, newOrder);
            orders.push_back(newOrder);
            // We immediately assign the order to a warehouse and a picker
            if (counter > k){
                chooseWarehouseForOrderNN(newOrder, n);
            }else{
                newOrder->assignedWarehouse = orders[newOrder->orderID]->assignedWarehouse;
                newOrder->accepted = true;
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
}

void Environment::simulate(char *argv[])
{   
    if(std::string(argv[5])== "nearestWarehouse"){
        basePolicy(0);
    }else if(std::string(argv[5])== "reassignmentPolicy"){
        basePolicy(1);
    }else if(std::string(argv[5])== "nnPolicy"){
        trainPolicy();
    }else{
        std::cerr<<"Method: " << argv[5] << " not found."<<std::endl;
    }

}