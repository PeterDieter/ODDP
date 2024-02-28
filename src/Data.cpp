#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

#include "Data.h"
#include "Matrix.h"
#include "xorshift128.h"


int Quadrant::num_quad = 0;

// defautl constructor
Quadrant::Quadrant() {
	quadrantID = num_quad++;
	latSW = lonSW = latNE = lonNE = 0;
	assignedToWarehouse = nullptr;
	clientsInQuadrant = std::vector<Client*>(0);
	contains_client = false;
	is_served = true;
	neighbours = std::vector< Quadrant*>(4);
	num_active_neighbours = 0;
};

// constructor (WARNING unused atm)
Quadrant::Quadrant(int _ID,double _latSW,double _lonSW,double _latNE,double _lonNE,Warehouse* _WH,std::vector< Client*> _CL)
									: quadrantID(_ID), latSW(_latSW), lonSW(_lonSW), latNE(_latNE), lonNE(_lonNE), assignedToWarehouse(_WH), clientsInQuadrant(_CL) {
	num_quad++;
	contains_client = ( _CL.size() > 0 );
	is_served = true;
	neighbours = std::vector< Quadrant*>(4);
	num_active_neighbours = 0;
};

// destructor
Quadrant::~Quadrant() {
	num_quad--;										// update static counter
	clientsInQuadrant.resize(0);	// clean up (probably not necessary for vector class)
	neighbours.resize(0);
};


/***********************************************************************/


Data::Data(char * argv[])
{
	rng = XorShift128(0);
	nbClients = 0;
	nbWarehouses = 0;
	nbCouriers = 0;
	nbPickers = 0;
	maxWaiting = std::stoi(argv[2]);
	meanCommissionTime = 180;
	meanServiceTimeAtClient = 60;
	paramClients = std::vector<Client>(40000); // 40000 is an upper limit, can be increase ofc
	paramWarehouses = std::vector<Warehouse>(30); // 30 is an upper limit, can be increased ofc
	//hourlyArrivalRates = {20,20,20,20,20,20,20,20,20,20,20};
	//hourlyArrivalRates = {48,48};
	hourlyArrivalRates = {24,25,24,13,14,19,21,15,11,15,18,24,24};
	//hourlyArrivalRates = {13,12,12,13,16,18,16,12,7,15,12,11,24};
	//hourlyArrivalRates = {33,32,29,29,28,26,25,23,21,17,20,25,28,30};
	std::string content, content2, content3;
	std::ifstream inputFile(argv[1]);
	if (!inputFile) throw std::runtime_error("Could not find file instance");
	if (inputFile.is_open())
	{
		for (inputFile >> content; content != "EOF"; inputFile >> content)
		{
			if (content == "NUMBER_CLIENTS")
				{
					inputFile >> content2 >> nbClients;
				}
			else if (content == "NUMBER_WAREHOUSES")
				{
					inputFile >> content2 >> nbWarehouses;
				}
			else if (content == "MEAN_COMMISSION_TIME")
				{
					inputFile >> content2 >> meanCommissionTime;
				}
			else if (content == "MEAN_SERVICE_AT_CLIENT_TIME")
				{
					inputFile >> content2 >> meanServiceTimeAtClient;
				}
			else if (content == "GRID_SW_COORDINATES")
				{
					inputFile >> content2 >> grid.lonSW >> grid.latSW;
				}
			else if (content == "GRID_NE_COORDINATES")
				{
					inputFile >> content2 >> grid.lonNE >> grid.latNE;
				}
			else if (content == "STEPSIZE_LAT")
				{
					inputFile >> content2 >> grid.stepLat;
				}
			else if (content == "STEPSIZE_LON")
				{
					inputFile >> content2 >> grid.stepLon;
				}
			else if (content == "WAREHOUSE_SECTION")
				{
					// Reading warehouse data
					for (int i = 0; i < nbWarehouses; i++)
					{
						inputFile >> paramWarehouses[i].wareID >> paramWarehouses[i].location.lon >> paramWarehouses[i].location.lat >> paramWarehouses[i].initialNbCouriers >> paramWarehouses[i].initialNbPickers;
						nbCouriers += paramWarehouses[i].initialNbCouriers;
						nbPickers += paramWarehouses[i].initialNbPickers;
					}
					
					// Reduce the size of the vector of warehouses if possible
					paramWarehouses.resize(nbWarehouses);
				}
			else if (content == "CLIENT_SECTION")
				{
					// Reading client data
					for (int i = 0; i < nbClients; i++)
					{
						inputFile >> paramClients[i].clientID >> paramClients[i].location.lon >> paramClients[i].location.lat;
						paramClients[i].nbOrders.assign(hourlyArrivalRates.size(), 0);
						paramClients[i].waitingTimes.assign(hourlyArrivalRates.size(), 0);
					}
								// Reduce the size of the vector of clients if possible
					paramClients.resize(nbClients);
				}
			else if (content == "EDGE_WEIGHT_SECTION")
				{
					travelTime = Matrix(nbClients, nbWarehouses);
					for (int i = 0; i < nbClients; i++)
					{
						for (int j = 0; j < nbWarehouses; j++)
						{	
							// Keep track of the largest distance between two clients (or the depot)
							int cost;
							inputFile >> cost;
							travelTime.set(i, j, cost);
						}
					}
				}
		}
	}

}


