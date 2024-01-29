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
	//hourlyArrivalRates = {23,24,17,14,13,12,15,18,21,24,24,230,240,17,14,13,12,15,18,21,24,24};
	//hourlyArrivalRates = {40,38,34,30,26,23,21,22,23,24,24,20,18,19,16,14,15,16,19,22,23};
	//hourlyArrivalRates = {24,25,24,14,19,24,23,24,14,12,18,24,24};
	hourlyArrivalRates = {18,18,18,18,17,16,15,14,15,16,17,18,18,18};
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
			else if (content == "WAREHOUSE_SECTION")
				{
					// Reading warehouse data
					for (int i = 0; i < nbWarehouses; i++)
					{
						inputFile >> paramWarehouses[i].wareID >> paramWarehouses[i].lon >> paramWarehouses[i].lat >> paramWarehouses[i].initialNbCouriers >> paramWarehouses[i].initialNbPickers;
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
						inputFile >> paramClients[i].clientID >> paramClients[i].lon >> paramClients[i].lat;
						paramClients[i].nbOrders.assign(hourlyArrivalRates.size(), 0);
						paramClients[i].nbRejected.assign(hourlyArrivalRates.size(), 0);
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


