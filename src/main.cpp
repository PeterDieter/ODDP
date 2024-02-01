#include <time.h>
#include <iostream>
#include <typeinfo>


#include "Data.h"
#include "Environment.h"
#include "WQAssign.h"

int main(int argc, char * argv[])
{
	int verbose = 1;				// output level
	double eps  = 1e-6;				// precision epsilon
	double travel_cap = 900;		// cap for reachability (seconds)
  int balancing = 1;				// balancing constraitns


  // Reading the data file and initializing some data structures
  std::cout << "----- READING DATA SET " << argv[1] << " -----" << std::endl;
  Data data(argv);
  std::cout << "----- Instance with " << data.nbClients << " Clients, " << data.nbWarehouses << " Warehouses -----"<< std::endl;
  
  if (std::string(argv[3])=="staticRebalancing"){
    WQAssign wqassign(&data, travel_cap, balancing, eps, true);
    int init_mod = wqassign.init();
  }

  // Creating the Environment
  Environment environment(&data);
  environment.simulate(argv);


  // return 0 upon succesful completion
  return 0;
}