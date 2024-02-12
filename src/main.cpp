#include <time.h>
#include <iostream>
#include <typeinfo>


#include "Data.h"
#include "Environment.h"
#include "WQAssign.h"

int main(int argc, char * argv[])
{
  // Reading the data file and initializing some data structures
  std::cout << "----- READING DATA SET " << argv[1] << " -----" << std::endl;
  Data data(argv);
  std::cout << "----- Instance with " << data.nbClients << " Clients, " << data.nbWarehouses << " Warehouses -----"<< std::endl;
  
  // In case of static Rebalancing, we use Gurobi to create such a static assignment
  if (std::string(argv[3])=="staticAssignment"){
    WQAssign wqassign(&data, 900, 1, 1e-6, false);
    int init_mod = wqassign.init();
  }

  // Creating the Environment
  Environment environment(&data);
  environment.simulate(argv);


  // return 0 upon succesful completion
  return 0;
}