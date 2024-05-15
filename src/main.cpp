#include <time.h>
#include <iostream>
#include <typeinfo>
#include <unordered_map>
#include <sys/time.h>

#include "Data.h"
#include "Environment.h"
#include "WQAssign.h"

double measureTime () {
	struct timeval tim; // timeval from "sys/time.h" also records microseconds ("usec")
	gettimeofday(&tim,NULL);
	return tim.tv_sec + (tim.tv_usec/1000000.0);	
}

int main(int argc, char * argv[])
{
  // Reading the data file and initializing some data structures
  std::unordered_map<std::string, std::string> arguments;

  // Parse arguments
  for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg.substr(0, 2) == "--") {  // Check if argument starts with '--'
          // Find the position of '=' to separate key and value
          size_t pos = arg.find('=');
          if (pos != std::string::npos) {
              std::string key = arg.substr(2, pos - 2);
              std::string value = arg.substr(pos + 1);
              arguments[key] = value;
          }
          // Handle case when argument has no value
          else {
              std::string key = arg.substr(2);
              arguments[key] = "";  // Set value to empty string
          }
      }
  }


  // Access the value of "city" if it exists
  if (arguments.find("AMethod") == arguments.end() || arguments.find("RMethod") == arguments.end() ||  arguments.find("instance") == arguments.end() || arguments.find("maxWaiting") == arguments.end()) {
      std::cout << "Assignment Method, Rebalancing Method, Instance, or maxWaiting not provided. Please read the ReadMe for more information." << std::endl;
      std::exit(-1);
  }

  if (arguments["AMethod"] == "w"){
    if (arguments.find("a") == arguments.end()|| arguments["a"]==""){
        std::cout << "Error: Need alpha value between 0 and 1 for the weighted assignment method (float)." << std::endl;
        std::exit(-1); 
    }
  }


  if (arguments["RMethod"] == "l"){
    if (arguments.find("beta") == arguments.end()|| arguments["beta"]==""){
        std::cout << "Error: Need beta value between 0 and 1 for the level rebalancing method (float)." << std::endl;
        std::exit(-1); 
    }
  }

    if (arguments["maxWaiting"]==""){
        std::cout << "Error: maxWaiting must be an integer above 0 (int)." << std::endl;
        std::exit(-1); 
    }


  std::cout << "----- READING DATA SET " << arguments["instance"]<< " -----" << std::endl;
  Data data(std::stoi(arguments["maxWaiting"]), arguments["instance"]);
  std::cout << "----- Instance with " << data.nbClients << " Clients, " << data.nbWarehouses << " Warehouses -----"<< std::endl;

  if (arguments["AMethod"] == "s"){
      if (arguments["instance"] == "instances/grid.txt"){
        std::cout << "Error: Cannot apply Gurobi partitioning to grid instance" << std::endl;
        std::exit(-1);
      }
       // In case of static Partitioning, we use Gurobi to create such a static partitioning
      WQAssign wqassign(&data, 900, 1, 1e-6, false);
      int init_mod = wqassign.init();
  }

  // Creating the Environment
  double time = -measureTime();
  
  Environment environment(&data);
  environment.simulate(arguments);
  
  time += measureTime();
	std::cout << "----- Time " << time << " -----" << std::endl;

  // return 0 upon succesful completion
  return 0;
}
