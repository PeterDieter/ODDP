#include <time.h>
#include <iostream>
#include <typeinfo>
#include <unordered_map>


#include "Data.h"
#include "Environment.h"
#include "WQAssign.h"

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
  if (arguments.find("method") == arguments.end() || arguments.find("instance") == arguments.end() || arguments.find("maxWaiting") == arguments.end()) {
      std::cout << "Method, Instance, or maxWaiting not provided. Please read the ReadMe for more information." << std::endl;
      std::exit(-1);
  }

  std::cout << "----- READING DATA SET " << arguments["instance"]<< " -----" << std::endl;
  Data data(std::stoi(arguments["maxWaiting"]), arguments["instance"]);
  std::cout << "----- Instance with " << data.nbClients << " Clients, " << data.nbWarehouses << " Warehouses -----"<< std::endl;

  if (arguments["method"] == "s"){
      if (arguments["instance"] == "instances/grid.txt"){
        std::cout << "Error: Cannot apply Gurobi partitioning to grid instance" << std::endl;
        std::exit(-1);
      }
       // In case of static Partitioning, we use Gurobi to create such a static partitioning
      WQAssign wqassign(&data, 900, 1, 1e-6, false);
      int init_mod = wqassign.init();
  }

  // Creating the Environment
  Environment environment(&data);
  environment.simulate(arguments);


  // return 0 upon succesful completion
  return 0;
}