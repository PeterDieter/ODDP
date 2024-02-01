#ifndef WHQASSIGN_H
#define WHQASSIGN_H

#include <vector>

#include <gurobi_c++.h>

#include "Data.h"
#include "Matrix.h"


// auxiliary functions
long int cputime();
void printtime(long int t);
std::string itos(int);


/********************************************************

  "Warehouse-Quadrant Assignment" class

********************************************************/


class WQAssign {
public:
	// Constructor
	WQAssign(Data*,double=900,int=2,double=1e-6,int=0);
	// Parameters:
	//  - params				Pointer to Parameters
	//  - travel_cap		max allowed travel time
	//  - balance				balancing type
	//  - eps						precision epsilon
	//  - verbose				output level
	
	// Destructor
	~WQAssign();
	
	// Initialize grid and build & solve initial model
	int init();
	
	// TODO Iteratively do stuff
	int reoptimize(int=-1);

protected:
	
	// darn warning if it's initialized first ...
	double m_traveltime_cap;
	
	// which balancing constraints are used
	// - 1	relative balance ( served_u < 4.5 * served_v )
	// - 2	absolute balance ( served_u - served_v < 1.1 * gamma ; where gamma is optimized )
	int m_balance_type;
	
	// stuff
	int m_o;			// output level (for verbose mode)
	double m_eps;	// precision epsilon
	
	// access to Params class
	Data* m_params;
	
	// Distance matrix from quadrant (farthest/mean client) clients to warehouses
	//  - full dimensions (nbWarehouses times nbQuadrants)
	//  - values depend on whether max or mean is used
	Matrix m_travelTime;
	
	// list of quadrants with at least one customer
	//  - WARNING index doesn't necessarily coincide with quadrantID
	//  - WARNING still contain quadrants removed after simulation (check member "is_served" instead)
	std::vector<Quadrant*> m_used_quads;
	
// 	// list per used Quadrant number of assigned neighbouring Quadrants
// 	//  - same dimension as m_used_quads
// 	//  - used to find quadrants to remove in calc_quadrant_to_remove()
// 	std::vector<int> m_num_active_neighbours;
	
	// list per WH with reachable quadrants
	//  - full dimension (nbWarehouses)
	//  - BUT only quadrants from "m_used_quads" are ever added (unused quadrants are never appended to lists)
	std::vector<std::vector<Quadrant*> > m_reachable_quadrants;
	
	// list per Quadrant with reachable WHs
	//  - full dimension (nbQuadrants)
	//  - BUT only quadrants from "m_used_quads" are ever considered (unused quadrants have empty lists)
	std::vector<std::vector<Warehouse*> > m_reachable_wh;
	
	
	//--------------------------------
	
	// gurobi stuff
	GRBEnv* m_env;
	GRBModel* m_model;
	GRBVar* m_assign;														// assignment variables for easier access
	int m_num_assign; 
	std::pair<int,int>* m_rev_id;								// same dimension as m_assign
	std::vector<std::vector<int> > m_whvar_id;	// list of variable indices per quadrant / full outer dimension (nbQuadrants)
	std::vector<std::vector<int> > m_quvar_id;	// list of variable indices per warehouse / full outer dimension (nbWarehouses)
	GRBLinExpr* m_serve_q;	// full dimension (nbQuadrants)
	GRBLinExpr* m_serve_w;	// full dimension (nbWarehouses)
	GRBConstr* m_served;		// full dimension (nbQuadrants)
	GRBConstr** m_balance;	// full dimension (nbWarehouses * nbWarehouses)
	
	//--------------------------------
	
	// computes initial assignment of clients to quadrants
	//  - also initializes
	//			+ quadrantID's
	//			+ quadrant corners (lat/lon)
	//			+ quadrant neighbours (set to "nullptr" on boundaries)
	//			+ clientsInQuadrant/contains_client
	//			+ containing quadrant per client
	//			+ m_used_quads
	//	- removes unreachable clients
	void initialize_quadrants();
	
	//--------------------------------
	
	// build gurobi model TODO
	void build_model();
	// Parameters:
	//  - balance_type				which balance constraints to use (0: none; 1: factor; 2: absolute minimum gamma)
	
	//--------------------------------
	
	// to save assignments per iteration (for visualization)
	int iter_cnt_write = 0;
	
	// optimize current model and write out solution
	//  - passes solution to m_params
	//  - writes out assignment to textfile
	int solve_check_writesol (bool=true, bool=true);
	// Parameters:
	//  - write_iter 	true if assignment after each (removal) iteration is written to a separate file
	//  - quad_level 	true if quad level (not customer level)
	// Return:
	//  - 0 	if a solution was found
	//  - -1 	else
	
	//--------------------------------
	
	// calculate distances WH->Quad
	//  - saves distances in m_travelTime (set to "-1" for unused quadrants)
	//  - values depend on max or mean
	void calc_dist(bool=false);
	// Parameters:
	//  - use_mean				use mean distance (maximum distance otherwise)
	
	//--------------------------------
	
	// calculate reachability
	//  - only called in init (so not again afterwards)
	//  - only based on m_travelTime which depends on max or mean (taken over customers in each quadrant)
	//    (but absolutely unrachable clients are removed in initialize_quadrants() )
	//  - update lists
	//			+ m_reachable_quadrants/m_reachable_wh
	bool calc_reachability();
	// Return:
	//  - true iff all quadrants can be reached (depending on max or mean value over clients therein)
	
	//--------------------------------
	
	// remove unreachable clients from Quadrants
	//  - only called in init (so not again afterwards)
	//  - removed clients are erased from clientsInQuadrant
	//  - no unreachable clients left after calling this function
	int remove_unreachable_clients(bool=false);
	// Parameters:
	//  - use_mean				use mean distance (to remove unreachable clients even if mean is used)
	//  - max_nb					maximum number of "active" neigbours when cosidering removal
	// Return:
	//  - number of removed clients
	
	//--------------------------------
	
	// calculates a Quadrant to remove based on the waiting time values from the simulation
	Quadrant* calc_quadrant_to_remove(bool=true,int=2);
	// Parameters:
	//  - use_mean				use mean distance (to remove unreachable clients even if mean is used)
	//  - max_nb					maximum number of "active" neigbours when cosidering removal
	// Return:
	//  - pointer to quadrant which should be removed

};

#endif
