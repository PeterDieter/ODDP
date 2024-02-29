#include <stdlib.h> 	// exit(1)
#include <cmath>
#include <iostream>
#include <iomanip>		// setprecision
#include <unistd.h>		// sysconf(_SC_CLK_TCK)
#include <sys/times.h>	// times(&...)
#include <numeric>
#include <cstring>		//strcat
//#include <fstream>
//#include <sstream>

#include <gurobi_c++.h>

#include "WQAssign.h"
#include "Data.h"


/**********************************************************************************

  elementary auxiliary functions

**********************************************************************************/

long int cputime() {
	struct tms now;
	times(&now);
	return now.tms_utime;
};

void printtime(long int t) {
	std::cout << (double)t/(double)sysconf(_SC_CLK_TCK);
	//std::cout << std::setprecision(10) << (double)t/(double)sysconf(_SC_CLK_TCK);
};

std::string itos(int i) {
	std::stringstream s;
	s << i;
	return s.str();
};

/**********************************************************************************

  Warehouse-Quadrant-Assignment class

**********************************************************************************/

WQAssign::WQAssign(Data* _params, double _travel_cap, int _balance, double _eps, int _verbose)
										: m_traveltime_cap(_travel_cap), m_balance_type(_balance), m_o(_verbose), m_eps(_eps), m_params(_params), m_travelTime() {
	
	m_used_quads.resize(0); // might be redundant
	//m_num_active_neighbours.resize(0);	
	m_reachable_quadrants.resize(m_params->nbWarehouses);
	for (int i=0; i < m_params->nbWarehouses; i++) { m_reachable_quadrants[i].resize(0); }
	m_reachable_wh.resize(0); // size determined in init()
	
	// init gurobi stuff
	try{
		m_env = new GRBEnv();
		m_env->set(GRB_IntParam_OutputFlag, 0);
		m_env->set(GRB_IntParam_Threads, 1);
		m_model = new GRBModel(m_env);
		double mipGapLimit = 0.05;
        m_model->set(GRB_DoubleParam_MIPGap, mipGapLimit);
	//} catch (const GRBException* e) {
	} catch (const GRBException& e) {
		std::cout << "Gurobi error code: " << e.getErrorCode()<<std::endl;
		std::cout << "Gurobi error message: " << e.getMessage()<<std::endl;
 	}// catch (...) {
// 		std::cout << "Error in GUROBI" << std::endl;
// 	}
	
	m_assign = NULL;
	m_num_assign = 0;
	m_rev_id = NULL;
	m_serve_q = NULL;
	m_serve_w = NULL;
	m_served = NULL;
	m_balance = NULL;
};


//===============================================================================


WQAssign::~WQAssign() {
	for (int v=0; v < m_params->nbWarehouses; v++) {
		delete[] m_balance[v];
	}
	delete[] m_balance;
	delete[] m_served;
	delete[] m_serve_w;
	delete[] m_serve_q;
	delete[] m_rev_id;
	delete[] m_assign;
	delete m_model;
	delete m_env;
};


//===============================================================================


void WQAssign::calc_dist(bool use_mean) {
#ifdef DEBUG
	std::cout << "call calc_dist()" << std::endl;
#endif
	
	// calculate mean or max distances
	
	// Initialize Quadrant->Warehouse distances (farthest client or mean)
	Warehouse* curr_wh;
	Quadrant* curr_quad;
	Client* curr_client;
	
	for (int i=0; i < m_params->nbWarehouses; i++) {
		curr_wh = &m_params->paramWarehouses[i];
		for (int j=0; j < m_params->nbQuadrants; j++) {
			curr_quad = &m_params->paramQuadrants[j];
			
			// set to -1 for empty quadrants
			if ( !curr_quad->is_served || !curr_quad->contains_client ) {
				m_travelTime.set(curr_quad->quadrantID,curr_wh->wareID,-1);
				continue;
			}
			int max_time = 0;
			// find farthest/average client
			for (size_t k=0; k<curr_quad->clientsInQuadrant.size(); k++) {
				curr_client = curr_quad->clientsInQuadrant[k];
				int tmp_time = m_params->travelTime.get(curr_client->clientID,curr_wh->wareID);
				if ( use_mean ) { max_time += tmp_time; }
				else if ( tmp_time - max_time > m_eps ) { max_time = tmp_time; }
			}
			if ( use_mean ) { max_time /= curr_quad->clientsInQuadrant.size(); }
			m_travelTime.set(curr_quad->quadrantID,curr_wh->wareID,max_time);
		}
	}
	
	return;
};


//===============================================================================


bool WQAssign::calc_reachability() {
#ifdef DEBUG
	std::cout << "call calc_reachability()" << std::endl;
#endif
	
	bool ret_flag = true;
	for (int i=0; i < m_params->nbWarehouses; i++) { m_reachable_quadrants[i].resize(0); }
	for (int j=0; j < m_params->nbQuadrants; j++) { m_reachable_wh[j].resize(0); }
	
	// check if all customers in each quadrant can be reached
	Warehouse* curr_wh;
	Quadrant* curr_quad;
	std::vector<bool> quad_reached(m_used_quads.size(),false);
	for (size_t j=0; j < m_used_quads.size(); j++) {
		curr_quad = m_used_quads[j];
		bool flag = false;
		for (int i=0; i < m_params->nbWarehouses; i++) {
			curr_wh = &m_params->paramWarehouses[i];
			int tmp_time = m_travelTime.get(curr_quad->quadrantID,curr_wh->wareID);
			if ( tmp_time > m_eps && tmp_time - m_traveltime_cap < m_eps ) {
				flag = true;
				m_reachable_quadrants[curr_wh->wareID].push_back(curr_quad);
				m_reachable_wh[curr_quad->quadrantID].push_back(curr_wh);
			}
		}
		if ( flag ) { quad_reached[j] = true; }
		else {
			ret_flag = false;
			if ( m_o > 0 ) {
				std::cout << "quad " << curr_quad->quadrantID << ": with " << curr_quad->clientsInQuadrant.size() << " not reachable" << std::endl;
			}
		}
	}
	
	return ret_flag;
};

//===============================================================================


int WQAssign::remove_unreachable_clients(bool use_mean) {
#ifdef DEBUG
	std::cout << "call remove_unreachable_clients()" << std::endl;
#endif
	
	// remove unreachable clients
	
	Warehouse* curr_wh;
	Quadrant* curr_quad;
	Client* curr_client;
	
	int cnt = 0;
	std::vector<int> remo (0); // list with clientID's for removed (output only)
	
	// check all quadrants (with clients) one by one
	for (size_t j=0; j < m_used_quads.size(); j++) {
		curr_quad = m_used_quads[j];
		// skip already reachable quadrants (unless mean is used per quadrant)
		if ( m_reachable_wh[curr_quad->quadrantID].size() < 1 || use_mean ) {
			bool flag = false;
			
			while ( !flag ) {
				
				// iteratively remove problem clients until all quadrants are reachable
				std::vector<int> max_val (m_params->nbWarehouses, 0);
				std::vector<int> max_ind (m_params->nbWarehouses, -1);
				
				// check each warehouse to get farthest client
				for (int w=0; w < m_params->nbWarehouses; w++) {
					curr_wh = &m_params->paramWarehouses[w];
					int tmp_time = 0;
					int tmp_ind = -1;
					for (size_t i=0; i<curr_quad->clientsInQuadrant.size(); i++) {
						curr_client = curr_quad->clientsInQuadrant[i];
						int tmp_val = m_params->travelTime.get(curr_client->clientID,curr_wh->wareID);
						// save "worst yet" client (value and index [not the ID])
						if ( tmp_val - tmp_time > m_eps ) { 
							tmp_time = tmp_val;
							tmp_ind = i;
						}
					}
					// save "worst" client for the warehouse
					max_val[w] = tmp_time;
					max_ind[w] = tmp_ind;
				}
				
				// find minimal value among the farthest clients
				// (not important for the index which warehouse yielded the value!)
				int min_val = 1e6;
				int min_ind = -1;
				for (int w=0; w < m_params->nbWarehouses; w++) {
					if ( min_val - max_val[w] > m_eps ) {
						min_val = max_val[w];
						min_ind = max_ind[w];
					}
				}
				
				// remove from clientsInQuadrant if client is unrachable
				if ( min_val - m_traveltime_cap > m_eps ) {
					remo.push_back(curr_quad->clientsInQuadrant[min_ind]->clientID);
					curr_quad->clientsInQuadrant[min_ind]->reached = false;
					curr_quad->clientsInQuadrant.erase(curr_quad->clientsInQuadrant.begin()+min_ind);
					cnt++;
				}
				// move on to next quadrant otherwise
				else { flag = true; }
			}
		}
	}
	
	if ( m_o > 0 ) {
		std::cout << "removed_clients";
		for (size_t k=0; k<remo.size(); k++) {
			std::cout << " " << remo[k];
		}
		std::cout << std::endl;
	}
	
	return cnt;
};


//===============================================================================


void WQAssign::initialize_quadrants() {
#ifdef DEBUG
	std::cout << "call initialize_quadrants()" << std::endl;
#endif
	std::cout << "using traveltime cap of "<< m_traveltime_cap << " sec" << std::endl;	
	
	// grid data
	double step_lon = m_params->grid.stepLon;
	double step_lat = m_params->grid.stepLat;
	double latNE = m_params->grid.latNE;
	double latSW = m_params->grid.latSW;
	double lonNE = m_params->grid.lonNE;
	double lonSW = m_params->grid.lonSW;
	
	// get number of rows/columns in grid
	// latitude: southpole->-90, equator->0, northpole->+90
	int grid_hei = std::ceil((latNE - latSW) / step_lat);
	int grid_wid = std::ceil((lonNE - lonSW) / step_lon);
	
	// init quadrants
	int q_x,q_y;
	Quadrant* curr_quad;
	m_params->paramQuadrants.resize(grid_hei*grid_wid);
	m_params->nbQuadrants = grid_hei*grid_wid;
	m_reachable_wh.resize(m_params->nbQuadrants);
	for (int i=0; i < m_params->nbQuadrants; i++) { m_reachable_wh[i].resize(0); }
	// set SW,NE corners of Quadrants
	for (int i=0; i < m_params->nbQuadrants; i++) {
		curr_quad = &m_params->paramQuadrants[i];
		q_x = i % grid_wid;
		q_y = std::floor(i / grid_wid);
		curr_quad->latSW = latSW + q_y*step_lat;
		curr_quad->lonSW = lonSW + q_x*step_lon;
		curr_quad->latNE = curr_quad->latSW + step_lat;
		curr_quad->lonNE = curr_quad->lonSW + step_lon;
		
		// set neighbours
		int indN = i+grid_wid;
		int indE = i+1;
		int indS = i-grid_wid;
		int indW = i-1;
		// lower border
		if ( i < grid_wid ) curr_quad->neighbours[2] = NULL; // no neigh to the south
		else curr_quad->neighbours[2] = &m_params->paramQuadrants[indS];
		// upper border
		if ( m_params->nbQuadrants - 1 - i < grid_wid ) curr_quad->neighbours[0] = NULL; // no neigh to the north
		else curr_quad->neighbours[0] = &m_params->paramQuadrants[indN];
		// left border
		if ( i % grid_wid == 0 ) curr_quad->neighbours[3] = NULL; // no neigh to the west
		else curr_quad->neighbours[3] = &m_params->paramQuadrants[indW];
		// right border
		if ( (i - 1) % grid_wid == 0 ) curr_quad->neighbours[1] = NULL; // no neigh to the east
		else curr_quad->neighbours[1] = &m_params->paramQuadrants[indE];
	}
	
	// calculate grid numbers for each customer
	// - left to right and bottom to top
	Client* curr_client;
	double lat,lon;
	for (int i=0; i < m_params->nbClients; i++) {
		// read position
		curr_client = &m_params->paramClients[i];
		curr_client->reached = true;		// init reached (can be changed to "false" by remove_unreachable_clients() )
		lat = curr_client->location.lat;					// Latitude
		lon = curr_client->location.lon;					// Longitude
		
		// calc x/y number in grid to get quadID
		q_x = std::floor((lon - lonSW) / step_lon);
		q_y = std::floor((lat - latSW) / step_lat);
		int q_id = q_y*grid_wid + q_x;
		curr_quad = &m_params->paramQuadrants[q_id];
		
		// update bookkeeping (what a beautiful word!)
		curr_quad->contains_client = true;
		curr_quad->clientsInQuadrant.push_back(curr_client); // is initialized as an empty vector by quadrand constructor!
		curr_client->inQuadrant = curr_quad;
	}
	
	// to avoid accessing empty quadrants (faster loops but be careful with idices!)
	for (int i=0; i < m_params->nbQuadrants; i++) {
		curr_quad = &m_params->paramQuadrants[i];
		if ( curr_quad->contains_client ) { 
			m_used_quads.push_back(curr_quad);
			for (size_t nb = 0; nb < curr_quad->neighbours.size(); nb++ ) {
				Quadrant* nb_quad = curr_quad->neighbours[nb];
				if ( nb_quad ) { nb_quad->num_active_neighbours++; }
			}
		}
	}
	
	// count number of single neighbour quadrants for removal
	// TODO JUST FOR DEBUGGING PURPOSES
	int cnt_one = 0;
	for (size_t i=0; i < m_used_quads.size(); i++) {
		curr_quad = m_used_quads[i];
		int num_nb = curr_quad->num_active_neighbours;
		//if ( curr_quad->contains_client && num_nb == 1 ) {
		if ( num_nb == 1 ) {
			//std::cout << "num_neigh[" << i << "]= " << m_params->paramQuadrants[i].num_active_neighbours << std::endl;
			cnt_one++;
		}
	}
	std::cout << "num_neigh==1 " << cnt_one << std::endl;
	
	
	m_travelTime.setDimension(m_params->nbWarehouses,m_params->nbQuadrants);
	bool use_mean = true;
	// calculate mean or max distances
	calc_dist(use_mean);
	
	// check if all customers in each quadrant can be reached
	bool all_reachable = calc_reachability();
	
	// make model feasible by removing clients who cannot be reached due to the quadrants
	if ( !all_reachable || use_mean ) {
		int removed_clients = remove_unreachable_clients(use_mean);
		if ( m_o > 0 ) {
			std::cout << "removed_clients = " << removed_clients << std::endl;
		}
		
		// recalculate distances
		calc_dist(use_mean);
		
		// check of all customers in each quadrant can be reached after cleanup
		//  - not redundant since m_reachable_quadrants and m_reachable_wh are updated here
		all_reachable = calc_reachability();
	}
	
	// some output if needed
	if ( m_o > 0 ) {
		std::cout << "Grid width : " << grid_wid << std::endl;
		std::cout << "Grid height: " << grid_hei << std::endl;
		std::cout << "Used quads : " << m_used_quads.size() << " (of " << m_params->nbQuadrants << ")" << std::endl;
	}
	
	return;
};


//===============================================================================


void WQAssign::build_model() {
#ifdef DEBUG
	std::cout << "call build_model()" << std::endl;
#endif
	std::cout << "using balancing constraint type " << m_balance_type << std::endl;
	
	// timer
	long int bm_time = - cputime();
	
	Warehouse* curr_wh;
	Quadrant* curr_quad;
	
	// number of variables
	int nn = 0;
	for (int i=0; i<m_params->nbWarehouses; i++) {
		nn += m_reachable_quadrants[i].size();
	}
	m_num_assign = nn;
	
	// add variable to the model
	m_assign = new GRBVar[nn];
	
	// reverse IDs to warehouse and quadrant since variables are simply numbered 0,1,2,...
	m_rev_id = new std::pair<int,int>[nn];
	
	// save variable pointers per warehouse/quadrant for faster access
	m_whvar_id.resize(m_params->nbQuadrants);
	for (int i=0; i < m_params->nbQuadrants; i++) { m_whvar_id[i].resize(0); }
	m_quvar_id.resize(m_params->nbWarehouses);
	for (int i=0; i < m_params->nbWarehouses; i++) { m_quvar_id[i].resize(0); }
	
	m_serve_q = new GRBLinExpr[m_params->nbQuadrants];	// sum of assigned warehouses per quadrant
	m_serve_w = new GRBLinExpr[m_params->nbWarehouses];	// sum of demand * distance * assigned quadrant per warehouse
	GRBLinExpr orig_objective = 0;	// objective function (will temporarily be changed if "m_balance_type == 2" is used)
	
	int cnt = 0;
	for (int i=0; i < m_params->nbQuadrants; i++) { m_serve_q[i] = 0; }
	for (int i=0; i < m_params->nbWarehouses; i++) { m_serve_w[i] = 0; }
	
	for (int i=0; i < m_params->nbWarehouses; i++) {
		curr_wh = &m_params->paramWarehouses[i];
		int whid = curr_wh->wareID;
		for (size_t j=0; j < m_reachable_quadrants[whid].size(); j++) {
			curr_quad = m_reachable_quadrants[whid][j];
			int quid = curr_quad->quadrantID;
			m_rev_id[cnt].first = whid;
			m_rev_id[cnt].second = quid;
			
			m_whvar_id[quid].push_back(cnt);
			m_quvar_id[whid].push_back(cnt);
			
			int dem = curr_quad->clientsInQuadrant.size();
			int dist =  m_travelTime.get(quid,whid);
			m_assign[cnt] = m_model->addVar(0.0, 1.0, dem*dist, GRB_BINARY, "x_"+std::to_string(whid)+"_"+std::to_string(quid));
			orig_objective += m_assign[cnt] * dem * dist;
			m_serve_q[quid] += m_assign[cnt];
			m_serve_w[whid] += m_assign[cnt] * dem * dist;
			cnt++;
		}
	}
	
	// constraints to ensure all quadrants are m_served
	m_served = new GRBConstr[m_params->nbQuadrants];
	for (int i=0; i < m_params->nbQuadrants; i++) {
		if ( m_params->paramQuadrants[i].contains_client ) {
			m_served[i] = m_model->addConstr(m_serve_q[i] == 1.0, "served_"+std::to_string(i));
		}
	}
	
	// balancing constraints
	
	// allocate balancing constraints
	m_balance = new GRBConstr*[m_params->nbWarehouses];
	for (int v=0; v < m_params->nbWarehouses; v++) { m_balance[v] = new GRBConstr[m_params->nbWarehouses]; }
	//double max_balance = 4.5;
	//GRBVar* balance_vw = new GRBVar[m_params->nbWarehouses*(m_params->nbWarehouses-1)/2];
 	
	// choose which kind of balancing
	//  - type 1: use balancing factor
	if ( m_balance_type == 1 ) {
		
		// get upper bound on factors (used as big-M later)
		//  - bound called "UB" on maximum balancing factor
		m_model->update();
		m_model->optimize();
		int min_serve = m_model->get(GRB_DoubleAttr_ObjVal) + 1.5;
		int max_serve = 0;
		for (int v=0; v < m_params->nbWarehouses; v++) {
			int lhs = m_serve_w[v].getValue() + 0.5;
			if ( lhs < min_serve - m_eps ) min_serve = lhs;
			if ( lhs > max_serve + m_eps ) max_serve = lhs;
		}
		int ub_factor = max_serve / min_serve + 1;
		
		// add temporary balance-factor variables and rhs terms to the model
		//  - In the balancing constraints
		//       SUM_i ( x_iw  * d_i * c_iw ) <= gamma * SUM_i ( x_iw'  * d_i * c_iw' )
		//    gamma_iw replaces "gamma * x_iw'" in the rhs
		GRBVar* tmp_factors = new GRBVar[nn];
		GRBLinExpr* tmp_serve_w = new GRBLinExpr[m_params->nbWarehouses];	// sum of demand * distance * assigned quadrant per warehouse
		cnt = 0;
		for (int i=0; i < m_params->nbWarehouses; i++) {
			curr_wh = &m_params->paramWarehouses[i];
			int whid = curr_wh->wareID;
			for (size_t j=0; j < m_reachable_quadrants[whid].size(); j++) {
				curr_quad = m_reachable_quadrants[whid][j];
				int quid = curr_quad->quadrantID;
				// add variables gammi_iw
				tmp_factors[cnt] = m_model->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "gamma_"+std::to_string(whid)+"_"+std::to_string(quid));
				// add balancing constraint with artificial variables
				tmp_serve_w[whid] += tmp_factors[cnt] * curr_quad->clientsInQuadrant.size() * m_travelTime.get(quid,whid);
				cnt++;
			}
		}
		
		// build auxiliary model to minimize balancing factor
		// auxiliary balancing constraints
		for (int v=0; v < m_params->nbWarehouses; v++) {
			for (int w=0; w < m_params->nbWarehouses; w++) {
				if ( v == w ) { continue; }
				m_balance[v][w] = m_model->addConstr(m_serve_w[v] <= tmp_serve_w[w], "balance_"+std::to_string(v)+"_"+std::to_string(w));
			}
		}
		// extra constraints
		//  - bind "gamma_iw" and "x_iw"
		//        gamma_iw <= x_iw * UB
		//  - make "gamma" the maximum of "gamma_iw" and minimize it
		//        gamma_iw <= gamma
		GRBConstr* tmp_bind_x = new GRBConstr[nn];
		GRBConstr* tmp_max_gamma = new GRBConstr[nn];
		GRBVar gamma = m_model->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "gamma_max");
		m_model->setObjective(1.0*gamma,GRB_MINIMIZE);
		for (int i=0; i < nn; i++) {
			int whid = m_rev_id[i].first;
			int quid = m_rev_id[i].second;
			tmp_bind_x[i] = m_model->addConstr(tmp_factors[i] <= ub_factor*m_assign[i], "bind_"+std::to_string(whid)+"_"+std::to_string(quid));
			tmp_max_gamma[i] = m_model->addConstr(tmp_factors[i] <= gamma, "bind_"+std::to_string(whid)+"_"+std::to_string(quid));
		}		
		
		// commit changes, solve auxiliary model and retrieve maximum balance factor
		m_model->update();
		m_model->optimize();
		double gamma_max = m_model->get(GRB_DoubleAttr_ObjVal);
		
		std::cout << "relative_balance = " << gamma_max << std::endl;
		
		// reset optimization to original objective with the value of the computed gamma
		//  - remove extra constraints and auxiliary variable
		//  - change balancing constraints back to original variables "x_iw" and use computed factor "gamma"
		//  - revert back to original objective function
		for (int i=0; i < nn; i++) {
			m_model->remove(tmp_bind_x[i]);
			m_model->remove(tmp_max_gamma[i]);
			m_model->remove(tmp_factors[i]);
		}			
		for (int v=0; v < m_params->nbWarehouses; v++) {
			for (int w=0; w < m_params->nbWarehouses; w++) {
				if ( v == w ) { continue; }
				m_model->remove(m_balance[v][w]);
				m_balance[v][w] = m_model->addConstr(m_serve_w[v] <= gamma_max*(1.1)*m_serve_w[w], "balance_"+std::to_string(v)+"_"+std::to_string(w));
			}
		}
		m_model->setObjective(orig_objective,GRB_MINIMIZE);
		m_model->remove(gamma);
		
		// clean up auxialiary stuff
		delete[] tmp_max_gamma;
		delete[] tmp_bind_x;
		delete[] tmp_serve_w;
		delete[] tmp_factors;
	}
	
	if ( m_balance_type == 2 ) {
		
		// optimize minimal balance value
		GRBVar gamma = m_model->addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "gamma_max");
		m_model->setObjective(1.0*gamma,GRB_MINIMIZE);
		for (int v=0; v < m_params->nbWarehouses; v++) {
			for (int w=0; w < m_params->nbWarehouses; w++) {
				if ( v == w ) { continue; }
				m_balance[v][w] = m_model->addConstr(m_serve_w[v] - m_serve_w[w] <= gamma, "balance_"+std::to_string(v)+"_"+std::to_string(w));
			}
		}
		m_model->update();
		m_model->optimize();
		
		// reset optimization to original objective with the value of the computed gamma
		double gamma_max = m_model->get(GRB_DoubleAttr_ObjVal);
		std::cout << "absolute_balance = " << gamma_max << std::endl;
		for (int v=0; v < m_params->nbWarehouses; v++) {
			for (int w=0; w < m_params->nbWarehouses; w++) {
				if ( v == w ) { continue; }
				m_model->remove(m_balance[v][w]);
				m_balance[v][w] = m_model->addConstr(m_serve_w[v] - m_serve_w[w] <= 1.1*gamma_max, "balance_"+std::to_string(v)+"_"+std::to_string(w));
			}
		}
		m_model->setObjective(orig_objective,GRB_MINIMIZE);
		m_model->remove(gamma);
	}
	
	
	m_model->update();
 	//m_model->write("assignment.lp");
 	
  // print elapsed time
  bm_time += cputime();
	std::cout << "bm_time= "; printtime(bm_time); std::cout << std::endl;
	
	return;
};


int WQAssign::solve_check_writesol (bool write_iter, bool quad_level) {
#ifdef DEBUG
	std::cout << "call solve_check_write()" << std::endl;
#endif
	
	// timer
	long int scw_time = - cputime();
	
	// solve current model
	m_model->optimize();
	int ret = -1;
	
	// process solution data
	Quadrant* curr_quad;
	Warehouse* curr_wh;
	
	// reset prior assignment
	for (size_t i=0; i < m_used_quads.size(); i++) {
		m_used_quads[i]->assignedToWarehouse = NULL;
	}
	for (int v=0; v < m_params->nbWarehouses; v++) {
		m_params->paramWarehouses[v].assignedQuadrants.resize(0);
	}
	
	// pass solution to m_params
	if ( m_model->get(GRB_IntAttr_SolCount) > 0.5 ) {
		ret = 0;
		double objval = m_model->get(GRB_DoubleAttr_ObjVal);
		
		if ( m_o > 0 ) {
			std::cout << "model_cost = " << objval << std::endl;
		}
		for (int j = 0; j < m_num_assign; j++) {
			if ( m_assign[j].get(GRB_DoubleAttr_X) > 0.5 ) {
				curr_quad = &m_params->paramQuadrants[m_rev_id[j].second];
				curr_wh = &m_params->paramWarehouses[m_rev_id[j].first];
				curr_quad->assignedToWarehouse = curr_wh;
				curr_wh->assignedQuadrants.push_back(curr_quad);
			}
		}
	}
  else {
    std::cout << "No solution found!" << std::endl;
    exit(1);
  }
	
  // save quadrant assignment to textfile for visualization
	//const char* ofile = 0;
	std::string of_content;
	if (quad_level) of_content = "./data/heatmapData/wh_quad_assignment.txt";
	else of_content = "/data/heatmapData/wh_cust_assignment.txt";
	const char* ofile = of_content.data();
	
	FILE* f = fopen(ofile,"w");
	if ( ofile ) {
		if ( !quad_level ) {
			for (Client c : m_params->paramClients) {		    
				double lat = c.location.lat;
				double lon = c.location.lon;
				double meanWaitingTime = c.averageWaitingTime;
				fprintf(f,"%f %f %.2f %i\n",lat, lon, meanWaitingTime, iter_cnt_write-1);
			}
		}
		else {
			for (int i=0; i < m_params->nbQuadrants; i++) {
				curr_quad = &m_params->paramQuadrants[i];
				int quid = curr_quad->quadrantID;
				double lasw = curr_quad->latSW;
				double losw = curr_quad->lonSW;
				double lane = curr_quad->latNE;
				double lone = curr_quad->lonNE;
				int whid = -1;
				if ( !curr_quad->is_served ) whid = -2; // empty quads are "true" by default, hence stay at "-1"
				if ( curr_quad->assignedToWarehouse ) {
					whid = curr_quad->assignedToWarehouse->wareID;
				}
				
				fprintf(f,"%d %f %f %f %f %d %f %f\n",quid,lasw,losw,lane,lone,whid, curr_quad->meanWaitingTime, curr_quad->maxWaitingTime);
			}
		}
		std::cout << "writing cust -> wh assignments to \"" << ofile << "\"" << std::endl; 
		fclose(f);
	}
	else {
		std::cout << "ERROR: Can't open \"" << ofile << "\"" << std::endl;
		exit(-1);
	}
	
  scw_time += cputime();
	std::cout << "scw_time= "; printtime(scw_time); std::cout << std::endl;
	
  return ret;
};

//===============================================================================


int WQAssign::init() {
#ifdef DEBUG
	std::cout << "call init()" << std::endl;
#endif
	
	initialize_quadrants();						// get the party started
	build_model();										// build initial model
	int ret = solve_check_writesol(); // solve it
	
	return ret;
};


//===============================================================================


Quadrant* WQAssign::calc_quadrant_to_remove(bool use_mean, int max_nb) {
#ifdef DEBUG
	std::cout << "call calc_quadrant_to_remove()" << std::endl;
#endif
	
	Client* curr_client;
	Quadrant* curr_quad;
	
	// calculate mean and max waiting time per quadrant
	//  (using simulation data "averageWaitingTime" and "visitedCount")
	std::vector<double> wait_avg (m_used_quads.size(),-1);
	std::vector<double> wait_max (m_used_quads.size(),-1);
	std::vector<double> num_order_quad (m_used_quads.size(),0);
	
	for (size_t j=0; j < m_used_quads.size(); j++) {
		curr_quad = m_used_quads[j];
		double mean_wait = 0;
		double max_wait = 0;
		bool quad_sim = false; // indicate whether a client in quadrant was simulated
		
		// check all clients in quadrant
		for (size_t i=0; i<curr_quad->clientsInQuadrant.size(); i++) {
			curr_client = curr_quad->clientsInQuadrant[i];
			// skip clients without orders
			if ( curr_client->visitedCount > m_eps ) {
				quad_sim = true;
				num_order_quad[j] += curr_client->visitedCount;
				double awt = curr_client->averageWaitingTime;
				mean_wait += curr_client->visitedCount*awt;
				if ( awt - max_wait > m_eps ) {
					max_wait = awt;
				}
			}
		}
		if ( quad_sim ) {
			mean_wait /= num_order_quad[j];
			wait_avg[j] = mean_wait;
			wait_max[j] = max_wait;
		}
		curr_quad->meanWaitingTime = mean_wait;
		curr_quad->maxWaitingTime = max_wait;
	}
	
	
	// find quadrants with largest avg and largest max waiting time
	Quadrant* max_avg_quad = NULL;
	Quadrant* max_max_quad = NULL;
	double max_wait_avg = -1;
	double max_wait_max = -1;
	bool got_one = false;
	for (int nb=0; nb < max_nb+1; nb++) {
#ifdef DEBUG
		std::cout << "checking quads with " << nb << " active neighbours" << std::endl;
#endif
		for (size_t j=0; j < m_used_quads.size(); j++) {
			
			curr_quad = m_used_quads[j];
			if ( !curr_quad->is_served || !num_order_quad[j] ) { continue; } // skip
			//std::cout << "bef acc" << j << " (" << m_used_quads.size() << ")"<< std::endl;
			// removes isolated quadrants
			if ( curr_quad->num_active_neighbours == nb ) {
				got_one = true;
				double curr_avg_wait = wait_avg[j];
				double curr_max_wait = wait_max[j];
				if ( curr_avg_wait - max_wait_avg > m_eps ) {
					max_wait_avg = curr_avg_wait;
					max_avg_quad = curr_quad;
				}
				if ( curr_max_wait - max_wait_max > m_eps ) {
					max_wait_max = curr_max_wait;
					max_max_quad = curr_quad;
				}
			}
			//std::cout << "aft acc" << j << " (" << m_used_quads.size() << ")" << std::endl;
		}
		if ( got_one ) { break; }
	}
	// return NULL if no quadrant is found
	if ( !got_one ) {
		return 0;
	}
	
	// output and return pointer to quadrant
	if ( use_mean ) {
		if ( m_o > 0 ) {
			std::cout << "quad_to_remove: ID= " << max_avg_quad->quadrantID << ", wait_avg= " << max_wait_avg << " (clientsInQuadrant= " << m_params->paramQuadrants[max_avg_quad->quadrantID].clientsInQuadrant.size() << ")" << std::endl;
			for (size_t i=0; i<max_avg_quad->clientsInQuadrant.size(); i++ ) {
				if ( max_avg_quad->clientsInQuadrant[i]->visitedCount ) {
					std::cout << "client: ID= " << max_avg_quad->clientsInQuadrant[i]->clientID << ", wait_avg= " << max_avg_quad->clientsInQuadrant[i]->averageWaitingTime << std::endl;
				}
			}
		}
		return max_avg_quad;
	}
	else {
		if ( m_o > 0 ) {
			std::cout << "quad_to_remove: ID= " << max_max_quad->quadrantID << ", wait_max= " << max_wait_max << " (clientsInQuadrant= " << m_params->paramQuadrants[max_max_quad->quadrantID].clientsInQuadrant.size() << ")" << std::endl;
			for (size_t i=0; i<max_max_quad->clientsInQuadrant.size(); i++ ) {
				if ( max_max_quad->clientsInQuadrant[i]->visitedCount ) {
					std::cout << "client: ID= " << max_max_quad->clientsInQuadrant[i]->clientID << ", wait_avg= " << max_max_quad->clientsInQuadrant[i]->averageWaitingTime << std::endl;
				}
			}
		}
		return max_max_quad;
	}
};


//===============================================================================


int WQAssign::reoptimize(int reps) {
#ifdef DEBUG
	std::cout << "call reoptimize()" << std::endl;
	
	double total_avg_wait = 0;
	double total_max_wait = 0;
	int cnt_orders = 0;
	
	for (size_t j=0; j < m_used_quads.size(); j++) {
		Quadrant* curr_quad = m_used_quads[j];
		if ( !curr_quad->is_served ) continue; // empty quads aren't looped over, hence no need to check "contains_client"
		for (size_t i=0; i<curr_quad->clientsInQuadrant.size(); i++) {
			Client* curr_client = curr_quad->clientsInQuadrant[i];
			int num_visited = curr_client->visitedCount;
			if ( num_visited < 0.5 ) continue; // skip non-simulated clients
			cnt_orders += num_visited;
			total_avg_wait += num_visited * curr_client->averageWaitingTime;
			if ( total_max_wait < curr_client->averageWaitingTime + m_eps ) {
				total_max_wait = curr_client->averageWaitingTime;
			}
		}
	}
	total_avg_wait /= cnt_orders;
	std::cout << "### orders= " << cnt_orders << " , avg= " << total_avg_wait << " , max= " << total_max_wait << std::endl;
	
	// ==================================
	std::string ostr = "./data/waitProgData/waitprog_"+itos(m_traveltime_cap)+"_"+itos(m_balance_type);
	if ( reps > m_eps ) {
		ostr += "_"+itos(reps)+".txt";
	}
	const char* ofile = ostr.data();
	FILE* f = fopen(ofile,"a");
	if ( ofile ) {				
		fprintf(f,"%.2f %.2f\n",total_avg_wait,total_max_wait);
	}
	else {
		std::cout << "ERROR: Can't open \"" << ofile << "\"" << std::endl;
		exit(-1);
	}
	fclose(f);
	// ==================================
	
#endif
	
	// get quadrant to remove
	Quadrant* quad_to_remove = calc_quadrant_to_remove(false);
	
	if ( !quad_to_remove ) {
		std::cout << "End reoptimize: No suitable quadrant found to remove" << std::endl;
		return 1;
	}
		
	// mark removed quadrant and clients
	Client* client_to_remove;
	int quid = quad_to_remove->quadrantID;
	quad_to_remove->is_served = false;
	//quad_to_remove->assignedToWarehouse = NULL;
	for ( size_t i = 0; i<quad_to_remove->clientsInQuadrant.size(); i++ ) {
		client_to_remove = quad_to_remove->clientsInQuadrant[i];
		client_to_remove->reached = false;
		//quad_to_remove->clientsInQuadrant.resize(0);
	}
	// update neughbours
	for (size_t nb = 0; nb < quad_to_remove->neighbours.size(); nb++ ) {
		Quadrant* nb_quad = quad_to_remove->neighbours[nb];
		if ( nb_quad && nb_quad->num_active_neighbours > -1 ) { nb_quad->num_active_neighbours--; }
	}
	
	// update gurobi model and commit changes
	for (size_t i = 0; i<m_whvar_id[quid].size(); i++) {
		int var_index = m_whvar_id[quid][i];
		m_assign[var_index].set(GRB_DoubleAttr_UB, 0.0);
	}
	m_model->remove(m_served[quid]);
	m_model->update();
	
	// resove the modek
	int ret = solve_check_writesol();
	
	return ret;
};
