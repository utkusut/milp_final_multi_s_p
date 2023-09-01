#include "ProductionPlan.h"

//remove "#include <ilcplex/ilocplex.h>" from the fist line and add to new class.

using namespace std; 

int main(){ 
 
	std::cout << "In main" << std::endl;

	ifstream file("00_baseParameters.txt", ios::in);
	ofstream outfile("progTime.txt", ios::out);

	int periods = 0; // Number of periods
	int blocks = 0; // Number of blocks
	int arcs = 0; // Number of precedence arcs
	int units = 0; // Number of mining units
	int facilities = 0; // Number of processing facilities
	double metal_price = 0.0, refining_cost = 0.0, mining_cost = 0.0, discount_rate = 0.0, handling_cost = 0.0; 
	double mining_cap_low = 0.0, mining_cap_up = 0.0;
	int stockpiles = 0; // Number of stockpiles

	file >> periods; // Read the number of periods from the input file
	file >> blocks; // Read the number of blocks from the input file
	file >> arcs; // Read the number of precedence arcs from the input file
	file >> units; // Read the number of mining units from the input file
	file >> facilities; // Read the number of processing facilities from the input file
		double* processing_cost = new double[facilities]; // Create an array for processing costs
		double* processing_cap_low = new double[facilities]; // Create an array for processing capacity lower bounds
		double* processing_cap_up = new double[facilities]; // Create an array for processing capacity upper bounds
	file >> metal_price >> refining_cost >> mining_cost >> handling_cost >> discount_rate; // Read the metal price, refining cost, mining cost, handling cost and discount rate from the input file
	file >> mining_cap_low >> mining_cap_up; // Read the mining capacity lower and upper bounds from the input file
	
    for (int p = 0; p < facilities; p++) 
	{ 
        file >> processing_cost[p]; // Read processing costs for each facility from the input file
		file >> processing_cap_low[p] >> processing_cap_up[p]; // Read processing capacities for each facility from the input file
    }
	file >> stockpiles; // Read the number of stockpiles from the input file
		double* stockpile_capacity = new double[stockpiles]; // Create an array for stockpile capacities
		double** stockpile_recovery = new double*[stockpiles]; // Create a two-dimensional array for stockpile recoveries
		double* stockpile_grade_low = new double[stockpiles]; // Create an array for stockpile grade lower bounds
		double* stockpile_grade_up = new double[stockpiles]; // Create an array for stockpile grade upper bounds
		for (int s = 0; s < stockpiles; s++) 
		{
			file >> stockpile_capacity[s]; // Read stockpile capacities from the input file
			file >> stockpile_grade_low[s];
			file >> stockpile_grade_up[s]; // Read stockpile grade lower and upper bounds from the input file
			stockpile_recovery[s] = new double[facilities]; // Create a two-dimensional array for stockpile recoveries
		for (int p = 0; p < facilities; p++) 
			file >> stockpile_recovery[s][p]; // Read stockpile recoveries for each facility from the input file
		}

	outfile << "Number of Blocks: " << blocks << endl; 
	outfile << "Number of precedence arcs: " << arcs << endl;
	outfile << "Number of processing facilities: " << facilities << endl;
	outfile << "Number of stockpiles: " << stockpiles << endl;
	outfile << "Number of periods: " << periods << endl;

	ProductionPlan pp(periods, blocks, arcs, units, facilities, stockpiles); // Update the constructor with the number of processing facilities and stockpiles

	pp.ReadData();
	pp.AssignBlockIDS();

	pp.AssignData(metal_price, refining_cost, mining_cost, handling_cost, discount_rate, mining_cap_low, mining_cap_up, processing_cap_low, processing_cap_up, processing_cost, stockpile_capacity, stockpile_grade_low, stockpile_grade_up, stockpile_recovery);

	pp.CalculateBlockValues();

	pp.CreateModel();
	pp.AnalysisOfResults();

	// delete[] processing_cost;  Delete the processing cost array
	//delete[] processing_cap_low;  Delete the processing capacity lower bound array
	// delete[] processing_cap_up; Delete the processing capacity upper bound array 


return 0;
}
