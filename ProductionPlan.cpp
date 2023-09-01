#include "ProductionPlan.h"

ProductionPlan::ProductionPlan(int periods, int blocks, int arcs, int units, int facilities, int stockpiles)
{
	system("copy 05_debugger.txt predebugger.txt");
	remove("debugger.txt");

	debugFile.open("05_debugger.txt", ios::app);
	debugFile.setf(ios::fixed, ios::floatfield);
	debugFile.precision(2);

	nPeriods = periods;
	nBlocks = blocks;
	nArcs = arcs;
	uIndicator = units;
	nFacilities = facilities;
	nStockpiles = stockpiles;

	metalPrice = refineCost = mining_cost = discount_rate = mining_cap_low = mining_cap_up = 0.0 ;
	processing_cost = nullptr;
    stockpile_capacity = nullptr;
    processing_cap_low = nullptr;
    processing_cap_up = nullptr;
    processRecovery = nullptr;
    stockpile_recovery = nullptr;

	mdouble.newMTX(bX, nBlocks);
	mdouble.newMTX(bY, nBlocks);
	mdouble.newMTX(bZ, nBlocks);
	mdouble.newMTX(bTons, nBlocks);
	mdouble.newMTX(bGrade, nBlocks);
	mdouble.newMTX(blockValue, nBlocks);
	mdouble.newMTX(wasteValue, nBlocks);
	mdouble.newMTX(processing_cost, nFacilities); // Initialize processing_cost with the correct dimensions
	mdouble.newMTX(processRecovery, nBlocks, nFacilities); // Initialize processRecovery with the correct dimensions
	mdouble.newMTX(stockpile_capacity, nStockpiles); // Initialize stockpile_capacity with the correct dimensions
	mdouble.newMTX(stockpile_recovery, nStockpiles, nFacilities); // Initialize stockpile_recovery with the correct dimensions
	mdouble.newMTX(processing_cap_low, nFacilities); // Initialize processing_cap_low with the correct dimensions
	mdouble.newMTX(processing_cap_up, nFacilities); // Initialize processing_cap_up with the correct dimensions
	

	mint.newMTX(blockID, nBlocks);
	mint.newMTX(fromBlock, nArcs);
	mint.newMTX(toBlock, nArcs);

}

ProductionPlan::~ProductionPlan(void)
{
	mdouble.delMTX(bX, nBlocks);
	mdouble.delMTX(bY, nBlocks);
	mdouble.delMTX(bZ, nBlocks);
	mdouble.delMTX(bTons, nBlocks);
	mdouble.delMTX(bGrade, nBlocks);
	mdouble.delMTX(blockValue, nBlocks);
	mdouble.delMTX(wasteValue, nBlocks);
	mdouble.delMTX(processing_cost, nFacilities); // Delete processingCosts with the correct dimensions
	mdouble.delMTX(processRecovery, nBlocks, nFacilities); // Delete processRecovery with the correct dimensions
	mdouble.delMTX(stockpile_capacity, nStockpiles); // Delete stockpile_capacity with the correct dimensions
	mdouble.delMTX(stockpile_recovery, nStockpiles, nFacilities); // Delete stockpile_recovery with the correct dimensions
	mdouble.delMTX(processing_cap_low, nFacilities); // Delete processing_cap_low with the correct dimensions
	mdouble.delMTX(processing_cap_up, nFacilities); // Delete processing_cap_up with the correct dimensions

	mint.delMTX(blockID, nBlocks);
	mint.delMTX(fromBlock, nArcs);
	mint.delMTX(toBlock, nArcs);
		debugFile.close();
}

void ProductionPlan::ReadData()
{
	debugFile << "Starting to read XYZ and block grades from 01_blockmodel.txt..."
		<< endl << endl << endl;

	std::cout << "Starting to read XYZ and block grades from 01_blockmodel.txt..."
		<< endl << endl << endl;

	ifstream myInput("01_blockmodel.txt", ios::in);
	ifstream myPrecedence("02_precedence.txt", ios::in);

	if (!myInput.is_open())
	{
		std::cerr << "Error: Unable to open file 01_blockmodel.txt for reading" << endl;
		return;
	}
		for (int i = 0; i < nBlocks; i++)
		{
			myInput >> bX[i] >> bY[i] >> bZ[i] >> bGrade[i] >> bTons[i];
			for (int p = 0; p < nFacilities; p++) // Read the processing recovery for each block and each processing facility
            {
                myInput >> processRecovery[i][p]; // where i is the block and p is the corresponding Processing Facility
            }
			//std::cout << (i + 1) << endl;
		}
}

	myInput.close();

	if (!myPrecedence.is_open())
	{
		std::cerr << "Error: Unable to open 02_precedence.txt for reading." << endl;
        return;
    }
		for (int i = 0; i < nArcs; i++)
		{
			myPrecedence >> fromBlock[i] >> toBlock[i];
			//std::cout << (i + 1) << endl;
		}
		myPrecedence.close();

		debugFile << "Finished reading XYZ and block grades from 01_blockmodel.txt."
		<< endl << endl << endl;

		std::cout << "Finished reading XYZ and block grades from 01_blockmodel.txt."
		<< endl << endl << endl;
}

void ProductionPlan::AssignBlockIDS()
{
	debugFile << "Assign IDs ..."
		<< endl << endl << endl;

	std::cout << "Assign IDs ....."
		<< endl << endl << endl;

	int k = 0;

	for (int j = 0; j < nBlocks; j++)
	{
		blockID[k] = k + 1;
		k++;
	}

	ofstream myOutput("03_InitialData.txt", ios::out);
	myOutput.setf(ios::fixed, ios::floatfield);

	    for (int i = 0; i < nBlocks; i++)
    {
        myOutput << setiosflags(ios::left) << setprecision(2)
            << setw(10) << blockID[i]
            << setw(10) << bX[i]
            << setw(10) << bY[i]
            << setw(10) << bZ[i]
            << setw(10) << bTons[i]
            << setw(10) << bGrade[i];
        
        for (int p = 0; p < nFacilities; p++)
        {
            myOutput << setw(10) << processRecovery[i][p];
        }

        myOutput << endl;
    }

    for (int i = 0; i < nArcs; i++)
    {
        myOutput << setiosflags(ios::left)
            << setw(10) << fromBlock[i]
            << setw(10) << toBlock[i] << endl;
    }
	myOutput.close(); 
	debugFile << "Assign IDs ...."
		<< endl << endl << endl;

	std::cout << "Assign IDs ....."
		<< endl << endl << endl;
}

void ProductionPlan::AssignData(double mP, double rC, double mC, double hC, double dR, double mLC, double mUC, double *pLC, double *pUC, double* pC, double* sC, double* sGL, double* sGU, double** sR, double** pR)
	metalPrice = mP;
	refineCost = rC;
	mining_cost = mC;
	handling_cost = hC;
	discount_rate = dR;
	mining_cap_low = mLC;
	mining_cap_up = mUC;

		// Assign stockpile capacities and recoveries
	for (int s = 0; s < nStockpiles; s++)
	{
		stockpile_capacity[s] = sC[s];
		stockpile_grade_low[s] = sGL[s];
		stockpile_grade_up[s] = sGU[s];

		for (int p = 0; p < nFacilities; p++)
		{
			stockpile_recovery[s][p] = sR[s][p];
		}
	}

	// Assign processing costs, capacities, and recoveries
	for (int p = 0; p < nFacilities; p++)
	{
		processing_cost[p] = pC[p];
		processing_cap_low[p] = pLC[p];
		processing_cap_up[p] = pUC[p];

		for (int i = 0; i < nBlocks; i++)
		{
			processRecovery[i][p] = pR[i][p];
		}
	}


void ProductionPlan::CalculateBlockValues()
{
    for (int i = 0; i < nBlocks; i++)
    {
        double maxValue = -std::numeric_limits<double>::max(); // Initialize with a very low value
        double wasteValue = (-1) * mining_cost * bTons[i]; // Value if treated as waste

        // Calculate value for each processing facility
        for (int p = 0; p < nFacilities; p++)
        {
            double value = 0.0;
            if (uIndicator == 1)
            {
                value = ((metalPrice - refineCost) * (bGrade[i] / 100) * processRecovery[i][p] * bTons[i]) - (mining_cost * bTons[i]) - (processing_cost[p] * bTons[i]);
            }
            else if (uIndicator == 2)
            {
                value = ((metalPrice - refineCost) * bGrade[i] * processRecovery[i][p] * bTons[i]) - (mining_cost * bTons[i]) - (processing_cost[p] * bTons[i]);
            }
            if (value > maxValue)
            {
                maxValue = value;
            }
        }

        // Calculate value for each stockpile
        for (int s = 0; s < nStockpiles; s++)
        {
            double stockValue = -mining_cost * bTons[i]; // Value if sent to stockpile
            if (stockValue > maxValue)
            {
                maxValue = stockValue;
            }

            // Calculate value of material from stockpile to each processing facility
            for (int p = 0; p < nFacilities; p++)
            {
                double stockpile_ProcessValue = ((metalPrice - refineCost) * stockpile_grade_low[s] * stockpile_recovery[s][p]) - handling_cost - processing_cost[p];
                if (stockpile_ProcessValue > maxValue)
                {
                    maxValue = stockpile_ProcessValue;
                }
            }
        }

        // Compare with wasteValue and store the maximum value
        if (maxValue <= 0.0)
        {
            blockValue[i] = wasteValue;
        }
        else
        {
            blockValue[i] = maxValue;
        }

        debugFile << setiosflags(ios::left) << setw(10) << (i + 1) << blockValue[i] << endl;
    }

    debugFile << endl;
}


void ProductionPlan::CreateModel()
{
    IloEnv env;
    IloModel mod(env);
    IloTimer solutiontime(env);

    typedef IloArray<IloArray<IloNumVarArray> > Array3D;

    /* Decision variable X - Block to Processing Facility 'p'
	X𝑏p𝑡: This is a binary variable that indicates if block b is mined in period t and moves to processing facility p. 
	It's a 3D array indexed by block, facility, and period. */
    Array3D x(env, nBlocks);
    for(int i = 0; i < nBlocks; i++) 
    {
        x[i] = IloArray<IloNumVarArray>(env, nFacilities);
        for (int p = 0; p < nFacilities; p++)
        {
            x[i][p] = IloNumVarArray(env, nPeriods, 0, 1, ILOINT);
        }
    }

    /* Decision variable K - Block to Stockpile bin 's'
	K𝑏s𝑡: This is a binary variable that indicates if block b is mined in period t and stored in stockpile bin s. 
	It's a 3D array indexed by block, stockpile, and period.*/
    Array3D K(env, nBlocks);
    for(int i = 0; i < nBlocks; i++) 
    {
        K[i] = IloArray<IloNumVarArray>(env, nStockpiles);
        for (int s = 0; s < nStockpiles; s++)
        {
            K[i][s] = IloNumVarArray(env, nPeriods, 0, 1, ILOINT);
        }
    }

    typedef IloArray<IloNumVarArray> Array2D;

    /* Decision variable W - Block to Waste
	W𝑏𝑡: This is a binary variable that indicates if block b is mined in period t and moves to waste. 
	It's a 2D array indexed by block and period. */
    Array2D W(env, nBlocks);
    for(int i = 0; i < nBlocks; i++) 
    {
        W[i] = IloNumVarArray(env, nPeriods, 0, 1, ILOINT);
    }

    /* Continuous variable Q - Material stored in Stockpile bin 's'
	qs𝑡: This is a continuous variable representing the quantity of material stored in stockpile bin s at the end of period t. 
	It's a 2D array indexed by stockpile and period. */
    Array2D q(env, nStockpiles);
    for(int s = 0; s < nStockpiles; s++) 
    {
        q[s] = IloNumVarArray(env, nPeriods, 0, IloInfinity, ILOFLOAT);
    }

    /* Continuous variable Y - Material retrieved from Stockpile bin 's' to Processing Facility 'p'
	qs𝑡: This is a continuous variable representing the quantity of material stored in stockpile bin s at the end of period t. 
	It's a 2D array indexed by stockpile and period. */
    Array3D Y(env, nStockpiles);
    for(int s = 0; s < nStockpiles; s++) 
    {
        Y[s] = IloArray<IloNumVarArray>(env, nFacilities);
        for (int p = 0; p < nFacilities; p++)
        {
            Y[s][p] = IloNumVarArray(env, nPeriods, 0, IloInfinity, ILOFLOAT);
        }
    }

    // Debug output for X variable
    for (int i = 0; i < nBlocks; i++)
    {
        for (int p = 0; p < nFacilities; p++)
        {
            for (int t = 0; t < nPeriods; t++) 
            {
                debugFile << i << "  " << p << " " << t << "  " << x[i][p][t] << endl; 
            }
        }
    }
	// Debug output for K variable
	for (int i = 0; i < nBlocks; i++)
	{
		for (int s = 0; s < nStockpiles; s++)
		{
			for (int t = 0; t < nPeriods; t++) 
			{
				debugFile << "K: " << i << "  " << s << " " << t << "  " << K[i][s][t] << endl; 
			}
		}
	}

	// Debug output for W variable
	for (int i = 0; i < nBlocks; i++)
	{
		for (int t = 0; t < nPeriods; t++) 
		{
			debugFile << "W: " << i << " " << t << "  " << W[i][t] << endl; 
		}
	}

	// Debug output for q variable
	for (int s = 0; s < nStockpiles; s++)
	{
		for (int t = 0; t < nPeriods; t++) 
			{
				debugFile << "q: " << s << "  " << q[s][t] << endl; 
			}
		
	}

	// Debug output for Y variable
	for (int s = 0; s < nStockpiles; s++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			for (int t = 0; t < nPeriods; t++) 
			{
				debugFile << "Y: " << s << "  " << p << " " << t << "  " << Y[s][p][t] << endl; 
			}
		}
	}
    debugFile << endl << endl;
}


IloExpr objExpression(env); //Objective Function Expression

int n = 0; //Loop integer, Loops until the objective completed.

// Objective Function - Discounted value of the production plan (maximize) 
for (int i = 0; i < nBlocks; i++) {
    for (int p = 0; p < nFacilities; p++) {
        for (int t = 0; t < nPeriods; t++) {
            n++;
            objExpression += (blockValue[i] * x[i][p][t]) / (pow((1 + (discount_rate / 100)), (t + 1)));
            debugFile << n << "  " << i << "  " << p << "  " << t << "  " << x[i][p][t] << endl;
        }
    }		
}

for (int i = 0; i < nBlocks; i++) {
    for (int t = 0; t < nPeriods; t++) {
        n++;
        objExpression += (wasteValue * W[i][t]) / (pow((1 + (discount_rate / 100)), (t + 1)));
        debugFile << n << "  " << i << "  WASTE  " << t << "  " << W[i][t] << endl;
    }
}

// Add the stockValue and stock_ProcessValue to the objective function
for (int i = 0; i < nBlocks; i++) {
    for (int s = 0; s < nStockpiles; s++) {
        for (int t = 0; t < nPeriods; t++) {
            n++;
            objExpression += (stockValue * K[i][s][t]) / (pow((1 + (discount_rate / 100)), (t + 1)));
            debugFile << n << "  " << i << "  STOCKPILED  " << s << "  " << t << "  " << K[i][s][t] << endl;
        }
    }
}

for (int s = 0; s < nStockpiles; s++) {
    for (int p = 0; p < nFacilities; p++) {
        for (int t = 0; t < nPeriods; t++) {
            n++;
            objExpression += (stock_ProcessValue * Y[s][p][t]) / (pow((1 + (discount_rate / 100)), (t + 1)));
            debugFile << n << "  " << s << "  STOCKPILE TO PROCESS  " << p << "  " << t << "  " << Y[s][p][t] << endl;
        }
    }
}

debugFile << endl << endl;

IloObjective objFunction(env, objExpression, IloObjective::Maximize);


mod.add(objFunction);

/*The given formulation ensures that each block is either sent to a processing facility, 
sent to a stockpile, or treated as waste only once. The constraint is applied for each block.*/	
for (int i = 0; i < nBlocks; i++) 
{
    IloExpr reserveExpression(env);

    for (int t = 0; t < nPeriods; t++) 
    {
        // Sum over all processing facilities
        for (int p = 0; p < nFacilities; p++) 
        {
            reserveExpression += x[i][p][t];
        }

        // Sum over all stockpiles
        for (int s = 0; s < nStockpiles; s++)
        {
            reserveExpression += K[i][s][t];
        }

        // Add the waste term
        reserveExpression += W[i][t];
    }

    mod.add(reserveExpression <= 1); // Reserve Constraint - Each Block can be mined at most once!
}

	 //Processing Facility Selection Constraint - Each block can assign only one processing facility
	for (int i = 0; i < nBlocks; i++)
	{
		for (int t = 0; t < nPeriods; t++)
		{
			IloExpr facilitySelection(env);

			for (int p = 0; p < nFacilities; p++)
			{
				facilitySelection += x[i][p][t];
			}

			mod.add(facilitySelection <= 1);
		}
	}
	// Stockpile Selection Constraint - If a block is assigned to a stockpile, it should be assigned to only one stockpile bin in that period.
	for (int i = 0; i < nBlocks; i++)
	{
		for (int t = 0; t < nPeriods; t++)
		{
			IloExpr stockpileSelection(env);

			for (int s = 0; s < nStockpiles; s++)
			{
				stockpileSelection += K[i][s][t];
			}

			mod.add(stockpileSelection <= 1);
		}
	}

	

	for (int t = 0; t < nPeriods; t++)
{
    IloExpr miningExpression(env);

    for (int i = 0; i < nBlocks; i++)
    {
        for (int p = 0; p < nFacilities; p++) 
        {
            miningExpression += (bTons[i] * x[i][p][t]); 
        }
        	miningExpression += (bTons[i] * W[i][t]); // Adding the tons treated as waste
		
		for (int s = 0; s < nStockpiles; s++)
		{
			miningExpression += (bTons[i] * K[i][s][t]); // Adding the tons sent to stockpile
		} 
	}

    mod.add(miningExpression >= mining_cap_low); // Lower Mining Capacity Constraint
    mod.add(miningExpression <= mining_cap_up); // Upper Mining Capacity Constraint
}

	
	for (int t = 0; t < nPeriods; t++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			IloExpr processExpression(env);

        for (int i = 0; i < nBlocks; i++)
        {
            if (blockValue[i] > 0)
            {
                processExpression += (bTons[i] * x[i][p][t]);
            }
        }
		for (int s = 0; s < nStockpiles; s++)
		{
			processExpression += (Y[s][p][t]); // Adding the tons sent to processing facility from stockpile bin s.
		}

        mod.add(processExpression >= processing_cap_low[p]);
		mod.add(processExpression <= processing_cap_up[p]);
		}
	}
	/*Stockpile Creation Constraints */
	// Stockpile creation constraint for t=1
	for (int s = 0; s < nStockpiles; s++)
	{
		IloExpr stockpileExpression(env);

		for (int i = 0; i < nBlocks; i++)
		{
			stockpileExpression += (bTons[i] * K[i][s][0]); // Material added to stockpile from blocks
		}

		for (int p = 0; p < nFacilities; p++)
		{
			stockpileExpression -= Y[s][p][0]; // Material removed from stockpile to processing facilities
		}

		stockpileExpression -= q[s][0]; // Material remaining in stockpile at the end of period 1

		mod.add(stockpileExpression == 0);
	}

	// Stockpile creation constraint for t>=2
	for (int s = 0; s < nStockpiles; s++)
	{
		for (int t = 1; t < nPeriods; t++)
		{
			IloExpr stockpileExpression(env);

			for (int i = 0; i < nBlocks; i++)
			{
				stockpileExpression += (bTons[i] * K[i][s][t]); // Material added to stockpile from blocks
			}

			for (int p = 0; p < nFacilities; p++)
			{
				stockpileExpression -= Y[s][p][t]; // Material removed from stockpile to processing facilities
			}

			stockpileExpression -= q[s][t]; // Material remaining in stockpile at the end of period t
			stockpileExpression += q[s][t-1]; // Material remaining from previous period (t-1)

			mod.add(stockpileExpression == 0);
		}
	}
	// Stockpile Capacity Constraints
	for (int s = 0; s < nStockpiles; s++)
	{
		for (int t = 0; t < nPeriods; t++)
		{
			// Ensure quantity in stockpile is non-negative
			mod.add(q[s][t] >= 0);

			// Ensure quantity in stockpile does not exceed its capacity
			mod.add(q[s][t] <= stockpile_capacity[s]);
		}
	}
		// Stockpile Grade Constraints
	for (int s = 0; s < nStockpiles; s++)
	{
		// For the first time period t=1
		IloExpr numeratorGrade1(env);
		IloExpr denominatorGrade1(env);
		for (int i = 0; i < nBlocks; i++)
		{
			numeratorGrade1 += (bTons[i] * bGrade[i] * K[i][s][0]);
			denominatorGrade1 += (bTons[i] * K[i][s][0]);
		}
		for (int p = 0; p < nFacilities; p++)
		{
			denominatorGrade1 -= Y[s][p][0];
		}
		mod.add(stockpile_grade_low[s] <= numeratorGrade1 / denominatorGrade1);
		mod.add(numeratorGrade1 / denominatorGrade1 <= stockpile_grade_up[s]);

		// For the time periods t>=2
		for (int t = 1; t < nPeriods; t++)
		{
			IloExpr numeratorGradeT(env);
			IloExpr denominatorGradeT(env);
			for (int i = 0; i < nBlocks; i++)
			{
				numeratorGradeT += (bTons[i] * bGrade[i] * K[i][s][t]);
				denominatorGradeT += (bTons[i] * K[i][s][t]);
			}
			numeratorGradeT += (q[s][t-1] * stockpile_grade_low[s]); // Assuming stockpile_grade_low stores the grade of stockpile at previous time period
			denominatorGradeT += q[s][t-1];
			for (int p = 0; p < nFacilities; p++)
			{
				denominatorGradeT -= Y[s][p][t];
			}
			mod.add(stockpile_grade_low[s] <= numeratorGradeT / denominatorGradeT);
			mod.add(numeratorGradeT / denominatorGradeT <= stockpile_grade_up[s]);
		}
}
	// Determine the number of precedence blocks for each block
	vector<int> precedenceCount(nBlocks, 0); // Initialize to zero

	for (int j = 0; j < nArcs; j++)
	{
		tBlock = toBlock[j];
		precedenceCount[tBlock - 1]++; // Assuming block IDs start from 1
	}

	// Precedence Constraint	   
	int fBlock = 0, tBlock = 0;

	for (int j = 0; j < nArcs; j++)
	{
		fBlock = fromBlock[j]; // This is the preceding block b'
		tBlock = toBlock[j];   // This is the current block b

		int n = precedenceCount[tBlock - 1]; // Number of preceding blocks for the current block

		debugFile << "fBlock = " << fBlock << "  " << fBlock - 1 << "   tBlock = " << tBlock << "  " << tBlock - 1 << endl;

		for (int t = 0; t < nPeriods; t++)
		{
			IloExpr sumExpression(env);

			// Sum over all preceding blocks and all facilities for the current block
			for (int tPrime = 0; tPrime <= t; tPrime++)
			{
				for (int p = 0; p < nFacilities; p++)
				{
					sumExpression += x[fBlock - 1][p][tPrime];
				}
				for (int s = 0; s < nStockpiles; s++) 
				{
					sumExpression += K[fBlock - 1][s][tPrime]; 
				}
				sumExpression += W[fBlock - 1][tPrime];
			}

			for (int p = 0; p < nFacilities; p++) // Loop over facilities
			{
				// Add the constraint for the block going to processing facility
				mod.add(n * x[tBlock - 1][p][t] - sumExpression <= 0);
			}

			for (int s = 0; s < nStockpiles; s++) // Loop over stockpile bins
			{
				// Add the constraint for the block going to stockpile bin
				mod.add(n * K[tBlock - 1][s][t] - sumExpression <= 0); // Assuming K is the decision variable for stockpile bins
			}

			// Add the constraint for the block going to waste
			mod.add(n * W[tBlock - 1][t] - sumExpression <= 0);
		}
	}

	debugFile << endl;

		IloCplex cplex(env);
	cplex.extract(mod);
	cplex.exportModel("04_Formulation.lp");

	cplex.setParam(IloCplex::ClockType, 1);
	cplex.setParam(IloCplex::EpGap, 0.01);
	cplex.setParam(IloCplex::TiLim, 2592000);

	auto start = std::chrono::system_clock::now();
	solutiontime.start();
	cplex.solve();
	solutiontime.stop();
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	ofstream myOutput("05_OptimalPP.txt", ios::out);
	myOutput.setf(ios::fixed, ios::floatfield);

	if (myOutput.is_open())
	{
		IloNum objValue = cplex.getObjValue();

		myOutput << "Discounted value of the production plan = " << setprecision(0) << objValue << endl;
		myOutput << "Solution time = " << setprecision(5) << solutiontime.getTime() << " seconds" << endl << endl;
		myOutput << "Solution time = " << setprecision(5) << elapsed_seconds.count() << " seconds" << endl << endl;

		// Output for blocks sent to processing
		for (int t = 0; t < nPeriods; t++)
		{
			for (int i = 0; i < nBlocks; i++)
			{
				for (int p = 0; p < nFacilities; p++)
				{
					if (cplex.getValue(x[i][p][t]) > 0)
					{
						myOutput << setiosflags(ios::left) << setprecision(0)
							<< setw(10) << bX[i]
							<< setw(10) << bY[i]
							<< setw(10) << bZ[i]
							<< setw(10) << bTons[i]
							<< setprecision(5)
							<< setw(15) << bGrade[i]
							<< setw(10) << setprecision(0) << (t+1)
							<< setw(10) << p + 1
							<< setw(10) << cplex.getValue(x[i][p][t]) << endl;
					}
				}
			}
		}

		// Output for blocks sent to stockpile
		for (int t = 0; t < nPeriods; t++)
		{
			for (int i = 0; i < nBlocks; i++)
			{
				for (int s = 0; s < nStockpiles; s++)
				{
					if (cplex.getValue(K[i][s][t]) > 0)
					{
						myOutput << setiosflags(ios::left) << setprecision(0)
							<< setw(10) << bX[i]
							<< setw(10) << bY[i]
							<< setw(10) << bZ[i]
							<< setw(10) << bTons[i]
							<< setprecision(5)
							<< setw(15) << bGrade[i]
							<< setw(10) << setprecision(0) << (t+1)
							<< setw(10) << "Stockpile " << s + 1
							<< setw(10) << cplex.getValue(K[i][s][t]) << endl;
					}
				}
			}
		}

		// Output for blocks sent to waste
		for (int t = 0; t < nPeriods; t++)
		{
			for (int i = 0; i < nBlocks; i++)
			{
				if (cplex.getValue(W[i][t]) > 0)
				{
					myOutput << setiosflags(ios::left) << setprecision(0)
						<< setw(10) << bX[i]
						<< setw(10) << bY[i]
						<< setw(10) << bZ[i]
						<< setw(10) << bTons[i]
						<< setprecision(5)
						<< setw(15) << bGrade[i]
						<< setw(10) << setprecision(0) << (t+1)
						<< setw(10) << "Waste"
						<< setw(10) << cplex.getValue(W[i][t]) << endl;
				}
			}
		}
			// Output for material stored in each stockpile bin for each period
		myOutput << "\nMaterial stored in Stockpile bins:\n";
		for (int s = 0; s < nStockpiles; s++)
		{
			for (int t = 0; t < nPeriods; t++)
			{
				myOutput << "Stockpile " << s+1 << ", Period " << t+1 << ": " << cplex.getValue(q[s][t]) << " tons" << endl;
			}
		}

		// Output for material retrieved from stockpile bins and processed
		myOutput << "\nMaterial retrieved from Stockpile bins and processed:\n";
		for (int s = 0; s < nStockpiles; s++)
		{
			for (int p = 0; p < nFacilities; p++)
			{
				for (int t = 0; t < nPeriods; t++)
				{
					if (cplex.getValue(Y[s][p][t]) > 0)
					{
						myOutput << "Stockpile " << s+1 << " to Facility " << p+1 << ", Period " << t+1 << ": " << cplex.getValue(Y[s][p][t]) << " tons" << endl;
					}
				}
			}
		}
		myOutput.close();
	

		ofstream myOutput1("06_PitProcess.txt", ios::out); // Print out the blocks that are in the UPL and goes processing
		myOutput1.setf(ios::fixed, ios::floatfield);

		for (int i = 0; i < nBlocks; i++)
		{
			for (int p = 0; p < nFacilities; p++) // Added this loop for processing facilities
			{
				for (int t = 0; t < nPeriods; t++)
				{
					if (cplex.getValue(x[i][p][t]) > 0)
					{
						myOutput1 << setiosflags(ios::left) << setprecision(0)
							<< setw(10) << blockID[i]
							<< setw(10) << (p)
							<< setw(10) << (t) << endl;
					}
				}
			}			
		}
	
		myOutput1.close();
	}

void ProductionPlan::AnalysisOfResults()
{
	ofstream resultsFile("07_results.txt", ios::out);
	resultsFile.setf(ios::fixed, ios::floatfield);

	resultsFile << setiosflags(ios::left) << setprecision(0);

	ofstream modelFile("08_blockmodelSchedule.txt", ios::out);
	modelFile.setf(ios::fixed, ios::floatfield);

	ofstream blockToFacilityFile("09_blockToFacility.txt", ios::out);
	blockToFacilityFile.setf(ios::fixed, ios::floatfield);

	ofstream processedMaterialFile("10_processedMaterial.txt", ios::out);
	processedMaterialFile.setf(ios::fixed, ios::floatfield);

	double *qtWaste, *qtOre, *qtMetal, *cutG, *avgG, *cashF, *netPV, *totalProcessCost;
	double **totalProcessed, **avgG_Stockpile, **qt_Stockpiled, **g_Stockpiled;

	mdouble.newMTX(qtWaste, nPeriods);
	mdouble.newMTX(qtOre, nPeriods);
	mdouble.newMTX(qtMetal, nPeriods);
	mdouble.newMTX(cutG, nPeriods);
	mdouble.newMTX(avgG, nPeriods);
	mdouble.newMTX(cashF, nPeriods);
	mdouble.newMTX(netPV, nPeriods);
	mdouble.newMTX(totalProcessCost, nPeriods);
	mdouble.newMTX(totalProcessed, nFacilities, nPeriods);
	mdouble.newMTX(avgG_Stockpile, nStockpiles, nPeriods);
	mdouble.newMTX(qt_Stockpiled, nStockpiles, nPeriods);
	mdouble.newMTX(g_Stockpiled, nStockpiles, nPeriods);

	for (int t = 0; t < nPeriods; t++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			totalProcessed[p][t] = 0.0;

			for (int i = 0; i < nBlocks; i++)
			{
				totalProcessed[p][t] += bTons[i] * cplex.getValue(X[i][p][t]);
			}
			for (int s = 0; s < nStockpiles; s++)
			{
				totalProcessed[p][t] += cplex.getValue(Y[s][p][t]);
			}
		}
	}
	for (int t=0; t<nPeriods; t++)
	{
		totalProcessCost[t] = 0.0;

		for (int p=0; p<nFacilities; p++)
		{
			totalProcessCost[t] += totalProcessed[p][t] * processCost[p];
		}
	}

	for (int p = 0; p < nFacilities; p++)
	{
		for (int t = 0; t < nPeriods; t++)
		{
			processedMaterialFile << setiosflags(ios::left) << setprecision(0)
				<< setw(10) << "Processing Facility " << (p + 1) << " processes in period " << (t + 1) << ": "
				<< setw(10) << totalProcessed[p][t] << " tons of material" << endl;
		}
	}

	for (int t = 0; t < nPeriods; t++)
	{
			double sumG = 0.0;
			qtWaste[t] = 0.0;
			qtOre[t] = 0.0;
			qtMetal[t] = 0.0;
			cutG[t] = 1000000000.0;
			avgG[t] = 0.0;	
			netPV[t] = 0.0;
			cashF[t] = 0.0;

	for (int i = 0; i < nBlocks; i++)  
	{
		qtWaste[t] += bTons[i] * cplex.getValue(W[i][t]);
	}

	for (int i = 0; i < nBlocks; i++)  // Loop over all blocks
	{
		for (int p = 0; p < nFacilities; p++)  // Loop over all facilities
		{
			qtOre[t] += bTons[i] * cplex.getValue(X[i][p][t]);
			qtMetal[t] += processRecovery[i][p] * bTons[i] * cplex.getValue(X[i][p][t]);
			sumG += bGrade[i] * bTons[i] * cplex.getValue(X[i][p][t]);
		}
	}

	double cutG_bp = 1000000000.0;  // Initialize with a large value
	double cutG_sp = 1000000000.0;  // Initialize with a large value

	for (int p = 0; p < nFacilities; p++)
	{
		for (int i = 0; i < nBlocks; i++)
		{
			double tempCutG = (metalPrice - refineCost - processCost[p]) / processRecovery[i][p];
			if (tempCutG < cutG[t])
			{
				cutG[t] = tempCutG;
			}
		}	
	}

	for (int s = 0; s < nStockpiles; s++)
	{
		double tempCutG = (metalPrice - refineCost - processCost[s]) / stockpile_recovery[s];
		if (tempCutG < cutG[t])
		{
			cutG[t] = tempCutG;
		}
	}


	for (int s = 0; s < nStockpiles; s++)  
	{	
		qt_Stockpiled[s][t] = 0.0;
		g_Stockpiled[s][t] = 0.0;
		for (int i = 0; i < nBlocks; i++)  
		{
			qt_Stockpiled[s][t] += bTons[i] * cplex.getValue(K[i][s][t]);
			g_Stockpiled[s][t] += bGrade[i] * bTons[i] * cplex.getValue(K[i][s][t]);
		}
	if (t == 0)
	{
		avgG_Stockpile[s][t] = (g_Stockpiled[s][t] * (1 - (cplex.getValue(Y[s][p][t]) / qt_Stockpiled[s][t]))) / cplex.getValue(q[s][t]);
	}
	else
	{
		avgG_Stockpile[s][t] = (g_Stockpiled[s][t] - avgG_Stockpile[s][t-1] * cplex.getValue(Y[s][p][t])) / cplex.getValue(q[s][t]);
	}

	for (int p = 0; p < nFacilities; p++)  
	{
		qtMetal[t] += avgG_Stockpile[s][t] * stockpile_recovery[s][p] * cplex.getValue(Y[s][p][t]);
		sumG += stockpile_grade[s] * cplex.getValue(Y[s][p][t]);
	}
	}
	if (qtOre[t] > 0)
	{
		avgG[t] = (sumG / qtOre[t]);
	}
		else
		{
			avgG[t] = 0.0;
		}
		// Apply the uIndicator logic
	if (uIndicator == 1)
	{
		qtMetal[t] = (qtMetal[t] * (avgG[t] / 100));
	}

	if (uIndicator == 2)
	{
		qtMetal[t] = (qtMetal[t] * avgG[t]);
	}
	}

	for (int t=0; t < nPeriods; t++ )
	{
		cashF[t] = ((metalPrice - refineCost) * qtMetal[t]) - mininig_cost * (qtWaste[t] + qtOre[t]) - totalProcessCost[t];
	}

	for (int t=0; t < nPeriods; t++ )
	{
		for (int t1=t; t1 < nPeriods; t1++ )
		{
			if (t == 0)
			{
				netPV[t] += (cashF[t1] / (pow((1 + (discount_rate / 100)), (t1 + 1))));
			}
			else
			{
				netPV[t] += (cashF[t1] / (pow((1 + (discount_rate / 100)), (t1 - t + 1))));
			}
		}
	}

	resultsFile << endl;
	resultsFile << "Overall Production Schedule Results:" << endl;
	resultsFile << setw(10) << "Period (Year)" 
				<< setw(20) << "Quantity of Waste (tonnes)"
				<< setw(20) << "Cut-off Grade (%)"
				<< setw(20) << "Average Grade (%)"
				<< setw(20) << "Quantity of Ore (tonnes)"
				<< setw(20) << "Quantity of Metal (tonnes)"
				<< setw(20) << "Cash Flow ($)"
				<< setw(20) << "NPV ($)" << endl;

	resultsFile << setw(10) << "" << setw(20) << "" << setw(20) << "" << setw(20) << "" << setw(20) << "" << setw(20) << "" << setw(20) << "" << setw(20) << "" << endl;

	for (int t = 0; t < nPeriods; t++)
	{
		resultsFile << setw(10) << (t+1)
					<< setw(20) << setprecision(0) << qtWaste[t]
					<< setw(20) << setprecision(5) << cutG[t]
					<< setw(20) << setprecision(5) << avgG[t]
					<< setw(20) << setprecision(0) << qtOre[t] 
					<< setw(20) << setprecision(0) << qtMetal[t] 
					<< setw(20) << setprecision(2) << cashF[t]  // Assuming cash flow might have cents
					<< setw(20) << setprecision(2) << netPV[t] << endl;  // Assuming NPV might have cents
	}
	resultsFile << "\nProcessing Facility Information:\n";
	resultsFile << setw(10) << "Period (Year)" 
				<< setw(15) << "Facility"
				<< setw(20) << "Total Processed (tonnes)"
				<< setw(25) << "Total Processing Cost ($)" << endl;

	for (int t = 0; t < nPeriods; t++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			resultsFile << setw(10) << (t+1)
						<< setw(15) << (p+1)
						<< setw(20) << setprecision(0) << totalProcessed[p][t]
						<< setw(25) << setprecision(2) << totalProcessed[p][t] * processCost[p] << endl;
		}
	}

	// Information for each stockpile
	resultsFile << "\nStockpile Information:\n";
	resultsFile << setw(10) << "Period (Year)" 
				<< setw(15) << "Stockpile Bin"
				<< setw(25) << "Average Grade (%)"
				<< setw(25) << "Quantity Stockpiled (tonnes)"
				<< setw(25) << "Grade Stockpiled (%)" << endl;

	for (int t = 0; t < nPeriods; t++)
	{
		for (int s = 0; s < nStockpiles; s++)
		{
			resultsFile << setw(10) << (t+1)
						<< setw(15) << (s+1)
						<< setw(25) << setprecision(5) << avgG_Stockpile[s][t]
						<< setw(25) << setprecision(0) << qt_Stockpiled[s][t]
						<< setw(25) << setprecision(5) << (g_Stockpiled[s][t] / qt_Stockpiled[s][t]) * 100 << endl;  // Assuming grade is a percentage
		}
	}

	mdouble.delMTX(qtWaste, nPeriods);
	mdouble.delMTX(qtOre, nPeriods);
	mdouble.delMTX(qtMetal, nPeriods);
	mdouble.delMTX(cutG, nPeriods);
	mdouble.delMTX(avgG, nPeriods);
	mdouble.delMTX(cashF, nPeriods);
	mdouble.delMTX(netPV, nPeriods);  
	mdouble.delMTX(totalProcessed, nFacilities);
	mdouble.delMTX(totalProcessCost, nPeriods);  
	mdouble.delMTX(avgG_Stockpile, nStockpiles);  
	mdouble.delMTX(qt_Stockpiled, nStockpiles);  
	mdouble.delMTX(g_Stockpiled, nStockpiles);  

}
