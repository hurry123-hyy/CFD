﻿#include "2D_Euler_Solver.h"
#include "Global.h"

int main(int argc, char ** argv )
{
	Simulation * two_dim_Euler_Solver = new Simulation();

	two_dim_Euler_Solver->Run();

	delete two_dim_Euler_Solver;

	return 0;
}

void Simulation::Run()
{
	Init_Global_Param();

	Generate_Mesh();

	Init_Flow();

	for (current_step = 0; current_step < max_num_of_steps; ++current_step)
	{
		Set_Solve_Direction('x');
		Load_Q();
		Compute_Boundary();
		Time_Integration();

		Set_Solve_Direction('y');
		Time_Integration();

		Post_Solve();

		if (Need_Stop_Iteration())
		{
			break;
		}
	}

	Output_Flowfield();	
}

void Test()
{
	vector < vector< double > > A = { {1,2,3,0.5},{3,2,1,1.2},{1,3,2,2.1} };
	//vector < vector< double > > B = { {2,3,4},{3,4,2},{4,2,3},{2,4,7} };
	//vector < vector< double > > B = { {2,3},{3,4},{4,2},{2,4} };

	//vector < vector< double > > C;
	//Allocate_2D_Vector(C,3,2);
	//MatrixMultiply(A,B,C,3,4,2);

	vector< double > B = { 2, 3, 9, 12 };
	vector< double > C(3);
	MatrixMultiply(A, B, C, 3, 4);
}
