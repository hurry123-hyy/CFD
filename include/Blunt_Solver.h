//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#pragma once
#include <string>
class Simulation {
   public:
    Simulation() {
    }
    ~Simulation() {
    }

   public:
    void Run();
};

void InitializeOpenMP(char** argv);

void Init_Global_Param();
void Init_Flow_Blunt_Body();
void Compute_Boundary_Blunt_Body();
void Load_Q();
void Set_Solve_Direction(char direction);
void Solve_QlQr();
void Solve_Flux();
void Solve_Spatial_Derivative();
void Solve_Time_Step();
void Time_Integration();
void Post_Solve();
void Compute_Residual();
bool Stop_by_Residual();
void Output_Flowfield();
void Input_Parameters();
void Set_Field();
void Generate_Mesh_Blunt();
void Init_Flow_Blunt_Body();
void Compute_Boundary_Blunt_Body();