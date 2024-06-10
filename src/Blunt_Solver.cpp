//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#include "Blunt_Solver.h"

#include <chrono>
#include <ctime>

#include "Global.h"

#ifdef _OPENMPH
#include <omp.h>

#include <cstdlib>
#include <iostream>
#endif

int main(int argc, char** argv) {
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
    InitializeOpenMP(argv);
#endif
    Simulation* blunt_Solver = new Simulation();

    blunt_Solver->Run();

    delete blunt_Solver;

    // 记录结束时间
    auto end = std::chrono::high_resolution_clock::now();

    // 计算持续时间
    std::chrono::duration<double> duration = end - start;

    // 输出运行时间
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    return 0;
}

void Simulation::Run() {
    Init_Global_Param();
    // return;

    Generate_Mesh_Blunt();

    Init_Flow_Blunt_Body();

    for (current_step = 1; current_step <= max_num_of_steps; current_step++) {  // 默认迭代5k
        Solve_Time_Step();

        Compute_Boundary_Blunt_Body();

        Set_Solve_Direction('x');
        Time_Integration();

        Set_Solve_Direction('y');
        Time_Integration();

        Post_Solve();

        if (Need_Stop_Iteration()) {
            break;
        }

    }
        Set_Field();

    Output_Flowfield();
}

#ifdef _OPENMP
void InitializeOpenMP(char** argv) {
    int num_threads = 20;  // 无命令行参数，默认为单线程
    if (argv[1] != NULL) {
        num_threads = atoi(argv[1]);
    }
    omp_set_num_threads(num_threads);

#pragma omp parallel
    {
        // cout << omp_get_num_threads() << endl;
        if (omp_get_thread_num() == 0) {
            cout << "OpenMP并行已开启，线程总数为: " << omp_get_num_threads() << endl;
        }
    }
}
#endif
