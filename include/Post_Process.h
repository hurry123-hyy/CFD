//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#pragma once
#include <vector>

#include "Global.h"
using namespace std;

class Residual {
   public:
    Residual();
    ~Residual() {
    }

   protected:
    VDouble res_L1;
    VDouble res_L2;
    VDouble res_Loo;
    VInt2D max_index;
    int ist, ied, jst, jed;

   public:
    void Compute_Residual();

   protected:
    void OutputResidual();
    bool Stop_by_Residual();
};