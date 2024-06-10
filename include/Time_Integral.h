//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#pragma once
#include "Global.h"

class Time_Marching_Solver {
   public:
    Time_Marching_Solver();
    ~Time_Marching_Solver() {
    }

   public:
    void Time_Marching();

   protected:
};

void Update_Flowfield(int iStage);
void Update_Flowfield_X(int iStage);
void Update_Flowfield_Y(int iStage);
void SolutionFix(VDouble& primitiveVector, int i, int j);