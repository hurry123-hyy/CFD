//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#include <iostream>
#include "Spatial_Derivative.h"

#include "Blunt_Solver.h"
#include "Flux_Solver.h"
#include "Geometry.h"
#include "Global.h"

VDouble3D rhs;
void Solve_Spatial_Derivative() {
    Spatial_Derivative* spatial_derivative = new Spatial_Derivative();

    spatial_derivative->Compute_Spatial_Derivative();

    delete spatial_derivative;
}

Spatial_Derivative::Spatial_Derivative() {
    Get_IJK_Region(ist, ied, jst, jed);
    Allocate_3D_Vector(rhs, num_half_point_x, num_half_point_y, num_of_primitive_vars);
}

void Spatial_Derivative::Compute_Spatial_Derivative() {
    if (solve_direction == 'x') {
        this->Spatial_Derivative_X();
    } else if (solve_direction == 'y') {
        this->Spatial_Derivative_Y();
    } else {
        cout << "出错，请检查！" << endl;
    }
}

void Spatial_Derivative::Spatial_Derivative_X() {
    VInt2D& marker = mesh->Get_Marker_Q();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 1; i <= ied + 1; i++)  // i=0和i=ied+2的空间导数项没有值
    {
        for (int j = jst; j <= jed; j++) {
            if (marker[i][j] == 0)
                continue;

            VDouble rhsVector(num_of_primitive_vars);
            for (int iVar = 0; iVar < num_of_primitive_vars; iVar++) {
                rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i - 1][j][iVar]) / dx;
            }
            rhs[i][j] = rhsVector;
        }
    }
}

void Spatial_Derivative::Spatial_Derivative_Y() {
    VInt2D& marker = mesh->Get_Marker_Q();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = ist; i <= ied; i++) {
        for (int j = 1; j <= jed + 1; j++)  // j=0和j=jed+2的空间导数项没有值
        {
            if (marker[i][j] == 0)
                continue;

            VDouble rhsVector(num_of_primitive_vars);
            for (int iVar = 0; iVar < num_of_primitive_vars; iVar++) {
                rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i][j - 1][iVar]) / dy;
            }
            rhs[i][j] = rhsVector;
        }
    }
}
