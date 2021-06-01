#include "2D_Euler_Solver.h"
#include "Spatial_Derivative.h"
#include "Geometry.h"
#include "Global.h"
#include "Flux_Solver.h"
#include <iostream>

VDouble3D rhs;
void Solve_Spatial_Derivative()
{
	auto* spatial_derivative = new Spatial_Derivative();

	spatial_derivative->Compute_Spatial_Derivative();

	delete spatial_derivative;
}

Spatial_Derivative::Spatial_Derivative()
{
	Get_IJK_Region(ist, ied, jst, jed);
	Allocate_3D_Vector(rhs, num_half_point_x, num_half_point_y, num_of_prim_vars);
}

void Spatial_Derivative::Compute_Spatial_Derivative()
{
	if (method_of_flux == 4)
	{
		if (solve_direction == 'x')
		{
			this->Spatial_Derivative_WCNS_X();
		}
		else if (solve_direction == 'y')
		{
			this->Spatial_Derivative_WCNS_Y();
		}
	}
	else
	{
		if (solve_direction == 'x')
		{
			this->Spatial_Derivative_X();
		}
		else if (solve_direction == 'y')
		{
			this->Spatial_Derivative_Y();
		}
		else
		{
			cout << "�������飡" << endl;
		}
	}

}

void Spatial_Derivative::Spatial_Derivative_X()
{
	VInt2D& marker = mesh->Get_Marker();

	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i - 1][j][iVar]) / dx;
			}
			rhs[i][j] = rhsVector;
		}
	}
}

void Spatial_Derivative::Spatial_Derivative_Y()
{
	VInt2D& marker = mesh->Get_Marker();

	for (int i = ist; i < ied - 1; i++)
	{
		for (int j = jst; j < jed - 1; j++)
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = -(fluxVector[i][j][iVar] - fluxVector[i][j - 1][iVar]) / dy;
			}
			rhs[i][j] = rhsVector;
		}
	}
}

void Spatial_Derivative::Spatial_Derivative_WCNS_X()
{
	double ds = dy;
	double a  = 75.0 / 64.0, b = -25.0 / 384.0, c = 3.0 / 640;

	VInt2D& marker = mesh->Get_Marker();
	for (int j = jst; j < jed - 1; j++)
	{
		for (int i = ist; i < ied - 1; i++)
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = a / ds * (fluxVector[i    ][j][iVar] - fluxVector[i - 1][j][iVar])
								+ b / ds * (fluxVector[i + 1][j][iVar] - fluxVector[i - 2][j][iVar])
								+ c / ds * (fluxVector[i + 2][j][iVar] - fluxVector[i - 3][j][iVar]);

				rhsVector[iVar] = - rhsVector[iVar];
			}
			rhs[i][j] = rhsVector;
		}
	}
}

void Spatial_Derivative::Spatial_Derivative_WCNS_Y()
{
	double ds = dy;
	double a  = 75.0 / 64.0, b = -25.0 / 384.0, c = 3.0 / 640;
	
	VInt2D& marker = mesh->Get_Marker();
	for (int i = ist; i < ied - 1; i++)
	{
		for (int j = jst; j < jed - 1; j++)
		{
			if (marker[i][j] == 0) continue;

			VDouble rhsVector(num_of_prim_vars);
			for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
			{
				rhsVector[iVar] = a / ds * (fluxVector[i][j    ][iVar] - fluxVector[i][j - 1][iVar])
								+ b / ds * (fluxVector[i][j + 1][iVar] - fluxVector[i][j - 2][iVar])
								+ c / ds * (fluxVector[i][j + 2][iVar] - fluxVector[i][j - 3][iVar]);

				rhsVector[iVar] = - rhsVector[iVar];
			}
			rhs[i][j] = rhsVector;
		}
	}
}
