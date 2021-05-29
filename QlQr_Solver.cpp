#include <iostream>
#include "QlQr_Solver.h"
#include "Global.h"
#include "Geometry.h"
#include "2D_Euler_Solver.h"

using namespace GLOBAL;

void GLOBAL::Solve_QlQr()
{
	auto * half_node_q = new QlQr_Solver();

	half_node_q->Solve_QlQr();

	delete half_node_q;
}

QlQr_Solver::QlQr_Solver()
{
	qField1.resize(num_of_prim_vars);
	qField2.resize(num_of_prim_vars);

	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		Allocate_2D_Vector(qField1[iVar], num_half_point_x, num_half_point_y);
		Allocate_2D_Vector(qField2[iVar], num_half_point_x, num_half_point_y);
	}
}

void QlQr_Solver::Solve_QlQr()
{
	if (method_of_half_q == 1)
	{
		this->QlQr_MUSCL();
	}
	else if (method_of_half_q == 2)
	{
		//WENO������а�ڵ��ֵ
	}
	else if (method_of_half_q == 3)
	{
		this->QlQr_WCNS();
	}
	else
	{
		cout << "��ڵ����ֵ��ֵ�����������飡" << endl;
	}
}

void QlQr_Solver::QlQr_MUSCL()
{
	//����x������в�ֵ
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		vector< vector< double > >& qv = qField[iVar];
		for (int j = 0; j < num_half_point_y; j++)
		{
			for (int i = 0; i < num_half_point_x; i++)
			{
				double du_p1 = qv[i + 1][j] - qv[i    ][j];
				double du_m1 = qv[i    ][j] - qv[i - 1][j];
				double du_p3 = qv[i + 2][j] - qv[i + 1][j];

				double ita_m1_p = du_p1 / du_m1;
				double ita_p1_m = du_m1 / du_p1;
				double ita_p3_m = du_p1 / du_p3;
				double ita_p1_p = du_p3 / du_p1;
				
				double fai1 = Limiter_Function(ita_m1_p);
				double fai2 = Limiter_Function(ita_p1_m);
				double fai3 = Limiter_Function(ita_p3_m);
				double fai4 = Limiter_Function(ita_p1_p);

				qField1[iVar][i][j] = qv[i    ][j] + 1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1
																+ (1 + muscl_k) * fai2 * du_p1);

				qField2[iVar][i][j] = qv[i + 1][j] - 1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3
																+ (1 + muscl_k) * fai4 * du_p1);
			}
		}
	}
}

void QlQr_Solver::QlQr_MUSCL_Y()
{
	//��y������в�ֵ
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		vector< vector< double > >& qv = qField[iVar];
		for (int i = 0; i < num_half_point_x; i++)
		{
			for (int j = 0; j < num_half_point_y; j++)
			{
				double du_p1 = qv[i][j + 1] - qv[i][j    ];
				double du_m1 = qv[i][j    ] - qv[i][j - 1];
				double du_p3 = qv[i][j + 2] - qv[i][j + 1];

				double ita_m1_p = du_p1 / du_m1;
				double ita_p1_m = du_m1 / du_p1;
				double ita_p3_m = du_p1 / du_p3;
				double ita_p1_p = du_p3 / du_p1;
				
				double fai1 = Limiter_Function(ita_m1_p);
				double fai2 = Limiter_Function(ita_p1_m);
				double fai3 = Limiter_Function(ita_p3_m);
				double fai4 = Limiter_Function(ita_p1_p);

				qField1[iVar][i][j] = qv[i][j	] + 1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1
															   + (1 + muscl_k) * fai2 * du_p1);

				qField2[iVar][i][j] = qv[i][j + 1] - 1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3
																+ (1 + muscl_k) * fai4 * du_p1);
			}
		}
	}
}

void QlQr_Solver::QlQr_WCNS()
{

}

double QlQr_Solver::Limiter_Function( double ita )
{
	return vanleer_limiter(ita, 1.0);
}

//����������
double GLOBAL::minmod_limiter(double a, double b)
{
	if (a * b <= 0)
		return 0;
	else
	{
		if ((fabs(a) - fabs(b)) > 0)
			return b;
		else
			return a;
	}
}

double GLOBAL::vanleer_limiter(double a, double)
{
	return (a + abs(a)) / (1.0 + abs(a));
}

double GLOBAL::superbee_limiter(double a, double)
{
	double tmp1 = min(2.0 * a, 1.0);
	double tmp2 = min(a, 2.0);

	return max(0.0, max(tmp1, tmp2));
}
