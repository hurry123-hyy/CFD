#include "2D_Euler_Solver.h"
#include "QlQr_Solver.h"
#include "Global.h"
#include "Geometry.h"

int num_of_prim_vars;
int current_step, max_num_of_steps;
double cfl_num, time_step, physical_time, max_simu_time;
int method_of_half_q;
int method_of_limiter;
int method_of_flux;
double muscl_k;
double entropy_fix_coeff;
char solve_direction;
int residual_output_steps;
int flow_save_steps;
double converge_criterion;
string tec_file_name;

void Init_Global_Param()
{
	num_of_prim_vars = 4;		//ԭʼ�������������Ʒ��̸���

	max_num_of_steps = 10000;

	cfl_num   = 0.5;
	time_step = 0.0;			//ʱ�䲽��Ҫ�����������ֵȷ��������ֻ�ǳ�ʼ��
	physical_time     = 0.0;
	max_simu_time	  = 0.2;
	method_of_half_q  = 1;		//1-MUSCL,		2-WENO,		3-WCNS
	muscl_k			  = 0.0;	//0.0-����ӭ��ƫ�ã�		1/3-����ӭ��ƫ��
	method_of_limiter = 1;		//1-vanleer,	2-minmod,	3-superbee	
	method_of_flux    = 1;		//1-Roe,		2-WENO,		3-WCNS
	entropy_fix_coeff = 0.01;	//Roe��ʽ������ϵ��epsilon

	num_grid_point_x = 241;
	num_grid_point_y = 61;

	solve_direction  = 'x';

	residual_output_steps = 20;		//�в�����������
	flow_save_steps		  = 200;	//��������������
	converge_criterion	  = 1e-8;	//�в�������׼
	tec_file_name		  = "../../flow.plt";
}

void Init_Flow()
{
	//������ʼ��
	Allocate_3D_Vector(qField,	  total_points_x, total_points_y, num_of_prim_vars);
	Allocate_3D_Vector(qField_N1, total_points_x, total_points_y, num_of_prim_vars);
	//Allocate_3D_Vector(qField_N2, total_points_x, total_points_y, num_of_prim_vars);
	//Allocate_3D_Vector(qField_N3, total_points_x, total_points_y, num_of_prim_vars);

	//��������ֵ
	VInt2D& marker = mesh->Get_Marker();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			if (marker[i][j] == 0) continue;

			qField[i][j][IR] = 1.0;
			qField[i][j][IU] = 3.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = 0.71429;
		}
	}

	qField_N1 = qField;
}

void Init_Flow_Double_Mach()
{
	//������ʼ��
	Allocate_3D_Vector(qField, total_points_x, total_points_y, num_of_prim_vars);
	Allocate_3D_Vector(qField_N1, total_points_x, total_points_y, num_of_prim_vars);

	//��������ֵ
	VInt2D& marker = mesh->Get_Marker();
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();
	
	double gama = 1.4;
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			double x_node, y_node;

			grid_points[i][j].Get_Point_Coord(x_node, y_node);
			if (x_node <= (1.0 / 6.0 + y_node / tan(PI / 3)))
			{
				qField[i][j][IR] =  8.0;
				qField[i][j][IU] =  8.25 * cos(PI / 6);
				qField[i][j][IV] = -8.25 * sin(PI / 6);
				qField[i][j][IP] = 116.5;
			}
			else //if (x_node > (1.0 / 6.0 + y_node / tan(PI / 3)))
			{
				qField[i][j][IR] = 1.4;
				qField[i][j][IU] = 0.0;
				qField[i][j][IV] = 0.0;
				qField[i][j][IP] = 1.0;
			}
		}
	}

	qField_N1 = qField;
}

void Compute_Boundary()
{
	VInt2D& marker = mesh->Get_Marker();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	//��߽磺���������
	for (int i = 0; i < ist; i++)
	{
		for (int j = jst; j < jed; j++)
		{
			qField[i][j][IR] = 1.0;
			qField[i][j][IU] = 3.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = 0.71429;
		}
	}

	//�ϱ߽磺outflow
	for (int i = ist; i < ied; i++)
	{
		for (int j = jed; j < total_points_y; j++)
		{
			qField[i][j][IR] = qField[i][j - 1][IR];
			qField[i][j][IU] = qField[i][j - 1][IU];
			qField[i][j][IV] = qField[i][j - 1][IV];
			qField[i][j][IP] = qField[i][j - 1][IP];
		}
	}

	//�ұ߽磺outflow
	for (int j = jst; j < jed; j++)
	{
		for (int i = ied; i < total_points_x; i++)
		{
			if (marker[i - 1][j] == 0) continue;

			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//�±߽�, no-slip wall
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst - 1; j >= 0; j--)
		{
			qField[i][j][IR] = qField[i][j + 1][IR];
			qField[i][j][IU] = 0.0; 
			qField[i][j][IV] = 0.0; 
			qField[i][j][IP] = qField[i][j + 1][IP];
		}
	}

	//������, no-slip wall
	for (int j = jst + Jw1; j < jst + Jw2; j++)
	{
		for (int i = ist + Iw; i < ist + Iw + num_ghost_point; i++)
		{
			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//������, no-slip wall
	for (int i = ist + Iw; i < ied; i++)
	{
		//for (int j = Jw2; j < num_ghost_point + Jw2; j++)
		for (int j = jst + Jw2 - 1; j >= Jw2; j--)
		{
			qField[i][j][IR] = qField[i][j + 1][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j + 1][IP];
		}
	}

	//������, no-slip wall
	for (int i = ist + Iw; i < ied; i++)
	{
		for (int j = jst + Jw1; j < jst + Jw1 + num_ghost_point; j++)
		{
			qField[i][j][IR] = qField[i][j - 1][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j - 1][IP];
		}
	}
}

void Compute_Boundary_Double_Mach()
{
	VInt2D& marker = mesh->Get_Marker();
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	double ssw = 10.0 / sin(PI / 3);
	double xsw = 1.0 / 6 + 1.0 / tan(PI / 3) + ssw * physical_time;

	//��߽磺inflow
	//for (int i = 0; i < ist; i++)
	for (int i = ist; i <= ist; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			qField[i][j][IR] =  8.0;
			qField[i][j][IU] =  8.25 * cos(PI / 6);
			qField[i][j][IV] = -8.25 * sin(PI / 6);
			qField[i][j][IP] = 116.5;
		}
	}

	//�ϱ߽�
	//for (int j = jed; j < total_points_y; j++)
	for (int j = jed - 1; j <= jed - 1; j++)
	{
		for (int i = ist; i < ied; i++)
		{
			double x_node, y_node;
			grid_points[i][j].Get_Point_Coord(x_node, y_node);
			if (x_node >= 0 && x_node <= xsw)
			{
				qField[i][j][IR] =  8.0;
				qField[i][j][IU] =  8.25 * cos(PI / 6);
				qField[i][j][IV] = -8.25 * sin(PI / 6);
				qField[i][j][IP] = 116.5;
			}
			else if(x_node > xsw && x_node <= 4)
			{
				qField[i][j][IR] = 1.4;
				qField[i][j][IU] = 0.0;
				qField[i][j][IV] = 0.0;
				qField[i][j][IP] = 1.0;
			}
		}
	}

	//�ұ߽磺outflow
	for (int j = jst; j < jed; j++)
	{
		//for (int i = ied; i < total_points_x; i++)
		for (int i = ied - 1; i <= ied - 1; i++)
		{
			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//�±߽�
	for (int i = ist; i < ied; i++)
	{
		//for (int j = jst - 1; j >= 0; j--)
		for (int j = jst; j >= jst; j--)
		{
			double x_node, y_node;
			grid_points[i][j].Get_Point_Coord(x_node, y_node);
			if (x_node >= 0 && x_node <= 1.0 / 6)
			{
				qField[i][j][IR] = 8.0;
				qField[i][j][IU] = 8.25 * cos(PI / 6);
				qField[i][j][IV] = -8.25 * sin(PI / 6);
				qField[i][j][IP] = 116.5;
			}
			else if (x_node > 1.0 / 6 && x_node <= 4)
			{				
				qField[i][j][IR] = qField[i][j + 1][IR];
				qField[i][j][IU] = qField[i][j + 1][IU];
				qField[i][j][IV] = 0.0;
				qField[i][j][IP] = qField[i][j + 1][IP];
			}
		}
	}
}

void Load_Q()
{
	qField = qField_N1;
}

void Set_Solve_Direction(char direction)
{
	solve_direction = direction;
}

void Primitive_To_Conservative( VDouble & primitive, VDouble &conservative )
{
	double gama = 1.4;
	double rho, u, v, p;
	rho = primitive[IR];
	u = primitive[IU];
	v = primitive[IV];
	p = primitive[IP];

	conservative[IR] = rho;
	conservative[IU] = rho * u;
	conservative[IV] = rho * v;
	conservative[IP] = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
}

void Conservative_To_Primitive(VDouble& conservative, VDouble &primitive)
{
	double gama = 1.4;
	double rho		= conservative[IR];
	double rho_u	= conservative[IU];
	double rho_v	= conservative[IV];
	double E		= conservative[IP];

	double u = rho_u / rho;
	double v = rho_v / rho;
	double p = (gama - 1) * (E - 0.5 * rho * (u * u + v * v));

	primitive[IR] = rho;
	primitive[IU] = u;
	primitive[IV] = v;
	primitive[IP] = p;
}

//�����ı�ŷ�Χ
void Get_IJK_Region(int& ist, int& ied, int& jst, int& jed)
{
	ist = num_ghost_point;
	ied = num_ghost_point + num_grid_point_x;

	jst = num_ghost_point;
	jed = num_ghost_point + num_grid_point_y;
}