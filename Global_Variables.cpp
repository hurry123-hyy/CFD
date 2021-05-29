#include "Global_Variables.h"
using namespace GLOBAL;

void GLOBAL::Init_Global_Param()
{
	num_of_prim_vars = 4;		//ԭʼ�������������Ʒ��̸���

	max_num_of_steps = 10000;

	cfl_num   = 0.1;
	time_step = 0.0;			//ʱ�䲽��Ҫ�����������ֵȷ��������ֻ�ǳ�ʼ��

	method_of_half_q  = 1;		//1-MUSCL,		2-WENO,		3-WCNS
	muscl_k			  = 0.0;	//0.0-����ӭ��ƫ�ã�		1/3-����ӭ��ƫ��
	method_of_limiter = 1;		//1-vanleer,	2-minmod,	3-superbee	
	method_of_flux    = 1;		//1-Roe,		2-WENO,		3-WCNS
	entropy_fix_coeff = 0.01;	//Roe��ʽ������ϵ��epsilon

	num_grid_point_x = 101;
	num_grid_point_y = 51;
	num_ghost_point  = 2;

	total_points_x = num_grid_point_x + 2 * num_ghost_point;
	total_points_y = num_grid_point_y + 2 * num_ghost_point;

	num_half_point_x = num_grid_point_x - 1;   //������������Ԫ����������1
	num_half_point_y = num_grid_point_y - 1;   //������������Ԫ����������1

	//������ʼ��
	qField.resize(num_of_prim_vars);
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		Allocate_2D_Vector(qField[iVar], total_points_x, total_points_y);
	}
}

void GLOBAL::Generate_Mesh()
{
	double hx  = 4.0, hy  = 2.0;
	double hx1 = 0.6, hy1 = 0.8, hy2 = 1.2;

	Iw  = hx1 / hx * (num_grid_point_x - 1); //������߽�ı��I����ʼ���Ϊ0�������������
	Jw1 = hy1 / hy * (num_grid_point_y - 1); //�����±߽�ı��J
	Jw2 = hy2 / hy * (num_grid_point_y - 1); //�����ϱ߽�ı��J

	dx = hx / (num_grid_point_x - 1);
	dy = hy / (num_grid_point_y - 1);

	mesh = new Structured_Mesh();
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();
	vector< vector< int > >& marker = mesh->Get_Marker();

	for (int i = 0; i < num_grid_point_x; i++)
	{
		for (int j = 0; j < num_grid_point_y; j++)
		{
			double x_node = i * dx;
			double y_node = j * dy;

			grid_points[i][j].Set_Point_Coord(x_node, y_node);

			marker[i][j] = 1;
			if (x_node > hx1 && y_node > hy1 && y_node < hy2)
			{
				marker[i][j] = 0;		//�����Щ�������������棬���μӼ���
			}
		}
	}
}

Structured_Mesh::Structured_Mesh()
{
	Allocate_2D_Vector(grid_points, num_grid_point_x, num_grid_point_y);
	Allocate_2D_Vector(marker,      num_grid_point_x, num_grid_point_y);
}

void GLOBAL::Flow_Initialization()
{
	//��������ֵ
	vector< vector< int > >& marker = mesh->Get_Marker();
	for (int i = 0; i < total_points_x; i++)
	{
		for (int j = 0; j < total_points_y; j++)
		{
			if (marker[i][j] == 0) continue;

			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}
}

void GLOBAL::Compute_Boundary()
{
	vector< vector< int > >& marker = mesh->Get_Marker();
	//��߽磺���������
	for (int i = 0; i < num_ghost_point; i++)
	{
		for (int j = num_ghost_point; j < num_ghost_point + num_grid_point_y; j++)
		{
			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}

	//�ϱ߽磺outflow
	for (int i = num_ghost_point; i < num_ghost_point + num_grid_point_x; i++)
	{
		for (int j = num_ghost_point + num_grid_point_y; j < total_points_y; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = qField[IU][i][j - 1];
			qField[IV][i][j] = qField[IV][i][j - 1];
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}

	//�ұ߽磺outflow
	for (int i = num_ghost_point + num_grid_point_x; i < total_points_x; i++)
	{
		for (int j = num_ghost_point; j < num_ghost_point + num_grid_point_y; j++)
		{
			if (marker[i - 1][j] == 0) continue;

			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = qField[IU][i - 1][j];
			qField[IV][i][j] = qField[IV][i - 1][j];
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//�±߽�, no-slip wall
	for (int i = num_ghost_point; i < num_ghost_point + num_grid_point_x; i++)
	{
		for (int j = num_ghost_point - 1; j >= 0; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//������, no-slip wall
	for (int i = num_ghost_point + Iw; i < num_ghost_point + Iw + num_ghost_point; i++)
	{
		for (int j = num_ghost_point + Jw1; j < num_ghost_point + Jw2; j++)
		{
			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//������, no-slip wall
	for (int i = num_ghost_point + Iw; i < num_ghost_point + num_grid_point_x; i++)
	{
		//for (int j = Jw2; j < num_ghost_point + Jw2; j++)
		for (int j = num_ghost_point + Jw2 - 1; j >= Jw2; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//������, no-slip wall
	for (int i = num_ghost_point + Iw; i < num_ghost_point + num_grid_point_x; i++)
	{
		for (int j = num_ghost_point + Jw1; j < num_ghost_point + Jw1 + num_ghost_point; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}
}