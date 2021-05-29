#include "Global_Variables.h"
using namespace GLOBAL;

void GLOBAL::Init_Global_Param()
{
	max_num_of_steps = 10000;

	cfl_num = 0.1;
	time_step = 0.0; //ʱ�䲽��Ҫ�����������ֵȷ��������ֻ�ǳ�ʼ��

	grid_point_num_x = 101;
	grid_point_num_y = 51;
	ghost_point_num = 2;
	total_points_x = grid_point_num_x + 2 * ghost_point_num;
	total_points_y = grid_point_num_y + 2 * ghost_point_num;

	num_of_prim_vars = 4; //ԭʼ�������������Ʒ��̸���

	//������ʼ��
	qField.resize(num_of_prim_vars);
	for (int iVar = 0; iVar < num_of_prim_vars; iVar++)
	{
		qField[iVar].resize(total_points_x);
		for (int i = 0; i < total_points_x; i++)
		{
			qField[iVar][i].resize(total_points_y);
		}
	}
}

void GLOBAL::Generate_Mesh()
{
	double hx = 4.0, hy = 2.0;
	double hx1 = 0.6, hy1 = 0.8, hy2 = 1.2;

	Iw = hx1 / hx * (grid_point_num_x - 1); //������߽�ı��I����ʼ���Ϊ0�������������
	Jw1 = hy1 / hy * (grid_point_num_y - 1); //�����±߽�ı��J
	Jw2 = hy2 / hy * (grid_point_num_y - 1); //�����ϱ߽�ı��J

	//Iw  = Iw  + ghost_point_num;
	//Jw1 = Jw1 + ghost_point_num; 
	//Jw2 = Jw2 + ghost_point_num; 

	dx = hx / (grid_point_num_x - 1);
	dy = hy / (grid_point_num_y - 1);

	x_coord.resize(grid_point_num_x);
	y_coord.resize(grid_point_num_y);

	for (int iNode = 0; iNode < grid_point_num_x; ++iNode)
	{
		x_coord[iNode] = dx * iNode;
	}

	for (int jNode = 0; jNode < grid_point_num_y; ++jNode)
	{
		y_coord[jNode] = dy * jNode;
	}

	//���������Ƿ���������
	marker.resize(grid_point_num_x);
	for (int i = 0; i < grid_point_num_x; i++)
	{
		marker[i].resize(grid_point_num_y);
	}

	for (int i = 0; i < grid_point_num_x; i++)
	{
		for (int j = 0; j < grid_point_num_y; j++)
		{
			marker[i][j] = 1;

			double x_node = x_coord[i];
			double y_node = y_coord[j];

			if (x_node > hx1 && y_node > hy1 && y_node < hy2)
			{
				marker[i][j] = 0;		//�����Щ�������������棬���μӼ���
			}
		}
	}
}

void GLOBAL::Flow_Initialization()
{
	//��������ֵ
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
	//��߽磺���������
	for (int i = 0; i < ghost_point_num; i++)
	{
		for (int j = ghost_point_num; j < ghost_point_num + grid_point_num_y; j++)
		{
			qField[IR][i][j] = 1.0;
			qField[IU][i][j] = 3.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = 0.71429;
		}
	}

	//�ϱ߽磺outflow
	for (int i = ghost_point_num; i < ghost_point_num + grid_point_num_x; i++)
	{
		for (int j = ghost_point_num + grid_point_num_y; j < total_points_y; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = qField[IU][i][j - 1];
			qField[IV][i][j] = qField[IV][i][j - 1];
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}

	//�ұ߽磺outflow
	for (int i = ghost_point_num + grid_point_num_x; i < total_points_x; i++)
	{
		for (int j = ghost_point_num; j < ghost_point_num + grid_point_num_y; j++)
		{
			if (marker[i - 1][j] == 0) continue;

			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = qField[IU][i - 1][j];
			qField[IV][i][j] = qField[IV][i - 1][j];
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//�±߽�, no-slip wall
	for (int i = ghost_point_num; i < ghost_point_num + grid_point_num_x; i++)
	{
		for (int j = ghost_point_num - 1; j >= 0; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//������, no-slip wall
	for (int i = ghost_point_num + Iw; i < ghost_point_num + Iw + ghost_point_num; i++)
	{
		for (int j = ghost_point_num + Jw1; j < ghost_point_num + Jw2; j++)
		{
			qField[IR][i][j] = qField[IR][i - 1][j];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i - 1][j];
		}
	}

	//������, no-slip wall
	for (int i = ghost_point_num + Iw; i < ghost_point_num + grid_point_num_x; i++)
	{
		//for (int j = Jw2; j < ghost_point_num + Jw2; j++)
		for (int j = ghost_point_num + Jw2 - 1; j >= Jw2; j--)
		{
			qField[IR][i][j] = qField[IR][i][j + 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j + 1];
		}
	}

	//������, no-slip wall
	for (int i = ghost_point_num + Iw; i < ghost_point_num + grid_point_num_x; i++)
	{
		for (int j = ghost_point_num + Jw1; j < ghost_point_num + Jw1 + ghost_point_num; j++)
		{
			qField[IR][i][j] = qField[IR][i][j - 1];
			qField[IU][i][j] = 0.0;
			qField[IV][i][j] = 0.0;
			qField[IP][i][j] = qField[IP][i][j - 1];
		}
	}
}