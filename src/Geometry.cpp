#include "Geometry.h"
#include "Global.h"
#include "2D_Euler_Solver.h"

Structured_Mesh* mesh;
int num_ghost_point;
int num_grid_point_x, num_grid_point_y;
int total_points_x,   total_points_y;
int num_half_point_x, num_half_point_y;
int Iw, Jw1, Jw2;
double dx, dy;

void Generate_Mesh()
{
//===================================================================================
	double hx  = 4.0, hy  = 2.0;
	double hx1 = 0.6, hy1 = 0.8, hy2 = 1.2;
	Iw  = hx1 / hx * (num_grid_point_x - 1); //������߽�ı��I����ʼ���Ϊ0�������������
	Jw1 = hy1 / hy * (num_grid_point_y - 1); //�����±߽�ı��J
	Jw2 = hy2 / hy * (num_grid_point_y - 1); //�����ϱ߽�ı��J
//===================================================================================
	num_ghost_point = 2;
	total_points_x = num_grid_point_x + 2 * num_ghost_point;
	total_points_y = num_grid_point_y + 2 * num_ghost_point;

	//num_half_point_x = num_grid_point_x - 1;   //������������Ԫ����������1
	//num_half_point_y = num_grid_point_y - 1;   //������������Ԫ����������1
	num_half_point_x = total_points_x - 1;   //������������Ԫ����������1
	num_half_point_y = total_points_x - 1;
	
	dx = hx / (num_grid_point_x - 1);
	dy = hy / (num_grid_point_y - 1);

	mesh = new Structured_Mesh(total_points_x, total_points_y);
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();
	vector< vector< int > >& marker = mesh->Get_Marker();

	mesh->Set_Dxdy(dx, dy);

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);
	for (int i = ist; i < ied; i++)
	{
		for (int j = jst; j < jed; j++)
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

void Set_Mesh_Dimension(int NI, int NJ)
{
	num_grid_point_x = NI; num_grid_point_y = NJ;
}

Structured_Mesh::Structured_Mesh(int NI, int NJ)
{
	this->dx = 0.0;
	this->dy = 0.0;
	this->NI = NI;
	this->NJ = NJ;
	Allocate_2D_Vector(grid_points, NI, NJ);
	Allocate_2D_Vector(marker, NI, NJ);
}

Point::Point()
{
	xPoint = 0.0;
	yPoint = 0.0;
}