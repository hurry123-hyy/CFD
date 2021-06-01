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
int num_of_RK_stages;
VDouble2D RK_Coeff;

void Init_Global_Param()
{
	num_of_prim_vars = 4;		//ԭʼ�������������Ʒ��̸���

	max_num_of_steps = 10000;

	cfl_num   = 0.6;
	time_step = 0.0;			//ʱ�䲽��Ҫ�����������ֵȷ��������ֻ�ǳ�ʼ��
	physical_time     = 0.0;
	max_simu_time	  = 0.2;
	method_of_half_q  = 2;		//1-MUSCL,	  2-WENO(����ֵ),   3-WCNS
	muscl_k			  = 0.0;	//0.0-����ӭ��ƫ�ã�		    1/3-����ӭ��ƫ��
	method_of_limiter = 1;		//1-vanleer,  2-minmod,		    3-superbee	
	method_of_flux    = 3;		//1-Roe,	  2-Steger Warming  3-WENO,		  4-WCNS
	entropy_fix_coeff = 0.01;	//Roe��ʽ������ϵ��epsilon

	num_grid_point_x = 241;
	num_grid_point_y = 61;

	solve_direction  = 'x';

	residual_output_steps = 2;		//�в�����������
	flow_save_steps		  = 200;	//��������������
	converge_criterion	  = 1e-8;	//�в�������׼
	tec_file_name		  = "../../flow.plt";
	
	num_of_RK_stages	= 3;
	RK_Coeff			= { {1.0, 0.0, 1.0},{3.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0},{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0} };
}

void Load_Q()
{
	qField = qField_N1;
}

void Set_Solve_Direction(char direction)
{
	solve_direction = direction;
	if (direction=='x')
	{		
		qField_N0 = qField;		//RK��ʽ���һ���ǰʱ�䲽Q
		qField_N1 = qField;		//RK��ʽ��ڶ����һstage��Q
	}
	else if(direction == 'y')
	{
		//��ؼ��ĵ㣺y��������qField��x��������qField����ͬ�ģ����������ӷ��ѷ��Ĺؼ���
		//Ҳ����Ϊ��ˣ��ڼ���y����ʱ���������¼���߽�����
		qField	  = qField_N0;		//��qField��ԭΪ����x����֮ǰ��Qֵ��y��������ԭ����qField������rhs
		qField_N0 = qField_N1;		//RK��ʽ���һ���x������stage�ƽ�����������Qֵ
		qField_N1 = qField;			//qField��qField_N1��RK�ƽ��еĹؼ���������Ȼ��Ϊ����x����֮ǰ��Qֵ
	}
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

bool Need_Stop_Iteration()
{
	return stop_by_residual || physical_time >= max_simu_time || current_step == max_num_of_steps;
}

bool IsNaN(VDouble& data)
{
	bool flag = 0;
	for (int i = 0; i < data.size(); i++)
	{
		if (data[i] != data[i])
		{
			flag = 1;
			break;
		}
	}
	return flag;
}