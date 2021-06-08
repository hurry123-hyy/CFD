#include <cmath>
#include "2D_Euler_Solver.h"
#include "Geometry.h"
#include "QlQr_Solver.h"

void Compute_Boundary()
{
	if (global_case_id == 1)
	{
		Compute_Boundary_Double_Mach();
	}
	else if (global_case_id == 2)
	{
		Compute_Boundary_Blunt_Body();
	}
}

void Compute_Boundary_Blunt_Body()
{
	VInt2D& marker = mesh->Get_Marker_Q();
	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	//�ϱ߽磺outflow
	for (int i = ist; i <= ied; i++)
	{
		for (int j = jed; j <= jed + 2; j++)
		{
			qField[i][j][IR] = qField[i][j - 1][IR];
			qField[i][j][IU] = qField[i][j - 1][IU];
			qField[i][j][IV] = qField[i][j - 1][IV];
			qField[i][j][IP] = qField[i][j - 1][IP];
		}
	}

	//�ұ߽磺outflow
	for (int j = jst; j <= jed; j++)
	{
		for (int i = ied; i <= ied + 2; i++)
		{
			if (marker[i][j] == 0) continue;

			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//�±߽�, no-slip wall
	for (int i = ist; i <= ied; i++)
	{
		for (int j = jst; j >= 0; j--)
		{
			qField[i][j][IR] = qField[i][j + 1][IR];
			//qField[i][j][IU] = qField[i][j + 1][IU];
			//qField[i][j][IV] = qField[i][j + 1][IV];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j + 1][IP];
		}

		//qField[i][jst    ][IU] = 0.0;
		//qField[i][jst - 1][IU] = -qField[i][jst + 1][IU];
		//qField[i][jst - 2][IU] = -qField[i][jst + 1][IU];

		//qField[i][jst    ][IV] = 0.0;
		//qField[i][jst - 1][IV] = -qField[i][jst + 1][IV];
		//qField[i][jst - 2][IV] = -qField[i][jst + 1][IV];
	}

	//��߽磺���������    ��������߽��������÷����±߽���棬ȷ�����½ǵ�����Ϊ��ڱ߽�
	for (int i = 0; i <= ist; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			qField[i][j][IR] = 1.0;
			qField[i][j][IU] = 3.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = 0.71429;
		}
	}

	//������, no-slip wall
	for (int j = jst + Jw1; j <= jst + Jw2; j++)
	{
		for (int i = ist + Iw; i <= ist + Iw + num_ghost_point; i++)
		{
			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i - 1][j][IP];
		}

		//qField[ist + Iw    ][j][IU] = 0.0;
		//qField[ist + Iw + 1][j][IU] = -qField[ist + Iw - 1][j][IU];
		//qField[ist + Iw + 2][j][IU] = -qField[ist + Iw - 1][j][IU];

		//qField[ist + Iw    ][j][IV] = 0.0;
		//qField[ist + Iw + 1][j][IV] = -qField[ist + Iw - 1][j][IV];
		//qField[ist + Iw + 2][j][IV] = -qField[ist + Iw - 1][j][IV];
	}

	//������, no-slip wall
	for (int i = ist + Iw; i <= ied; i++)
	{
		for (int j = jst + Jw2; j >= Jw2; j--)
		{
			qField[i][j][IR] = qField[i][j + 1][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j + 1][IP];
		}

		//qField[i][jst + Jw2    ][IU] = 0.0;
		//qField[i][jst + Jw2 - 1][IU] = -qField[i][jst + Jw2 + 1][IU];
		//qField[i][jst + Jw2 - 2][IU] = -qField[i][jst + Jw2 + 1][IU];

		//qField[i][jst + Jw2    ][IV] = 0.0;
		//qField[i][jst + Jw2 - 1][IV] = -qField[i][jst + Jw2 + 1][IV];
		//qField[i][jst + Jw2 - 2][IV] = -qField[i][jst + Jw2 + 1][IV];
	}

	//������, no-slip wall
	for (int i = ist + Iw; i <= ied; i++)
	{
		for (int j = jst + Jw1; j <= jst + Jw1 + num_ghost_point; j++)
		{
			qField[i][j][IR] = qField[i][j - 1][IR];
			qField[i][j][IU] = 0.0;
			qField[i][j][IV] = 0.0;
			qField[i][j][IP] = qField[i][j - 1][IP];
		}

		//qField[i][jst + Jw1    ][IU] = 0.0;
		//qField[i][jst + Jw1 + 1][IU] = -qField[i][jst + Jw1 - 1][IU];
		//qField[i][jst + Jw1 + 2][IU] = -qField[i][jst + Jw1 - 1][IU];

		//qField[i][jst + Jw1    ][IV] = 0.0;
		//qField[i][jst + Jw1 + 1][IV] = -qField[i][jst + Jw1 - 1][IV];
		//qField[i][jst + Jw1 + 2][IV] = -qField[i][jst + Jw1 - 1][IV];
	}

	//�ǵ㴦�����ϣ����òο��ǵ���
	int II = ist + Iw;
	int JJ = jst + Jw2;
	//1�㣬�ԽǷ���
	qField[II + 2][JJ - 2][IR] =  qField[II][JJ][IR];
	qField[II + 2][JJ - 2][IU] = -qField[II][JJ][IU];
	qField[II + 2][JJ - 2][IV] = -qField[II][JJ][IV];
	qField[II + 2][JJ - 2][IP] =  qField[II][JJ][IP];

	//2��, x����ȡֵ
	qField[II + 1][JJ - 2][IR] =  qField[II - 1][JJ - 2][IR];
	qField[II + 1][JJ - 2][IU] = -qField[II - 1][JJ - 2][IU];
	qField[II + 1][JJ - 2][IV] = -qField[II - 1][JJ - 2][IV];
	qField[II + 1][JJ - 2][IP] =  qField[II - 1][JJ - 2][IP];

	//3��, y����ȡֵ
	qField[II + 2][JJ - 1][IR] =  qField[II + 2][JJ + 1][IR];
	qField[II + 2][JJ - 1][IU] = -qField[II + 2][JJ + 1][IU];
	qField[II + 2][JJ - 1][IV] = -qField[II + 2][JJ + 1][IV];
	qField[II + 2][JJ - 1][IP] =  qField[II + 2][JJ + 1][IP];

	//4�㣬 x����ȡֵ��y����ȡֵ������ƽ��
	qField[II + 1][JJ - 1][IR] = 0.5 * (qField[II - 1][JJ - 1][IR] + qField[II + 1][JJ + 1][IR]);
	qField[II + 1][JJ - 1][IU] = 0.5 * (qField[II - 1][JJ - 1][IU] + qField[II + 1][JJ + 1][IU]);
	qField[II + 1][JJ - 1][IV] = 0.5 * (qField[II - 1][JJ - 1][IV] + qField[II + 1][JJ + 1][IV]);
	qField[II + 1][JJ - 1][IP] = 0.5 * (qField[II - 1][JJ - 1][IP] + qField[II + 1][JJ + 1][IP]);
	
	//=========================
	//�ǵ㴦�����£����òο��ǵ���
	II = ist + Iw;
	JJ = jst + Jw1;

	//1�㣬�ԽǷ���
	qField[II + 2][JJ + 2][IR] =  qField[II][JJ][IR];
	qField[II + 2][JJ + 2][IU] = -qField[II][JJ][IU];
	qField[II + 2][JJ + 2][IV] = -qField[II][JJ][IV];
	qField[II + 2][JJ + 2][IP] =  qField[II][JJ][IP];

	//2��, x����ȡֵ
	qField[II + 1][JJ + 2][IR] =  qField[II - 1][JJ + 2][IR];
	qField[II + 1][JJ + 2][IU] = -qField[II - 1][JJ + 2][IU];
	qField[II + 1][JJ + 2][IV] = -qField[II - 1][JJ + 2][IV];
	qField[II + 1][JJ + 2][IP] =  qField[II - 1][JJ + 2][IP];

	//3��, y����ȡֵ
	qField[II + 2][JJ + 1][IR] =  qField[II + 2][JJ - 1][IR];
	qField[II + 2][JJ + 1][IU] = -qField[II + 2][JJ - 1][IU];
	qField[II + 2][JJ + 1][IV] = -qField[II + 2][JJ - 1][IV];
	qField[II + 2][JJ + 1][IP] =  qField[II + 2][JJ - 1][IP];

	//4�㣬 x����ȡֵ��y����ȡֵ������ƽ��
	qField[II + 1][JJ + 1][IR] = 0.5 * (qField[II - 1][JJ + 1][IR] + qField[II + 1][JJ - 1][IR]);
	qField[II + 1][JJ + 1][IU] = 0.5 * (qField[II - 1][JJ + 1][IU] + qField[II + 1][JJ - 1][IU]);
	qField[II + 1][JJ + 1][IV] = 0.5 * (qField[II - 1][JJ + 1][IV] + qField[II + 1][JJ - 1][IV]);
	qField[II + 1][JJ + 1][IP] = 0.5 * (qField[II - 1][JJ + 1][IP] + qField[II + 1][JJ - 1][IP]);
}

void Compute_Boundary_Double_Mach()
{
	VInt2D& marker = mesh->Get_Marker_Q();
	vector< vector< Point > >& grid_points = mesh->Get_Grid_Points();

	int ist, ied, jst, jed;
	Get_IJK_Region(ist, ied, jst, jed);

	double ssw = 10.0 / sin(PI / 3);
	double xsw = 1.0 / 6 + 1.0 / tan(PI / 3) + ssw * physical_time;

	//��߽磺inflow
	//for (int i = 0; i < ist; i++)
	for (int i = 0; i <= ist; i++)
	{
		for (int j = jst; j <= jed; j++)
		{
			qField[i][j][IR] = 8.0;
			qField[i][j][IU] = 8.25 * cos(PI / 6);
			qField[i][j][IV] = -8.25 * sin(PI / 6);
			qField[i][j][IP] = 116.5;
		}
	}

	//�ϱ߽�
	//for (int j = jed + 1; j <= jed + 2; j++)
	for (int j = jed; j <= jed + 2; j++)
	{
		for (int i = ist; i <= ied; i++)
		{
			double x_node, y_node;
			grid_points[i][j].Get_Point_Coord(x_node, y_node);
			if (x_node >= 0 && x_node <= xsw)
			{
				qField[i][j][IR] = 8.0;
				qField[i][j][IU] = 8.25 * cos(PI / 6);
				qField[i][j][IV] = -8.25 * sin(PI / 6);
				qField[i][j][IP] = 116.5;
			}
			else if (x_node > xsw && x_node <= 4)
			{
				qField[i][j][IR] = 1.4;
				qField[i][j][IU] = 0.0;
				qField[i][j][IV] = 0.0;
				qField[i][j][IP] = 1.0;
			}
		}
	}

	//�ұ߽磺outflow
	for (int j = jst; j <= jed; j++)
	{
		//for (int i = ied + 1; i <= ied + 2; i++)
		for (int i = ied; i <= ied + 2; i++)
		{
			qField[i][j][IR] = qField[i - 1][j][IR];
			qField[i][j][IU] = qField[i - 1][j][IU];
			qField[i][j][IV] = qField[i - 1][j][IV];
			qField[i][j][IP] = qField[i - 1][j][IP];
		}
	}

	//�±߽�
	for (int i = ist; i <= ied; i++)
	{
		//for (int j = jst - 1; j >= 0; j--)
		for (int j = jst; j >= 0; j--)
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
