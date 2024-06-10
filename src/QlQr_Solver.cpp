//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>

#include "QlQr_Solver.h"
#include "Blunt_Solver.h"
#include "Geometry.h"

void Solve_QlQr() {
    QlQr_Solver* half_node_q = new QlQr_Solver();

    half_node_q->Solve_QlQr();

    delete half_node_q;
}

VDouble3D qField;
VDouble3D qField1;
VDouble3D qField2;
VDouble3D qField_N0;
VDouble3D qField_N1;

QlQr_Solver::QlQr_Solver() {
    Get_IJK_Region(ist, ied, jst, jed);
    Allocate_3D_Vector(qField1, num_half_point_x, num_half_point_y, num_of_primitive_vars);
    Allocate_3D_Vector(qField2, num_half_point_x, num_half_point_y, num_of_primitive_vars);
}

void QlQr_Solver::Solve_QlQr() {
    this->QlQr_MUSCL();
}

void QlQr_Solver::QlQr_MUSCL() {
    if (solve_direction == 'x') {
        this->QlQr_MUSCL_X();
        this->Boundary_QlQr_MUSCL_X();
    } else if (solve_direction == 'y') {
        this->QlQr_MUSCL_Y();
        this->Boundary_QlQr_MUSCL_Y();
    } else {
        cout << "MUSCL插值出错，请检查！" << endl;
    }
}

void QlQr_Solver::Boundary_QlQr_MUSCL_X() {
    // VInt2D& marker = mesh->Get_Marker();

    for (int j = jst; j <= jed; j++) {
        for (int iVar = 0; iVar < num_of_primitive_vars; iVar++) {
            qField1[0][j][iVar] = qField[0][j][iVar];  // 虚拟点i=0, ied + 1的值还没有
            qField2[0][j][iVar] = qField[1][j][iVar];
            // qField1[1][j][iVar] = qField[1][j][iVar];
            // qField2[1][j][iVar] = qField[2][j][iVar];

            // qField1[ied    ][j][iVar] = qField[ied    ][j][iVar];
            // qField2[ied    ][j][iVar] = qField[ied + 1][j][iVar];
            qField1[ied + 1][j][iVar] = qField[ied + 1][j][iVar];
            qField2[ied + 1][j][iVar] = qField[ied + 2][j][iVar];
        }
    }

    for (int j = jst + Jw1 + 1; j <= jst + Jw2 - 1; j++) {
        for (int iVar = 0; iVar < num_of_primitive_vars; iVar++) {
            qField1[ist + Iw][j][iVar] = qField[ist + Iw][j][iVar];
            qField2[ist + Iw][j][iVar] = qField[ist + Iw + 1][j][iVar];

            qField1[ist + Iw + 1][j][iVar] = qField[ist + Iw + 1][j][iVar];
            qField2[ist + Iw + 1][j][iVar] = qField[ist + Iw + 2][j][iVar];
        }
    }
}

void QlQr_Solver::Boundary_QlQr_MUSCL_Y() {
    // VInt2D& marker = mesh->Get_Marker();

    for (int i = ist; i <= ied; i++) {
        for (int iVar = 0; iVar < num_of_primitive_vars; iVar++)  // 虚拟点j=0, jed + 1的值还没有
        {
            qField1[i][0][iVar] = qField[i][0][iVar];
            qField2[i][0][iVar] = qField[i][1][iVar];
            // qField1[i][1][iVar] = qField[i][1][iVar];
            // qField2[i][1][iVar] = qField[i][2][iVar];

            // qField1[i][jed    ][iVar] = qField[i][jed    ][iVar];
            // qField2[i][jed    ][iVar] = qField[i][jed + 1][iVar];
            qField1[i][jed + 1][iVar] = qField[i][jed + 1][iVar];
            qField2[i][jed + 1][iVar] = qField[i][jed + 2][iVar];
        }
    }

    for (int i = ist + Iw + 1; i <= ied; i++) {
        for (int iVar = 0; iVar < num_of_primitive_vars; iVar++) {
            qField1[i][jst + Jw1][iVar] = qField[i][jst + Jw1][iVar];
            qField2[i][jst + Jw1][iVar] = qField[i][jst + Jw1 + 1][iVar];

            qField1[i][jst + Jw1 + 1][iVar] = qField[i][jst + Jw1 + 1][iVar];
            qField2[i][jst + Jw1 + 1][iVar] = qField[i][jst + Jw1 + 2][iVar];

            qField1[i][jst + Jw2 - 1][iVar] = qField[i][jst + Jw2 - 1][iVar];
            qField2[i][jst + Jw2 - 1][iVar] = qField[i][jst + Jw2][iVar];

            qField1[i][jst + Jw2 - 2][iVar] = qField[i][jst + Jw2 - 2][iVar];
            qField2[i][jst + Jw2 - 2][iVar] = qField[i][jst + Jw2 - 1][iVar];
        }
    }
}

void QlQr_Solver::QlQr_MUSCL_X() {
    // 在x方向进行插值
    VInt2D& marker = mesh->Get_Marker_F();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 1; i <= ied; i++)  // 虚拟点i=0, ied + 1的值还没有
    {
        for (int j = jst; j <= jed; j++) {
            if (marker[i][j] == 0)
                continue;

            // 需要用四个点进行插值
            VDouble qVector_m1 = qField[i - 1][j];
            VDouble qVector_c0 = qField[i][j];
            VDouble qVector_p1 = qField[i + 1][j];
            VDouble qVector_p2 = qField[i + 2][j];

            for (int iVar = 0; iVar < num_of_primitive_vars; iVar++) {
                double du_p1 = qVector_p1[iVar] - qVector_c0[iVar] + SMALL;
                double du_m1 = qVector_c0[iVar] - qVector_m1[iVar] + SMALL;
                double du_p3 = qVector_p2[iVar] - qVector_p1[iVar] + SMALL;

                double ita_m1_p = du_p1 / du_m1;
                double ita_p1_m = du_m1 / du_p1;
                double ita_p3_m = du_p1 / du_p3;
                double ita_p1_p = du_p3 / du_p1;

                double fai1 = Limiter_Function(ita_m1_p);
                double fai2 = Limiter_Function(ita_p1_m);
                double fai3 = Limiter_Function(ita_p3_m);
                double fai4 = Limiter_Function(ita_p1_p);

                qField1[i][j][iVar] =
                    qVector_c0[iVar] +
                    1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1 + (1 + muscl_k) * fai2 * du_p1);  // left

                qField2[i][j][iVar] =
                    qVector_p1[iVar] -
                    1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3 + (1 + muscl_k) * fai4 * du_p1);  // right
            }
        }
    }
}

void QlQr_Solver::QlQr_MUSCL_Y() {
    // 在y方向进行插值
    VInt2D& marker = mesh->Get_Marker_F();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = ist; i <= ied; i++) {
        for (int j = 1; j <= jed; j++)  // 虚拟点j=0, jed + 1的值还没有
        {
            if (marker[i][j] == 0)
                continue;

            VDouble qVector_m1 = qField[i][j - 1];
            VDouble qVector_c0 = qField[i][j];
            VDouble qVector_p1 = qField[i][j + 1];
            VDouble qVector_p2 = qField[i][j + 2];

            for (int iVar = 0; iVar < num_of_primitive_vars; iVar++) {
                double du_p1 = qVector_p1[iVar] - qVector_c0[iVar] + SMALL;
                double du_m1 = qVector_c0[iVar] - qVector_m1[iVar] + SMALL;
                double du_p3 = qVector_p2[iVar] - qVector_p1[iVar] + SMALL;

                double ita_m1_p = du_p1 / du_m1;
                double ita_p1_m = du_m1 / du_p1;
                double ita_p3_m = du_p1 / du_p3;
                double ita_p1_p = du_p3 / du_p1;

                double fai1 = Limiter_Function(ita_m1_p);
                double fai2 = Limiter_Function(ita_p1_m);
                double fai3 = Limiter_Function(ita_p3_m);
                double fai4 = Limiter_Function(ita_p1_p);

                qField1[i][j][iVar] =
                    qVector_c0[iVar] + 1.0 / 4.0 * ((1 - muscl_k) * fai1 * du_m1 + (1 + muscl_k) * fai2 * du_p1);

                qField2[i][j][iVar] =
                    qVector_p1[iVar] - 1.0 / 4.0 * ((1 - muscl_k) * fai3 * du_p3 + (1 + muscl_k) * fai4 * du_p1);
            }
        }
    }
}

double QlQr_Solver::Limiter_Function(double ita) {
    if (method_of_limiter == 0)  // nolim
    {
        return 1.0;
    } else if (method_of_limiter == 1)  // vanleer
    {
        return vanleer_limiter(ita, 1.0);
    } else if (method_of_limiter == 2)  // minmod
    {
        return minmod_limiter(ita, 1.0);
    } else if (method_of_limiter == 3)  // superbee
    {
        return superbee_limiter(ita, 1.0);
    }
}

// 限制器函数
double minmod_limiter(double a, double b) {
    if (a * b <= 0) {
        return 0.0;
    } else {
        if ((fabs(a) - fabs(b)) > 0) {
            return b;
        } else {
            return a;
        }
    }
}

double vanleer_limiter(double a, double) {
    return (a + fabs(a)) / (1.0 + fabs(a));
}

double superbee_limiter(double a, double) {
    double tmp1 = min(2.0 * a, 1.0);
    double tmp2 = min(a, 2.0);

    return max(0.0, max(tmp1, tmp2));
}
