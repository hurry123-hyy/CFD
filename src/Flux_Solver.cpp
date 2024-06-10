//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>

#include "Flux_Solver.h"
#include "Blunt_Solver.h"
#include "Geometry.h"
#include "Global.h"
#include "QlQr_Solver.h"

VDouble3D fluxVector;

void Solve_Flux() {
    Flux_Solver* fluxSolver = new Flux_Solver();

    fluxSolver->Solve_Flux();

    delete fluxSolver;
}

Flux_Solver::Flux_Solver() {
    this->gama = 1.4;

    Get_IJK_Region(ist, ied, jst, jed);

    Allocate_3D_Vector(fluxVector, num_half_point_x, num_half_point_y, num_of_primitive_vars);
    Allocate_3D_Vector(fluxVector1, total_points_x, total_points_y, num_of_primitive_vars);
    Allocate_3D_Vector(fluxVector2, total_points_x, total_points_y, num_of_primitive_vars);
}

void Flux_Solver::Solve_Flux() {
    Roe_Scheme();
}

void Flux_Solver::Roe_Scheme() {
    if (solve_direction == 'x') {
        this->Roe_Scheme_X();
    } else if (solve_direction == 'y') {
        this->Roe_Scheme_Y();
    } else {
        cout << "出错，请检查！" << endl;
    }
}

void Flux_Solver::Roe_Scheme_X() {
    VInt2D& marker = mesh->Get_Marker_F();

#ifndef _OPENMP
    VDouble fluxVector11(num_of_primitive_vars);
    VDouble fluxVector22(num_of_primitive_vars);
    VDouble2D Jacobian_A;
    Allocate_2D_Vector(Jacobian_A, num_of_primitive_vars, num_of_primitive_vars);
#endif

#ifdef _OPENMP
#pragma omp parallel
    {
        VDouble fluxVector11(num_of_primitive_vars);
        VDouble fluxVector22(num_of_primitive_vars);
        VDouble2D Jacobian_A;
        Allocate_2D_Vector(Jacobian_A, num_of_primitive_vars, num_of_primitive_vars);
#pragma omp for
#endif
        for (int i = 0; i <= ied + 1; i++) {
            for (int j = jst; j <= jed; j++) {
                if (marker[i][j] == 0)
                    continue;

                double rho1, u1, v1, p1;
                // 左值通量
                ExtractValue(qField1[i][j], rho1, u1, v1, p1);
                double H1 = Enthalpy(rho1, u1, v1, p1, gama);

                Inviscid_Flux_F(fluxVector11, rho1, u1, v1, p1);

                // 右值通量
                double rho2, u2, v2, p2;
                ExtractValue(qField2[i][j], rho2, u2, v2, p2);
                double H2 = Enthalpy(rho2, u2, v2, p2, gama);

                Inviscid_Flux_F(fluxVector22, rho2, u2, v2, p2);

                double D = sqrt(rho2 / rho1);

                double rho_roe = sqrt(rho1 * rho2);
                double u_roe = (u1 + u2 * D) / (1 + D);
                double v_roe = (v1 + v2 * D) / (1 + D);
                double H_roe = (H1 + H2 * D) / (1 + D);
                double c2_roe =
                    (gama - 1) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe));  // 这里是否存在错误？u_roe^2即可
                double c_roe = sqrt(fabs(c2_roe));                                 // 声速取绝对值

                Compute_Jacobian(Jacobian_A, u_roe, v_roe, c_roe, H_roe);

                VDouble qPrimitive1 = qField1[i][j];
                VDouble qPrimitive2 = qField2[i][j];

                VDouble qConservative1(num_of_primitive_vars);
                VDouble qConservative2(num_of_primitive_vars);
                Primitive_To_Conservative(qPrimitive1, qConservative1);
                Primitive_To_Conservative(qPrimitive2, qConservative2);

                VDouble dq(num_of_primitive_vars);
                dq[IR] = qConservative2[IR] - qConservative1[IR];
                dq[IU] = qConservative2[IU] - qConservative1[IU];
                dq[IV] = qConservative2[IV] - qConservative1[IV];
                dq[IP] = qConservative2[IP] - qConservative1[IP];

                VDouble flux_add(num_of_primitive_vars);
                MatrixMultiply(Jacobian_A, dq, flux_add, 4, 4);

                fluxVector[i][j][IR] = 0.5 * (fluxVector11[IR] + fluxVector22[IR] - flux_add[IR]);
                fluxVector[i][j][IU] = 0.5 * (fluxVector11[IU] + fluxVector22[IU] - flux_add[IU]);
                fluxVector[i][j][IV] = 0.5 * (fluxVector11[IV] + fluxVector22[IV] - flux_add[IV]);
                fluxVector[i][j][IP] = 0.5 * (fluxVector11[IP] + fluxVector22[IP] - flux_add[IP]);
            }
        }
#ifdef _OPENMP
    }
#endif  // _OPENMP
}

void Flux_Solver::Roe_Scheme_Y() {
    VInt2D& marker = mesh->Get_Marker_F();
#ifndef _OPENMP
    VDouble2D Jacobian_A;
    Allocate_2D_Vector(Jacobian_A, num_of_primitive_vars, num_of_primitive_vars);
    VDouble fluxVector11(num_of_primitive_vars);
    VDouble fluxVector22(num_of_primitive_vars);
#endif

#ifdef _OPENMP
#pragma omp parallel
    {
        VDouble2D Jacobian_A;
        Allocate_2D_Vector(Jacobian_A, num_of_primitive_vars, num_of_primitive_vars);
        VDouble fluxVector11(num_of_primitive_vars);
        VDouble fluxVector22(num_of_primitive_vars);
#pragma omp for
#endif
        for (int i = ist; i <= ied; i++) {
            for (int j = 0; j <= jed + 1; j++) {
                if (marker[i][j] == 0)
                    continue;

                double rho1, u1, v1, p1;
                ExtractValue(qField1[i][j], rho1, u1, v1, p1);
                double H1 = Enthalpy(rho1, u1, v1, p1, gama);

                Inviscid_Flux_G(fluxVector11, rho1, u1, v1, p1);

                double rho2, u2, v2, p2;
                ExtractValue(qField2[i][j], rho2, u2, v2, p2);
                double H2 = Enthalpy(rho2, u2, v2, p2, gama);

                Inviscid_Flux_G(fluxVector22, rho2, u2, v2, p2);

                double D = sqrt(rho2 / rho1);

                double rho_roe = sqrt(rho1 * rho2);
                double u_roe = (u1 + u2 * D) / (1 + D);
                double v_roe = (v1 + v2 * D) / (1 + D);
                double H_roe = (H1 + H2 * D) / (1 + D);
                double c2_roe = (gama - 1) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe));
                double c_roe = sqrt(fabs(c2_roe));  // 声速取绝对值

                Compute_Jacobian(Jacobian_A, u_roe, v_roe, c_roe, H_roe);

                VDouble qPrimitive1 = qField1[i][j];
                VDouble qPrimitive2 = qField2[i][j];

                VDouble qConservative1(num_of_primitive_vars);
                VDouble qConservative2(num_of_primitive_vars);
                Primitive_To_Conservative(qPrimitive1, qConservative1);
                Primitive_To_Conservative(qPrimitive2, qConservative2);

                VDouble dq(num_of_primitive_vars);
                dq[IR] = qConservative2[IR] - qConservative1[IR];
                dq[IU] = qConservative2[IU] - qConservative1[IU];
                dq[IV] = qConservative2[IV] - qConservative1[IV];
                dq[IP] = qConservative2[IP] - qConservative1[IP];

                VDouble flux_add(num_of_primitive_vars);
                MatrixMultiply(Jacobian_A, dq, flux_add, 4, 4);

                fluxVector[i][j][IR] = 0.5 * (fluxVector11[IR] + fluxVector22[IR] - flux_add[IR]);
                fluxVector[i][j][IU] = 0.5 * (fluxVector11[IU] + fluxVector22[IU] - flux_add[IU]);
                fluxVector[i][j][IV] = 0.5 * (fluxVector11[IV] + fluxVector22[IV] - flux_add[IV]);
                fluxVector[i][j][IP] = 0.5 * (fluxVector11[IP] + fluxVector22[IP] - flux_add[IP]);
            }
        }
#ifdef _OPENMP
    }
#endif  // _OPENMP
}

void Flux_Solver::Inviscid_Flux_F(VDouble& fluxVector, double rho, double u, double v, double p) {
    double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
    fluxVector[IR] = rho * u;
    fluxVector[IU] = rho * u * u + p;
    fluxVector[IV] = rho * u * v;
    fluxVector[IP] = (E + p) * u;
}

void Flux_Solver::Inviscid_Flux_G(VDouble& fluxVector, double rho, double u, double v, double p) {
    double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
    fluxVector[IR] = rho * v;
    fluxVector[IU] = rho * u * v;
    fluxVector[IV] = rho * v * v + p;
    fluxVector[IP] = (E + p) * v;
}

double Flux_Solver::Enthalpy(double rho, double u, double v, double p, double gama) {
    double E = p / (gama - 1) + 0.5 * rho * (u * u + v * v);
    return (E + p) / rho;
}

void Flux_Solver::EntropyFix(double& lamda1, double& lamda2, double& lamda3, double& lamda4) {
    this->EntropyFix_Harten(lamda1);
    this->EntropyFix_Harten(lamda2);
    this->EntropyFix_Harten(lamda3);
    this->EntropyFix_Harten(lamda4);
}

void Flux_Solver::EntropyFix_Harten(double& lamda) {
    double eps = entropy_fix_coeff;  // 熵格式修正系数，0.1或者0.01
    lamda = fabs(lamda) > eps ? fabs(lamda) : ((lamda * lamda + eps * eps) / 2.0 / eps);
}

void Flux_Solver::Compute_Jacobian(VDouble2D& Jacobian, double u, double v, double a, double h) {
    // a为声速
    if (solve_direction == 'x') {
        double lamda1 = u;
        double lamda2 = u;
        double lamda3 = u - a;
        double lamda4 = u + a;

        EntropyFix(lamda1, lamda2, lamda3, lamda4);

        Compute_Jacobian_X(Jacobian, u, v, a, h, lamda1, lamda2, lamda3, lamda4);
    } else if (solve_direction == 'y') {
        double lamda1 = v;
        double lamda2 = v;
        double lamda3 = v - a;
        double lamda4 = v + a;

        EntropyFix(lamda1, lamda2, lamda3, lamda4);

        Compute_Jacobian_Y(Jacobian, u, v, a, h, lamda1, lamda2, lamda3, lamda4);
    }
}

void Flux_Solver::Compute_Jacobian_X(VDouble2D& Jacobian, double u, double v, double a, double h, double lamda1,
                                     double lamda2, double lamda3, double lamda4) {
    // 特征值对角矩阵
    VDouble2D Lmdx = {{lamda1, 0, 0, 0}, {0, lamda2, 0, 0}, {0, 0, lamda3, 0}, {0, 0, 0, lamda4}};

    // 右特征矩阵
    VDouble2D Rx = {{1, 0, 1, 1}, {u, 0, u - a, u + a}, {0, 1, v, v}, {(u * u - v * v) / 2, v, h - a * u, h + a * u}};

    double b2 = (gama - 1) / a / a;
    double b1 = b2 * (u * u + v * v) / 2.0;

    // 左特征矩阵
    VDouble2D Lx = {{1 - b1, b2 * u, b2 * v, -b2},
                    {-b1 * v, b2 * u * v, 1 + b2 * v * v, -b2 * v},
                    {(b1 + u / a) / 2, -(b2 * u + 1 / a) / 2, -b2 * v / 2, b2 / 2},
                    {(b1 - u / a) / 2, -(b2 * u - 1 / a) / 2, -b2 * v / 2, b2 / 2}};

    VDouble2D tmp1;
    Allocate_2D_Vector(tmp1, 4, 4);

    MatrixMultiply(Rx, Lmdx, tmp1, 4, 4, 4);
    MatrixMultiply(tmp1, Lx, Jacobian, 4, 4, 4);
}

void Flux_Solver::Compute_Jacobian_Y(VDouble2D& Jacobian, double u, double v, double a, double h, double lamda1,
                                     double lamda2, double lamda3, double lamda4) {
    VDouble2D Lmdy = {{lamda1, 0, 0, 0}, {0, lamda2, 0, 0}, {0, 0, lamda3, 0}, {0, 0, 0, lamda4}};

    VDouble2D Ry = {{0, 1, 1, 1}, {1, 0, u, u}, {0, v, v - a, v + a}, {u, (v * v - u * u) / 2, h - a * v, h + a * v}};

    double b2 = (gama - 1) / a / a;
    double b1 = b2 * (u * u + v * v) / 2.0;

    VDouble2D Ly = {{-b1 * u, 1 + b2 * u * u, b2 * u * v, -b2 * u},
                    {1 - b1, b2 * u, b2 * v, -b2},
                    {(b1 + v / a) / 2, -b2 * u / 2, -(b2 * v + 1 / a) / 2, b2 / 2},
                    {(b1 - v / a) / 2, -b2 * u / 2, -(b2 * v - 1 / a) / 2, b2 / 2}};

    VDouble2D tmp1;
    Allocate_2D_Vector(tmp1, 4, 4);

    MatrixMultiply(Ry, Lmdy, tmp1, 4, 4, 4);
    MatrixMultiply(tmp1, Ly, Jacobian, 4, 4, 4);
}
