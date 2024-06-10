//////////////////////////////////////////////////////////
//														//
//			Blunt Solver on Structured Grid			//
//					by Hu yiyue						//
//				Email: yiyuehuu@gmail.com				//
//					   2024.6.2							//
//////////////////////////////////////////////////////////
#pragma once
#include <Global.h>

#include <vector>
using namespace std;

extern VDouble3D fluxVector;
class Flux_Solver {
   public:
    Flux_Solver();
    ~Flux_Solver() {};

   protected:
    double gama;
    int ist, ied, jst, jed;
    VDouble3D fluxVector1;
    VDouble3D fluxVector2;

   public:
    void Solve_Flux();

   protected:
    void Roe_Scheme();
    void Roe_Scheme_X();
    void Roe_Scheme_Y();

   protected:
    void Inviscid_Flux_F(VDouble& fluxVector, double rho, double u, double v, double p);
    void Inviscid_Flux_G(VDouble& fluxVector, double rho, double u, double v, double p);

    double Enthalpy(double rho, double u, double v, double p, double gama);
    void EntropyFix(double& lamda1, double& lamda2, double& lamda3, double& lamda4);
    void EntropyFix_Harten(double& lamda);

    void Compute_Jacobian(VDouble2D& Jacobian, double u, double v, double c, double H);
    void Compute_Jacobian_X(VDouble2D& Jacobian, double u, double v, double c, double H, double lamda1, double lamda2,
                            double lamda3, double lamda4);
    void Compute_Jacobian_Y(VDouble2D& Jacobian, double u, double v, double c, double H, double lamda1, double lamda2,
                            double lamda3, double lamda4);
};
