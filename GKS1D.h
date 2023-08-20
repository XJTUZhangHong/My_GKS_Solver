#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

class Block1d
{
public:
    int ghost;
    int nodex;
    int nx;
    double dx;
    int stages;
    double timecoefficient[5][5][2]{};
    double dt;
    double tstop;
    double tnow;
    double CFL;
    int step = 0;
    double c1;
    double c2;
    string variable_type;
};

class Fluid
{
    /*
     * W = [rho, rhoU, rhoE]
     * WL = [rho_left, rhoU_left, rhoE_left]
     * WR = [rho_right, rhoU_right, rhoE_right]
     * DerL = [der_rho_left, der_rhoU_left, der_rhoE_left]
     * DerR = [der_rho_right, der_rhoU_right, der_rhoE_right]
     * prim = [rho, U, p]
     * time_coefficient[i][j][k]:
     * i -- refers the n stage
     * j -- refers the n-th coefficient at n stage
     * k -- refers flux, der1flux, der2flux
     * Note: the block left side correspond to the interface right side,
     * the block right side as well
     * */
public:
    double step = 0;
    double W[3]{};
    double W_old[3]{};
    double WL[3]{};
    double WR[3]{};
    double DerL[3]{};
    double DerR[3]{};
    double prim[3]{};
    double alpha = 1.0;
};

class Interface
{
public:
    double W[3]{};
    double WL[3]{};
    double WR[3]{};
    double WC[3]{};
    double DerL[3]{};
    double DerR[3]{};
    double DerC[3]{};
    double al[3]{};
    double ar[3]{};
    double ac[3]{};
    double AL[3]{};
    double AR[3]{};
    double AC[3]{};
    double Fluxes[2][3]{};
    double fluxes[2][3]{}; // flux[stage][rho_flux, rhoU_flux, rhoE_flux]
    double der1fluxes[2][3]{};
    double tau = 0;
};

class MMDF1d
{
public:
    double uwhole[10]{};
    double uplus[10]{};
    double uminus[10]{};
    double upxi[10][4]{};
    double unxi[10][4]{};
    double uxi[10][4]{};
    double xi2 = 0;
    double xi4 = 0;
};

void Riemann_problem_1d();

void Accuracy_test_1d();
