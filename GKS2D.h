#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

class Cell2d
{
public:
	int loc;
	int num_gauss;
	int stage;
	double dx;
	double dt;
};

class Block2d
{
public:
    int N;
    int step = 0;
    int stages;
    int ghost;
    int gausspoint;
    string variable_type;
    double x;
    double y;
    double dx;
    double dt;
    double CFL;
    double tnow;
    double tstop;
    double gauss_loc[4]{};
    double gauss_weight[4]{};
    double timecoefficient[5][5][3]{};
};

class Fluid2d
{
public:
	double W[4]{};
	double W_old[4]{};
	double prim[4]{0, 0, 0, 0};
	double WL[4]{};
	double WR[4]{};
	double WU[4]{};
	double WD[4]{};
	double Der1Lx[4]{};  // refers the left side normal first derivative
	double Der1Rx[4]{};
	double Der1Ux[4]{};
	double Der1Dx[4]{};
	/*
	* timecoefficient[i][j][k]
	* i refers the i stage
	* j refers the nth coefficient at n stage
	* k refers f, derf, der2f
	*/
	double GaussWL[4][4]{};
	double GaussWR[4][4]{};
	double GaussWD[4][4]{};
	double GaussWU[4][4]{};
	double GaussDer1Lx[4][4]{};
	double GaussDer1Ly[4][4]{};
	double GaussDer1Rx[4][4]{};
	double GaussDer1Ry[4][4]{};
	double GaussDer1Dx[4][4]{};
	double GaussDer1Dy[4][4]{};
	double GaussDer1Ux[4][4]{};
	double GaussDer1Uy[4][4]{};
	double alpha = 1;
};

class Interface2d
{
public:
	// Gauss[n][m]
	// n refers the nth gauss point
	// m refers the variables rho, rhoU, rhoV, rhoE
	double GaussWL[4][4]{};
	double GaussWR[4][4]{};
	double GaussWcenter[4][4]{};
	double GaussDer1Lx[4][4]{};
	double GaussDer1Ly[4][4]{};
	double GaussDer1Rx[4][4]{};
	double GaussDer1Ry[4][4]{};
	double GaussCenterDer1x[4][4]{};
	double GaussCenterDer1y[4][4]{};
	// flux[i][j][k]
	// i refers the ith stage, j refers the jth gauss point
	// k refers the variabels
	double flux[4][4][4]{};
	double der_flux[4][4][4]{};
	double Flux[4][4][4]{};
	double c1 = 0.05;
	double c2 = 1;
};

class MMDF
{
public:
	double uwhole[7]{};
	double vwhole[7]{};
	double uplus[7]{};
	double uminus[7]{};
	double vplus[7]{};
	double vminus[7]{};
	// u, half domain > 0; v, whole domain; xi, equals 0/2/4
	double upvxi[7][7][3]{};
	// u, half domain < 0; v, whole domain; xi, equals 0/2/4
	double unvxi[7][7][3]{};
	// u, whole domain; v, whole domain; xi, equals 0/2/4
	double uvxi[7][7][3]{};
	double overlambda;
	double xi2;
	double xi4;
};

void Riemann_problem_2d();

void Accuracy_test_2d();
