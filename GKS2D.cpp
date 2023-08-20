/*
 * Two dimension second order gas-kinetic scheme
 * Using uniform(structure) mesh
 * Multi-stages multi-derivatives
 * S1O1: one stage first derivatives
 * S1O2: one stage second derivatives
 * S2O4: two stages forth derivatives
 * Vanleer reconstruction: second order spatial accuracy
 * WENO5-AO reconstruction: fifth order spatial accuracy
 * Now only applied for Euler problems
 * Reference: https://www.osredm.com/p35462178/gks2d-str
 * */

#include "GKS2D.h"
#define Gamma 1.4
#define tol 1e-10
#define pi acos(-1)
#define K 3

double Alpha(double lambda, double u)
{
    return erfc(sqrt(lambda) * u);
}

double Beta(double lambda, double u)
{
    return exp(-lambda * u * u) / sqrt(pi * lambda);
}

void Initial_variables(Fluid2d* fluid, double zone1[], double zone2[], double zone3[], double zone4[], Block2d block)
{
    double dx = block.dx;
    int ghost = block.ghost;
    double x = block.x, y = block.y;
    int n = block.N + 2 * block.ghost;
#pragma omp parallel for
    for (int i = ghost; i < n - ghost; i++)
    {
        for (int j = ghost; j < n - ghost; j++)
        {
            if ((double)(i - ghost + 1) * dx <= x && (double)(j - ghost + 1) * dx <= y)
            {
                for (int k = 0; k < 4; k++)
                {
                    fluid[i * n + j].prim[k] = zone1[k];
                }
            }
            else if ((double)(i - ghost + 1) * dx > x && (double)(j - ghost + 1) * dx <= y)
            {
                for (int k = 0; k < 4; k++)
                {
                    fluid[i * n + j].prim[k] = zone2[k];
                }
            }
            else if ((double)(i - ghost + 1) * dx > x && (double)(j - ghost + 1) * dx > y)
            {
                for (int k = 0; k < 4; k++)
                {
                    fluid[i * n + j].prim[k] = zone3[k];
                }
            }
            else
            {
                for (int k = 0; k < 4; k++)
                {
                    fluid[i * n + j].prim[k] = zone4[k];
                }
            }
            fluid[i * n + j].W[0] = fluid[i * n + j].prim[0];
            fluid[i * n + j].W[1] = fluid[i * n + j].prim[1] * fluid[i * n + j].prim[0];
            fluid[i * n + j].W[2] = fluid[i * n + j].prim[2] * fluid[i * n + j].prim[0];
            fluid[i * n + j].W[3] = 2.5 * fluid[i * n + j].prim[3] +
                0.5 * fluid[i * n + j].prim[0] * (pow(fluid[i * n + j].prim[1], 2) + pow(fluid[i * n + j].prim[2], 2));
        }
    }

}

void SetGaussPoint(Block2d& block)
{
    double gauss = block.gausspoint;
#pragma omp parallel for
    if (gauss == 1)
    {
        block.gauss_loc[0] = 0.0;
        block.gauss_weight[0] = 1.0;
    }
    if (gauss == 2)
    {
        block.gauss_loc[0] = -sqrt(1.0 / 3.0);
        block.gauss_loc[0] = sqrt(1.0 / 3.0);
        block.gauss_weight[0] = 0.5;
        block.gauss_weight[1] = 0.5;
    }
}

void Copy_array(Fluid2d* fluid, Block2d block)
{
    int n = pow(block.N + 2 * block.ghost, 2);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            fluid[i].W_old[j] = fluid[i].W[j];
        }
    }
}

void YchangetoX(double* w, double* convar)
{
    w[0] = convar[0];
    w[1] = convar[2];
    w[2] = -convar[1];
    w[3] = convar[3];
}

double Lambda(double R, double RU, double RV, double RE)
{
    return(K + 2.0) * 0.25 * (R / (RE - 0.5 * (RU * RU + RV * RV) / R));
}

void Convar_to_ULambda_2d(double* prim, double* convar)
{
    prim[0] = convar[0];
    prim[1] = convar[1] / prim[0];
    prim[2] = convar[2] / prim[0];
    // prim[3] refers lambda
    prim[3] = Lambda(convar[0], convar[1], convar[2], convar[3]);
}

void VanLeer_2D_normal(int i, int j, Fluid2d* fluid, Block2d block)
{
    int n = block.N + 2 * block.ghost, ghost = block.ghost, N = block.N;
    // x direction normal derivatives
    if (i > ghost - 2 && i < N + ghost + 1)
    {
        double splus[4]{}, sminus[4]{}, w[4]{}, wp[4]{}, wn[4]{};
        for (int k = 0; k < 4; k++)
        {
            w[k] = fluid[i * n + j].W[k];
            wp[k] = fluid[(i + 1) * n + j].W[k];
            wn[k] = fluid[(i - 1) * n + j].W[k];
            splus[k] = (wp[k] - w[k]) / block.dx;
            sminus[k] = (w[k] - wn[k]) / block.dx;
        }
        for (int k = 0; k < 4; k++)
        {
            if ((splus[k] * sminus[k]) > 0)
            {
                fluid[i * n + j].Der1Lx[k] = 2 * splus[k] * sminus[k] / (splus[k] + sminus[k]);
                fluid[i * n + j].Der1Rx[k] = fluid[i * n + j].Der1Lx[k];
            }
            else
            {
                fluid[i * n + j].Der1Lx[k] = 0.0;
                fluid[i * n + j].Der1Rx[k] = 0.0;
            }
            fluid[i * n + j].WL[k] = w[k] - 0.5 * block.dx * fluid[i * n + j].Der1Lx[k];
            fluid[i * n + j].WR[k] = w[k] + 0.5 * block.dx * fluid[i * n + j].Der1Rx[k];
        }
        // if lambda < 0, then reduce to the first order
        double flag_l = fluid[i * n + j].WL[3] - 0.5 * (pow(fluid[i * n + j].WL[1], 2) + pow(fluid[i * n + j].WL[2], 2)) / fluid[i * n + j].WL[0];
        double flag_r = fluid[i * n + j].WR[3] - 0.5 * (pow(fluid[i * n + j].WR[1], 2) + pow(fluid[i * n + j].WR[2], 2)) / fluid[i * n + j].WR[0];
        if (flag_l <= 0 || flag_r <= 0)
        {
            for (int k = 0; k <= 3; k++)
            {
                fluid[i * n + j].WL[k] = w[k];
                fluid[i * n + j].WR[k] = w[k];
                fluid[i * n + j].Der1Lx[k] = 0.0;
                fluid[i * n + j].Der1Rx[k] = 0.0;
            }
        }
    }
    // y direction normal derivatives
    if (j > ghost - 2 && j < N + ghost + 1)
    {
        double splus[4]{}, sminus[4]{}, w[4]{}, wp[4]{}, wn[4]{};
        double convar[4]{}, convarp[4]{}, convarn[4]{};
        for (int k = 0; k < 4; k++)
        {
            convar[k] = fluid[i * n + j].W[k];
            convarp[k] = fluid[i * n + j + 1].W[k];
            convarn[k] = fluid[i * n + j - 1].W[k];
        }
        YchangetoX(w, convar); YchangetoX(wp, convarp); YchangetoX(wn, convarn);
        for (int k = 0; k < 4; k++)
        {
            splus[k] = (wp[k] - w[k]) / block.dx;
            sminus[k] = (w[k] - wn[k]) / block.dx;
        }
        for (int k = 0; k < 4; k++)
        {
            if ((splus[k] * sminus[k]) > 0)
            {
                fluid[i * n + j].Der1Dx[k] = 2 * splus[k] * sminus[k] / (splus[k] + sminus[k]);
                fluid[i * n + j].Der1Ux[k] = fluid[i * n + j].Der1Dx[k];
            }
            else
            {
                fluid[i * n + j].Der1Dx[k] = 0.0;
                fluid[i * n + j].Der1Ux[k] = 0.0;
            }
            fluid[i * n + j].WD[k] = w[k] - 0.5 * block.dx * fluid[i * n + j].Der1Dx[k];
            fluid[i * n + j].WU[k] = w[k] + 0.5 * block.dx * fluid[i * n + j].Der1Ux[k];
        }
        // if lambda < 0, then reduce to the first order
        double flag_d = fluid[i * n + j].WD[3] - 0.5 * (pow(fluid[i * n + j].WD[1], 2) + pow(fluid[i * n + j].WD[2], 2)) / fluid[i * n + j].WD[0];
        double flag_u = fluid[i * n + j].WU[3] - 0.5 * (pow(fluid[i * n + j].WU[1], 2) + pow(fluid[i * n + j].WU[2], 2)) / fluid[i * n + j].WU[0];
        if (flag_d <= 0 || flag_u <= 0)
        {
            for (int k = 0; k <= 3; k++)
            {
                fluid[i * n + j].WD[k] = w[k];
                fluid[i * n + j].WU[k] = w[k];
                fluid[i * n + j].Der1Dx[k] = 0.0;
                fluid[i * n + j].Der1Ux[k] = 0.0;
            }
        }
    }
}

void VanLeer_2D_tangential(double* coe, double wn, double w, double wp, double dx)
{
    if ((w - wn) * (wp - w) > 0)
    {
        double slope_left = (w - wn) / dx;
        double slope_right = (wp - w) / dx;
        coe[1] = 2 * slope_left * slope_right / (slope_left + slope_right);
        coe[0] = w;
    }
    else
    {
        coe[1] = 0.0;
        coe[0] = w;
    }
}

void VanLeer_2D_tangent(int i, int j, Fluid2d* fluid, Block2d block)
{
    int n = block.N + 2 * block.ghost;
    int gausspoint = block.gausspoint;
    double rho, u, v, rhoE;
    double w[4]{}, wn[4]{}, wp[4]{}, Derwn[4]{}, Derw[4]{}, Derwp[4]{};
    // x direction
    for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
    {
        // cell left side
        for (int var = 0; var < 4; var++)
        {
            double coe[2];
            wn[var] = fluid[i * n + j - 1].WL[var];
            w[var] = fluid[i * n + j].WL[var];
            wp[var] = fluid[i * n + j + 1].WL[var];
            Derwn[var] = fluid[i * n + j - 1].Der1Lx[var];
            Derw[var] = fluid[i * n + j].Der1Lx[var];
            Derwp[var] = fluid[i * n + j + 1].Der1Lx[var];
            VanLeer_2D_tangential(coe, wn[var], w[var], wp[var], block.dx);
            fluid[i * n + j].GaussWL[num_gauss][var] = coe[0] + 0.5 * block.gauss_loc[num_gauss] * block.dx * coe[1];
            fluid[i * n + j].GaussDer1Ly[num_gauss][var] = coe[1];
            VanLeer_2D_tangential(coe, Derwn[var], Derw[var], Derwp[var], block.dx);
            fluid[i * n + j].GaussDer1Lx[num_gauss][var] = coe[0] + 0.5 * block.dx * block.gauss_loc[num_gauss] * coe[1];
        }
        // check order reduce
        rho = fluid[i * n + j].GaussWL[num_gauss][0];
        u = fluid[i * n + j].GaussWL[num_gauss][1] / rho;
        v = fluid[i * n + j].GaussWL[num_gauss][2] / rho;
        rhoE = fluid[i * n + j].GaussWL[num_gauss][3];
        if (rhoE - 0.5 * rho * (u * u + v * v) <= 0)
        {
            for (int m = 0; m < 4; ++m)
            {
                fluid[i * n + j].GaussWL[num_gauss][m] = fluid[i * n + j].WL[m];
                fluid[i * n + j].GaussDer1Lx[num_gauss][m] = fluid[i * n + j].Der1Lx[m];
                fluid[i * n + j].GaussDer1Ly[num_gauss][m] = 0.0;
            }
        }
        // cell right side
        for (int var = 0; var < 4; var++)
        {
            double coe[2];
            wn[var] = fluid[i * n + j - 1].WR[var];
            w[var] = fluid[i * n + j].WR[var];
            wp[var] = fluid[i * n + j + 1].WR[var];
            Derwn[var] = fluid[i * n + j - 1].Der1Rx[var];
            Derw[var] = fluid[i * n + j].Der1Rx[var];
            Derwp[var] = fluid[i * n + j + 1].Der1Rx[var];
            VanLeer_2D_tangential(coe, wn[var], w[var], wp[var], block.dx);
            fluid[i * n + j].GaussWR[num_gauss][var] = coe[0] + 0.5 * block.gauss_loc[num_gauss] * block.dx * coe[1];
            fluid[i * n + j].GaussDer1Ry[num_gauss][var] = coe[1];
            VanLeer_2D_tangential(coe, Derwn[var], Derw[var], Derwp[var], block.dx);
            fluid[i * n + j].GaussDer1Rx[num_gauss][var] = coe[0] + 0.5 * block.dx * block.gauss_loc[num_gauss] * coe[1];
        }
        // check order reduce
        rho = fluid[i * n + j].GaussWR[num_gauss][0];
        u = fluid[i * n + j].GaussWR[num_gauss][1] / rho;
        v = fluid[i * n + j].GaussWR[num_gauss][2] / rho;
        rhoE = fluid[i * n + j].GaussWR[num_gauss][3];
        if (rhoE - 0.5 * rho * (u * u + v * v) <= 0)
        {
            for (int m = 0; m < 4; ++m)
            {
                fluid[i * n + j].GaussWR[num_gauss][m] = fluid[i * n + j].WR[m];
                fluid[i * n + j].GaussDer1Rx[num_gauss][m] = fluid[i * n + j].Der1Rx[m];
                fluid[i * n + j].GaussDer1Ry[num_gauss][m] = 0.0;
            }
        }
    }
    // y direction
    for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
    {
        // cell down side
        for (int var = 0; var < 4; var++)
        {
            double coe[2];
            wn[var] = fluid[(i + 1) * n + j].WD[var];
            w[var] = fluid[i * n + j].WD[var];
            wp[var] = fluid[(i - 1) * n + j].WD[var];
            Derwn[var] = fluid[(i + 1) * n + j].Der1Dx[var];
            Derw[var] = fluid[i * n + j].Der1Dx[var];
            Derwp[var] = fluid[(i - 1) * n + j].Der1Dx[var];
            VanLeer_2D_tangential(coe, wn[var], w[var], wp[var], block.dx);
            fluid[i * n + j].GaussWD[num_gauss][var] = coe[0] + 0.5 * block.gauss_loc[num_gauss] * block.dx * coe[1];
            fluid[i * n + j].GaussDer1Dy[num_gauss][var] = coe[1];
            VanLeer_2D_tangential(coe, Derwn[var], Derw[var], Derwp[var], block.dx);
            fluid[i * n + j].GaussDer1Dx[num_gauss][var] = coe[0] + 0.5 * block.dx * block.gauss_loc[num_gauss] * coe[1];
        }
        // check order reduce
        rho = fluid[i * n + j].GaussWD[num_gauss][0];
        u = fluid[i * n + j].GaussWD[num_gauss][1] / rho;
        v = fluid[i * n + j].GaussWD[num_gauss][2] / rho;
        rhoE = fluid[i * n + j].GaussWD[num_gauss][3];
        if (rhoE - 0.5 * rho * (u * u + v * v) <= 0)
        {
            for (int m = 0; m < 4; ++m)
            {
                fluid[i * n + j].GaussWD[num_gauss][m] = fluid[i * n + j].WD[m];
                fluid[i * n + j].GaussDer1Dx[num_gauss][m] = fluid[i * n + j].Der1Dx[m];
                fluid[i * n + j].GaussDer1Dy[num_gauss][m] = 0.0;
            }
        }
        // cell up side
        for (int var = 0; var < 4; var++)
        {
            double coe[2];
            wn[var] = fluid[(i + 1) * n + j].WU[var];
            w[var] = fluid[i * n + j].WU[var];
            wp[var] = fluid[(i - 1) * n + j].WU[var];
            Derwn[var] = fluid[(i + 1) * n + j].Der1Ux[var];
            Derw[var] = fluid[i * n + j].Der1Ux[var];
            Derwp[var] = fluid[(i - 1) * n + j].Der1Ux[var];
            VanLeer_2D_tangential(coe, wn[var], w[var], wp[var], block.dx);
            fluid[i * n + j].GaussWU[num_gauss][var] = coe[0] + 0.5 * block.gauss_loc[num_gauss] * block.dx * coe[1];
            fluid[i * n + j].GaussDer1Uy[num_gauss][var] = coe[1];
            VanLeer_2D_tangential(coe, Derwn[var], Derw[var], Derwp[var], block.dx);
            fluid[i * n + j].GaussDer1Ux[num_gauss][var] = coe[0] + 0.5 * block.dx * block.gauss_loc[num_gauss] * coe[1];
        }
        // check order reduce
        rho = fluid[i * n + j].GaussWU[num_gauss][0];
        u = fluid[i * n + j].GaussWU[num_gauss][1] / rho;
        v = fluid[i * n + j].GaussWU[num_gauss][2] / rho;
        rhoE = fluid[i * n + j].GaussWU[num_gauss][3];
        if (rhoE - 0.5 * rho * (u * u + v * v) <= 0)
        {
            for (int m = 0; m < 4; ++m)
            {
                fluid[i * n + j].GaussWU[num_gauss][m] = fluid[i * n + j].WU[m];
                fluid[i * n + j].GaussDer1Ux[num_gauss][m] = fluid[i * n + j].Der1Ux[m];
                fluid[i * n + j].GaussDer1Uy[num_gauss][m] = 0.0;
            }
        }
    }
}

void VanLeer(Interface2d* xinterface, Interface2d* yinterface, Fluid2d* fluid, Block2d block)
{
    int n = block.N + 2 * block.ghost, ghost = block.ghost;
#pragma omp parallel  for
    // normal direction
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            VanLeer_2D_normal(i, j, fluid, block);
        }
    }
#pragma omp parallel  for
    // tangential direction
    for (int i = ghost - 1; i < n - ghost + 1; i++)
    {
        for (int j = ghost - 1; j < n - ghost + 1; j++)
        {
            VanLeer_2D_tangent(i, j, fluid, block);
        }
    }
    // interface reconstruction
    int gausspoint = block.gausspoint;
#pragma omp parallel  for
    for (int i = ghost; i < n - ghost + 1; i++)
    {
        for (int j = ghost; j < n - ghost + 1; j++)
        {
            for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
            {
                for (int var = 0; var < 4; var++)
                {
                    xinterface[i * n + j].GaussWL[num_gauss][var] = fluid[(i - 1) * n + j].GaussWR[num_gauss][var];
                    xinterface[i * n + j].GaussDer1Lx[num_gauss][var] = fluid[(i - 1) * n + j].GaussDer1Rx[num_gauss][var];
                    xinterface[i * n + j].GaussDer1Ly[num_gauss][var] = fluid[(i - 1) * n + j].GaussDer1Ry[num_gauss][var];
                    xinterface[i * n + j].GaussWR[num_gauss][var] = fluid[i * n + j].GaussWL[num_gauss][var];
                    xinterface[i * n + j].GaussDer1Rx[num_gauss][var] = fluid[i * n + j].GaussDer1Lx[num_gauss][var];
                    xinterface[i * n + j].GaussDer1Ry[num_gauss][var] = fluid[i * n + j].GaussDer1Ly[num_gauss][var];
                    yinterface[i * n + j].GaussWL[num_gauss][var] = fluid[i * n + j - 1].GaussWU[num_gauss][var];
                    yinterface[i * n + j].GaussDer1Lx[num_gauss][var] = fluid[i * n + j - 1].GaussDer1Ux[num_gauss][var];
                    yinterface[i * n + j].GaussDer1Ly[num_gauss][var] = fluid[i * n + j - 1].GaussDer1Uy[num_gauss][var];
                    yinterface[i * n + j].GaussWR[num_gauss][var] = fluid[i * n + j].GaussWD[num_gauss][var];
                    yinterface[i * n + j].GaussDer1Rx[num_gauss][var] = fluid[i * n + j].GaussDer1Dx[num_gauss][var];
                    yinterface[i * n + j].GaussDer1Ry[num_gauss][var] = fluid[i * n + j].GaussDer1Dy[num_gauss][var];
                }
            }
        }
    }
}

void MMDF_calculate(MMDF* m, double* prim)
{
    double u = prim[1], v = prim[2], lambda = prim[3];
    double overlambda = 1.0 / lambda;
    m[0].uwhole[0] = 1;
    m[0].uwhole[1] = u;
    m[0].vwhole[0] = 1;
    m[0].vwhole[1] = v;
    m[0].uplus[0] = 0.5 * Alpha(lambda, -u);
    m[0].uminus[0] = 1.0 - m[0].uplus[0];
    m[0].uplus[1] = u * m[0].uplus[0] + 0.5 * Beta(lambda, u);
    m[0].uminus[1] = u - m[0].uplus[1];
    m[0].xi2 = 0.5 * K * overlambda;
    m[0].xi4 = 0.25 * (K * K + 2 * K) * overlambda * overlambda;
    for (int i = 2; i <= 6; i++)
    {
        m[0].uwhole[i] = u * m[0].uwhole[i - 1] + 0.5 * (i - 1) * overlambda * m[0].uwhole[i - 2];
        m[0].vwhole[i] = v * m[0].vwhole[i - 1] + 0.5 * (i - 1) * overlambda * m[0].vwhole[i - 2];
        m[0].uplus[i] = u * m[0].uplus[i - 1] + 0.5 * (i - 1) * overlambda * m[0].uplus[i - 2];
        m[0].uminus[i] = m[0].uwhole[i] - m[0].uplus[i];
    }
    // i refers power of u, j refers power of v, and 2k power of xi
#pragma omp parallel for
    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if ((i + j + 2 * k) <= 6)
                {
                    if (k == 0)  // refers power of xi equals to 0
                    {
                        m[0].upvxi[i][j][k] = m[0].uplus[i] * m[0].vwhole[j];
                        m[0].unvxi[i][j][k] = m[0].uminus[i] * m[0].vwhole[j];
                        m[0].uvxi[i][j][k] = m[0].uwhole[i] * m[0].vwhole[j];
                    }
                    if (k == 1) // refers power of xi equals to 2
                    {
                        m[0].upvxi[i][j][k] = m[0].uplus[i] * m[0].vwhole[j] * m[0].xi2;
                        m[0].unvxi[i][j][k] = m[0].uminus[i] * m[0].vwhole[j] * m[0].xi2;
                        m[0].uvxi[i][j][k] = m[0].uwhole[i] * m[0].vwhole[j] * m[0].xi2;
                    }
                    if (k == 2) // refers power of xi equals to 4
                    {
                        m[0].upvxi[i][j][k] = m[0].uplus[i] * m[0].vwhole[j] * m[0].xi4;
                        m[0].unvxi[i][j][k] = m[0].uminus[i] * m[0].vwhole[j] * m[0].xi4;
                        m[0].uvxi[i][j][k] = m[0].uwhole[i] * m[0].vwhole[j] * m[0].xi4;
                    }
                }
            }
        }
    }
}

void A_point(double* a, double* der, double* prim)
{
    double R4, R3, R2;
    double overden = 1.0 / prim[0];
    R4 = der[3] * overden - 0.5 * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]) * der[0] * overden;
    R3 = (der[2] - prim[2] * der[0]) * overden;
    R2 = (der[1] - prim[1] * der[0]) * overden;
    a[3] = (4.0 / (K + 2)) * prim[3] * prim[3] * (R4 - 2 * prim[1] * R2 - 2 * prim[2] * R3);
    a[2] = 2 * prim[3] * R3 - prim[2] * a[3];
    a[1] = 2 * prim[3] * R2 - prim[1] * a[3];
    a[0] = der[0] * overden - prim[1] * a[1] - prim[2] * a[2] - 0.5 * a[3] * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]);
}

void GL_address(int no_u, int no_v, int no_xi, double* psi, double* a, MMDF* m)
{
    psi[0] = a[0] * m[0].upvxi[no_u][no_v][no_xi] + a[1] * m[0].upvxi[no_u + 1][no_v][no_xi] + a[2] * m[0].upvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m[0].upvxi[no_u + 2][no_v][no_xi] + m[0].upvxi[no_u][no_v + 2][no_xi] + m[0].upvxi[no_u][no_v][no_xi + 1]);
    psi[1] = a[0] * m[0].upvxi[no_u + 1][no_v][no_xi] + a[1] * m[0].upvxi[no_u + 2][no_v][no_xi] + a[2] * m[0].upvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m[0].upvxi[no_u + 3][no_v][no_xi] + m[0].upvxi[no_u + 1][no_v + 2][no_xi] + m[0].upvxi[no_u + 1][no_v][no_xi + 1]);
    psi[2] = a[0] * m[0].upvxi[no_u][no_v + 1][no_xi] + a[1] * m[0].upvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m[0].upvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m[0].upvxi[no_u + 2][no_v + 1][no_xi] + m[0].upvxi[no_u][no_v + 3][no_xi] + m[0].upvxi[no_u][no_v + 1][no_xi + 1]);
    psi[3] = 0.5 * (a[0] * (m[0].upvxi[no_u + 2][no_v][no_xi] + m[0].upvxi[no_u][no_v + 2][no_xi] + m[0].upvxi[no_u][no_v][no_xi + 1]) +
        a[1] * (m[0].upvxi[no_u + 3][no_v][no_xi] + m[0].upvxi[no_u + 1][no_v + 2][no_xi] + m[0].upvxi[no_u + 1][no_v][no_xi + 1]) +
        a[2] * (m[0].upvxi[no_u + 2][no_v + 1][no_xi] + m[0].upvxi[no_u][no_v + 3][no_xi] + m[0].upvxi[no_u][no_v + 1][no_xi + 1]) +
        a[3] * 0.5 * (m[0].upvxi[no_u + 4][no_v][no_xi] + m[0].upvxi[no_u][no_v + 4][no_xi] + m[0].upvxi[no_u][no_v][no_xi + 2] + 2 * m[0].upvxi[no_u + 2][no_v + 2][no_xi] + 2 * m[0].upvxi[no_u + 2][no_v][no_xi + 1] + 2 * m[0].upvxi[no_u][no_v + 2][no_xi + 1]));
}

void GR_address(int no_u, int no_v, int no_xi, double* psi, double* a, MMDF* m)
{
    psi[0] = a[0] * m[0].unvxi[no_u][no_v][no_xi] + a[1] * m[0].unvxi[no_u + 1][no_v][no_xi] + a[2] * m[0].unvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m[0].unvxi[no_u + 2][no_v][no_xi] + m[0].unvxi[no_u][no_v + 2][no_xi] + m[0].unvxi[no_u][no_v][no_xi + 1]);
    psi[1] = a[0] * m[0].unvxi[no_u + 1][no_v][no_xi] + a[1] * m[0].unvxi[no_u + 2][no_v][no_xi] + a[2] * m[0].unvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m[0].unvxi[no_u + 3][no_v][no_xi] + m[0].unvxi[no_u + 1][no_v + 2][no_xi] + m[0].unvxi[no_u + 1][no_v][no_xi + 1]);
    psi[2] = a[0] * m[0].unvxi[no_u][no_v + 1][no_xi] + a[1] * m[0].unvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m[0].unvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m[0].unvxi[no_u + 2][no_v + 1][no_xi] + m[0].unvxi[no_u][no_v + 3][no_xi] + m[0].unvxi[no_u][no_v + 1][no_xi + 1]);
    psi[3] = 0.5 * (a[0] * (m[0].unvxi[no_u + 2][no_v][no_xi] + m[0].unvxi[no_u][no_v + 2][no_xi] + m[0].unvxi[no_u][no_v][no_xi + 1]) +
        a[1] * (m[0].unvxi[no_u + 3][no_v][no_xi] + m[0].unvxi[no_u + 1][no_v + 2][no_xi] + m[0].unvxi[no_u + 1][no_v][no_xi + 1]) +
        a[2] * (m[0].unvxi[no_u + 2][no_v + 1][no_xi] + m[0].unvxi[no_u][no_v + 3][no_xi] + m[0].unvxi[no_u][no_v + 1][no_xi + 1]) +
        a[3] * 0.5 * (m[0].unvxi[no_u + 4][no_v][no_xi] + m[0].unvxi[no_u][no_v + 4][no_xi] + m[0].unvxi[no_u][no_v][no_xi + 2] + 2 * m[0].unvxi[no_u + 2][no_v + 2][no_xi] + 2 * m[0].unvxi[no_u + 2][no_v][no_xi + 1] + 2 * m[0].unvxi[no_u][no_v + 2][no_xi + 1]));
}

void G_address(int no_u, int no_v, int no_xi, double* psi, double* a, MMDF* m)
{

    psi[0] = a[0] * m[0].uvxi[no_u][no_v][no_xi] + a[1] * m[0].uvxi[no_u + 1][no_v][no_xi] + a[2] * m[0].uvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m[0].uvxi[no_u + 2][no_v][no_xi] + m[0].uvxi[no_u][no_v + 2][no_xi] + m[0].uvxi[no_u][no_v][no_xi + 1]);
    psi[1] = a[0] * m[0].uvxi[no_u + 1][no_v][no_xi] + a[1] * m[0].uvxi[no_u + 2][no_v][no_xi] + a[2] * m[0].uvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m[0].uvxi[no_u + 3][no_v][no_xi] + m[0].uvxi[no_u + 1][no_v + 2][no_xi] + m[0].uvxi[no_u + 1][no_v][no_xi + 1]);
    psi[2] = a[0] * m[0].uvxi[no_u][no_v + 1][no_xi] + a[1] * m[0].uvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m[0].uvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m[0].uvxi[no_u + 2][no_v + 1][no_xi] + m[0].uvxi[no_u][no_v + 3][no_xi] + m[0].uvxi[no_u][no_v + 1][no_xi + 1]);
    psi[3] = 0.5 * (a[0] * (m[0].uvxi[no_u + 2][no_v][no_xi] + m[0].uvxi[no_u][no_v + 2][no_xi] + m[0].uvxi[no_u][no_v][no_xi + 1]) +
        a[1] * (m[0].uvxi[no_u + 3][no_v][no_xi] + m[0].uvxi[no_u + 1][no_v + 2][no_xi] + m[0].uvxi[no_u + 1][no_v][no_xi + 1]) +
        a[2] * (m[0].uvxi[no_u + 2][no_v + 1][no_xi] + m[0].uvxi[no_u][no_v + 3][no_xi] + m[0].uvxi[no_u][no_v + 1][no_xi + 1]) +
        a[3] * 0.5 * (m[0].uvxi[no_u + 4][no_v][no_xi] + m[0].uvxi[no_u][no_v + 4][no_xi] + m[0].uvxi[no_u][no_v][no_xi + 2] + 2 * m[0].uvxi[no_u + 2][no_v + 2][no_xi] + 2 * m[0].uvxi[no_u + 2][no_v][no_xi + 1] + 2 * m[0].uvxi[no_u][no_v + 2][no_xi + 1]));
}

void Collision(double* center, MMDF* ml, MMDF* mr, double left, double right)
{
    // get the equilibrium variables by collision
    center[0] = left * ml[0].uplus[0] + right * mr[0].uminus[0];
    center[1] = left * ml[0].uplus[1] + right * mr[0].uminus[1];
    center[2] = left * ml[0].vwhole[1] * ml[0].uplus[0] + right * mr[0].vwhole[1] * mr[0].uminus[0];
    center[3] = 0.5 * left * (ml[0].uplus[2] + ml[0].uplus[0] * ml[0].vwhole[2] + ml[0].uplus[0] * ml[0].xi2) +
        0.5 * right * (mr[0].uminus[2] + mr[0].uminus[0] * mr[0].vwhole[2] + mr[0].uminus[0] * mr[0].xi2);
}

void Center_do_nothing_normal(Interface2d* interface, Fluid2d* fluid, int N, int ghost)
{
    // do nothing
}

void Center_all_collision_2D_multi(Interface2d* interface, Block2d block, int i, int j, int num_gauss)
{
    int n = block.N + 2 * block.ghost;
    auto* ml = new MMDF;
    auto* mr = new MMDF;
    double prim_left[4]{}, prim_right[4]{}, convar_left[4]{}, convar_right[4]{}, convar_center[4]{};
    double der1x_left[4]{}, der1x_right[4]{}, der1y_left[4]{}, der1y_right[4]{};

    for (int var = 0; var < 4; var++)
    {
        convar_left[var] = interface[i * n + j].GaussWL[num_gauss][var];
        convar_right[var] = interface[i * n + j].GaussWR[num_gauss][var];
        der1x_left[var] = interface[i * n + j].GaussDer1Lx[num_gauss][var];
        der1y_left[var] = interface[i * n + j].GaussDer1Ly[num_gauss][var];
        der1x_right[var] = interface[i * n + j].GaussDer1Rx[num_gauss][var];
        der1y_right[var] = interface[i * n + j].GaussDer1Ry[num_gauss][var];
    }
    Convar_to_ULambda_2d(prim_left, convar_left);
    Convar_to_ULambda_2d(prim_right, convar_right);
    MMDF_calculate(ml, prim_left);
    MMDF_calculate(mr, prim_right);
    Collision(convar_center, ml, mr, prim_left[0], prim_right[0]);

    double alx[4]{}, arx[4]{};
    A_point(alx, der1x_left, prim_left);
    A_point(arx, der1x_right, prim_right);
    double al0x[4]{}, ar0x[4]{};
    GL_address(0, 0, 0, al0x, alx, ml);
    GR_address(0, 0, 0, ar0x, arx, mr);

    double aly[4]{}, ary[4]{};
    A_point(aly, der1y_left, prim_left);
    A_point(ary, der1y_right, prim_right);
    double al0y[4]{}, ar0y[4]{};
    GL_address(0, 0, 0, al0y, aly, ml);
    GR_address(0, 0, 0, ar0y, ary, mr);
    
    for (int var = 0; var < 4; var++)
    {
        interface[i * n + j].GaussWcenter[num_gauss][var] = convar_center[var];
        interface[i * n + j].GaussCenterDer1x[num_gauss][var] = prim_left[0] * al0x[var] + prim_right[0] * ar0x[var];
        interface[i * n + j].GaussCenterDer1y[num_gauss][var] = prim_left[0] * al0y[var] + prim_right[0] * ar0y[var];
    }

    delete ml;
    delete mr;
}

void Center_all_collision_multi(Interface2d* xinterface, Interface2d* yinterface, Fluid2d* fluid, Block2d block)
{
    int n = block.N + 2 * block.ghost, gausspoint = block.gausspoint;
    int ghost = block.ghost;
//#pragma omp parallel  for
    for (int i = ghost; i < n - ghost + 1; i++)
    {
        for (int j = ghost; j < n - ghost + 1; j++)
        {
            for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
            {
                Center_all_collision_2D_multi(xinterface, block, i, j, num_gauss);
                Center_all_collision_2D_multi(yinterface, block, i, j, num_gauss);
            }
        }
    }
}

double Get_Tau(Interface2d* interface, double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt)
{
    double c1_euler = interface[0].c1;
    double c2_euler = interface[0].c2;
    double C = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right);
    //C+= c2_euler * abs(density_left - density_right) / abs(density_left + density_right);
    //if (C < 10)
    //{
    return c1_euler * dt + dt * C;
}

void GKS2D_non_smooth(Interface2d* interface, int num_gauss, int stage, int i, int j, int n, double dt)
{
    int loc = i * n + j;
    double Flux[2][4];
    double convar_left[4], convar_right[4], convar0[4];
    double der1x_left[4]{}, der1x_right[4]{}, der1y_left[4]{}, der1y_right[4]{};
    auto* ml = new MMDF;
    auto* mr = new MMDF;
    auto* m0 = new MMDF;
    for (int var = 0; var < 4; var++)
    {
        convar_left[var] = interface[loc].GaussWL[num_gauss][var];
        convar_right[var] = interface[loc].GaussWR[num_gauss][var];
        der1x_left[var] = interface[loc].GaussDer1Lx[num_gauss][var];
        der1y_left[var] = interface[loc].GaussDer1Ly[num_gauss][var];
        der1x_right[var] = interface[loc].GaussDer1Rx[num_gauss][var];
        der1y_right[var] = interface[loc].GaussDer1Ry[num_gauss][var];
    }
    // if g0reconstruction_2D_tangent == Center_all_collision_multi
    double prim_left[4], prim_right[4], prim0[4];
    Convar_to_ULambda_2d(prim_left, convar_left);
    Convar_to_ULambda_2d(prim_right, convar_right);
    MMDF_calculate(ml, prim_left);
    MMDF_calculate(mr, prim_right);

    Collision(convar0, ml, mr, prim_left[0], prim_right[0]);
    double alx_t[4]{}, arx_t[4]{};
    A_point(alx_t, der1x_left, prim_left);
    A_point(arx_t, der1x_right, prim_right);
    double al0x_t[4], ar0x_t[4];
    GL_address(0, 0, 0, al0x_t, alx_t, ml);
    GR_address(0, 0, 0, ar0x_t, arx_t, mr);
    // w_y
    double aly_t[4], ary_t[4];
    A_point(aly_t, der1y_left, prim_left);
    A_point(ary_t, der1y_right, prim_right);
    double al0y_t[4], ar0y_t[4];
    GL_address(0, 0, 0, al0y_t, aly_t, ml);
    GR_address(0, 0, 0, ar0y_t, ary_t, mr);
    for (int var = 0; var < 4; ++var)
    {
        interface[loc].GaussCenterDer1x[num_gauss][var] = prim_left[0] * al0x_t[var] + prim_right[0] * ar0x_t[var];
        interface[loc].GaussCenterDer1y[num_gauss][var] = prim_left[0] * al0y_t[var] + prim_right[0] * ar0y_t[var];
    }
    // calculate flux
    Convar_to_ULambda_2d(prim0, convar0);
    // get the coefficient of time intergation factors
    double tau = 0, tau_num;
    tau_num = Get_Tau(interface, prim_left[0], prim_right[0], prim0[0], prim_left[3], prim_right[3], prim0[3], dt);
    double eta = exp(-dt / tau_num), t[10]{};
   
    t[0] = tau_num * (1 - eta); // this refers glu, gru part
    t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
    t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
    // then, equ part time coefficient for gks 2nd
    t[3] = tau_num * eta + dt - tau_num; //this refers g0u part
    t[4] = tau_num * (tau_num - eta * (dt + tau_num) - tau * (eta - 1)) - dt * tau; //this refers a0uu part
    t[5] = 0.5 * dt * dt - tau * tau_num * (eta - 1) - tau * dt;
    
    double unit[4] { 1, 0.0, 0.0, 0.0 }, glu[4], gru[4];
    GL_address(1, 0, 0, glu, unit, ml); // (1 0 0) get GLu
    GR_address(1, 0, 0, gru, unit, mr); // (1 0 0) get GRu
    for (int var = 0; var < 4; var++)
    {
        Flux[0][var] = prim_left[0] * t[0] * glu[var] + prim_right[0] * t[0] * gru[var];
    }
    //now the equ part added, m0 term added
    MMDF_calculate(m0, prim0);
    double g0u[4];
    G_address(1, 0, 0, g0u, unit, m0);
    for (int var = 0; var < 4; var++)
    {
        Flux[0][var] = Flux[0][var] + prim0[0] * t[3] * g0u[var];
    }

    double alx[4];
    A_point(alx, der1x_left, prim_left);
    double alxuul[4];
    GL_address(2, 0, 0, alxuul, alx, ml);

    double arx[4];
    A_point(arx, der1x_right, prim_right);
    double arxuur[4];
    GR_address(2, 0, 0, arxuur, arx, mr);

    double aly[4];
    A_point(aly, der1y_left, prim_left);
    double alyuvl[4];
    GL_address(1, 1, 0, alyuvl, aly, ml);

    double ary[4];
    A_point(ary, der1y_right, prim_right);
    double aryuvr[4];
    GR_address(1, 1, 0, aryuvr, ary, mr);
    // t1 part
    for (int var = 0; var < 4; var++)
    {
        Flux[0][var] = Flux[0][var] + t[1] * (prim_left[0] * (alxuul[var] + alyuvl[var]) + prim_right[0] * (arxuur[var] + aryuvr[var]));
    }
    //for t[2] Aru,Alu part
    double alxu[4]{}, alyv[4]{}, arxu[4]{}, aryv[4]{};

    //take <u> moment for al, ar
    G_address(1, 0, 0, alxu, alx, ml);
    G_address(1, 0, 0, arxu, arx, mr);
    G_address(0, 1, 0, alyv, aly, ml);
    G_address(0, 1, 0, aryv, ary, mr);

    double Al[4]{}, Ar[4]{}, der_AL[4]{}, der_AR[4]{};

    //using compatability condition to get the time derivative
    for (int var = 0; var < 4; var++)
    {
        der_AL[var] = -prim_left[0] * (alxu[var] + alyv[var]);
        der_AR[var] = -prim_right[0] * (arxu[var] + aryv[var]);
    }
    // solve the coefficient martix b=ma
    A_point(Al, der_AL, prim_left);
    A_point(Ar, der_AR, prim_right);

    //to obtain the Alu and Aru
    double Alul[4]{}, Arur[4]{};
    GL_address(1, 0, 0, Alul, Al, ml);
    GR_address(1, 0, 0, Arur, Ar, mr);
    // t2 part
    for (int var = 0; var < 4; var++)
    {
        Flux[0][var] = Flux[0][var] + t[2] * (prim_left[0] * Alul[var] + prim_right[0] * Arur[var]);
    }
    // for t[4] a0xuu part
    double a0x[4]{}, derx[4]{};
    for (int var = 0; var < 4; var++)
    {
        derx[var] = interface[loc].GaussCenterDer1x[num_gauss][var]; //Only the averaged value for g0
    }
    //solve the microslope
    A_point(a0x, derx, prim0);
    //a0x <u> moment
    double a0xu[4];
    G_address(1, 0, 0, a0xu, a0x, m0); //get a0xu, used for the following determination of derA0, and then A0
    //a0x <u^2> moment
    double a0xuu[4];
    G_address(2, 0, 0, a0xuu, a0x, m0);

    double a0y[4]{}, dery[4]{};
    for (int var = 0; var < 4; var++)
    {
        dery[var] = interface[loc].GaussCenterDer1y[num_gauss][var];
    }
    A_point(a0y, dery, prim0);
    double a0yv[4];
    G_address(0, 1, 0, a0yv, a0y, m0); //get a0yv, used for the following determination of derA0, and then A0
    //a0x <u^2> moment
    double a0yuv[4];
    G_address(1, 1, 0, a0yuv, a0y, m0);
    // t4 part
    for (int var = 0; var < 4; var++)
    {
        Flux[0][var] = Flux[0][var] + prim0[0] * t[4] * (a0xuu[var] + a0yuv[var]);
    }
    // for t[5] A0u part
    double derA0[4];
    for (int var = 0; var < 4; var++)
    {
        derA0[var] = -prim0[0] * (a0xu[var] + a0yv[var]);
    }
    double A0[4];
    A_point(A0, derA0, prim0);
    double A0u[4];
    G_address(1, 0, 0, A0u, A0, m0);
    // t5 part
    for (int var = 0; var < 4; var++)
    {
        Flux[0][var] = Flux[0][var] + prim0[0] * t[5] * (A0u[var]);
    }

    // half time increment
    double dt2 = 0.5 * dt; // the following is dt2
    tau_num = Get_Tau(interface, prim_left[0], prim_right[0], prim0[0], prim_left[3], prim_right[3], prim0[3], dt2);
    eta = exp(-dt2 / tau_num);
    // non equ part time coefficient for gks_2nd algorithm
    t[0] = tau_num * (1 - eta); // this refers glu, gru part
    t[1] = tau_num * (eta * (dt2 + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
    t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
    // then, equ part time coefficient for gks 2nd
    t[3] = tau_num * eta + dt2 - tau_num; //this refers g0u part
    t[4] = tau_num * (tau_num - eta * (dt2 + tau_num) - tau * (eta - 1)) - dt2 * tau; //this refers a0uu part
    t[5] = 0.5 * dt2 * dt2 - tau * tau_num * (eta - 1) - tau * dt2; //this refers A0u part

    for (int var = 0; var < 4; var++)
    {
        // t0 part
        Flux[1][var] = t[0] * (prim_left[0] * glu[var] + prim_right[0] * gru[var]);
        // t1 part
        Flux[1][var] = Flux[1][var] + t[1] * (prim_left[0] * (alxuul[var] + alyuvl[var]) + prim_right[0] * (arxuur[var] + aryuvr[var]));
        // t2 part
        Flux[1][var] = Flux[1][var] + t[2] * (prim_left[0] * Alul[var] + prim_right[0] * Arur[var]);
        // t3 part
        Flux[1][var] = Flux[1][var] + prim0[0] * t[3] * g0u[var];
        // t4 part
        Flux[1][var] = Flux[1][var] + prim0[0] * t[4] * (a0xuu[var] + a0yuv[var]);
        // t5 part
        Flux[1][var] = Flux[1][var] + prim0[0] * t[5] * (A0u[var]);
    }

    for (int var = 0; var < 4; var++)
    {
        interface[loc].flux[stage][num_gauss][var] = (4.0 * Flux[1][var] - Flux[0][var]);
        interface[loc].der_flux[stage][num_gauss][var] = 4.0 * (Flux[0][var] - 2.0 * Flux[1][var]);
    }

    delete ml;
    delete mr;
    delete m0;
}

void Calculate_flux(Interface2d* xinterface, Interface2d* yinterface, Fluid2d* fluid, Block2d block, int stage, void (*flux_function_2d)(Interface2d*, int, int, int, int, int, double))
{
    int n = block.N + 2 * block.ghost, N = block.N, ghost = block.ghost;
    int gausspoint = block.gausspoint;
#pragma omp parallel  for
    for (int i = ghost; i < N + ghost + 1; i++)
    {
        for (int j = ghost; j < N + ghost; j++)
        {
            for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
            {
                flux_function_2d(xinterface, num_gauss, stage, i, j, n, block.dt);
                // calculate the final flux in xfluxes,  by the reconstructed variables in xinterfaces; the same for y
            }
        }
    }
#pragma omp parallel  for
    for (int i = ghost; i < N + ghost; i++)
    {
        for (int j = ghost; j < N + ghost + 1; j++)
        {
            for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
            {
                flux_function_2d(yinterface, num_gauss, stage, i, j, n, block.dt);
            }
        }
    }
}

void Local_to_Global(Interface2d* interface, Cell2d* cell)
{
    int stage = cell[0].stage;
    int num_gauss = cell[0].num_gauss;
    int loc = cell[0].loc;
    double temp[2]{};
    temp[0] = interface[loc].Flux[stage][num_gauss][1];
    temp[1] = interface[loc].Flux[stage][num_gauss][2];

    interface[loc].Flux[stage][num_gauss][1] = -temp[1];
    interface[loc].Flux[stage][num_gauss][2] = temp[0];
}

void Update_with_gauss(Interface2d* xinterface, Interface2d* yinterface, Fluid2d* fluid, Block2d block, int stage)
{
    int N = block.N, ghost = block.ghost;
    int n = block.N + 2 * block.ghost, gausspoint = block.gausspoint;
    auto* cell = new Cell2d;
    for (int i = ghost; i < N + ghost + 1; i++)
    {
        for (int j = ghost; j < N + ghost + 1; j++)
        {
            int loc = i * n + j;
            cell[0].loc = loc;
            for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
            {
                cell[0].num_gauss = num_gauss;
                cell[0].stage = stage;
                Local_to_Global(yinterface, cell);
            }
        }
    }
    for (int i = ghost; i < N + ghost; i++)
    {
        for (int j = ghost; j < N + ghost; j++)
        {
            int loc = i * n + j;

            for (int var = 0; var < 4; var++)
            {
                fluid[loc].W[var] = fluid[loc].W_old[var]; //get the Wn from convar_old
                double total_flux = 0.0;
                for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
                {
                    total_flux += block.dx * yinterface[i * n + j].Flux[stage][num_gauss][var];
                    total_flux += -block.dx * yinterface[i * n + j + 1].Flux[stage][num_gauss][var];
                    total_flux += block.dx * xinterface[i * n + j].Flux[stage][num_gauss][var];
                    total_flux += -block.dx * xinterface[(i+1)*n + j].Flux[stage][num_gauss][var];
                }
                double area = block.dx * block.dx;
                fluid[loc].W[var] += total_flux / area;
                // calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
                // in 2d, the total flux be updated by the line-averaged flux of four interface
            }
            fluid[loc].prim[0] = fluid[loc].W[0];
            fluid[loc].prim[1] = fluid[loc].W[1] / fluid[loc].prim[0];
            fluid[loc].prim[2] = fluid[loc].W[2] / fluid[loc].prim[0];
            fluid[loc].prim[3] = 0.4 * (fluid[loc].W[3] - 0.5 * fluid[loc].prim[0] * (pow(fluid[loc].prim[1], 2) + pow(fluid[loc].prim[2], 2)));
        }
    }
}

void Update_variables(Interface2d* xinterface, Interface2d* yinterface, Fluid2d* fluid, Block2d block, int stage)
{
    int N = block.N, ghost = block.ghost;
    int gausspoint = block.gausspoint, n = N + 2 * ghost;
    for (int i = ghost; i < N + ghost + 1; i++)
    {
        for (int j = ghost; j < N + ghost; j++)
        {
            for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
            {
                for (int var = 0; var < 4; var++)
                {
                    double Flux = 0.0;
                    for (int k = 0; k < stage + 1; k++)
                    {

                        Flux = Flux
                            + block.gauss_weight[num_gauss] *
                            (block.timecoefficient[stage][k][0] * xinterface[i*n+j].flux[stage][num_gauss][var]
                                + block.timecoefficient[stage][k][1] * xinterface[i * n + j].der_flux[stage][num_gauss][var]);

                    }
                    xinterface[i * n + j].Flux[stage][num_gauss][var] = Flux;
                    // calculate the final flux of the interface, in x in xfluxes, by the obtained the flux and its derivative (f, def, der2f) at guass points, in xfluxes, and the corresponding weight factors
                    // calculate by several stages according to the time marching method. same for yfluxes
                }
            }
        }
    }

    for (int i = ghost; i < N + ghost; i++)
    {
        for (int j = ghost; j < N + ghost + 1; j++)
        {
            for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
            {
                for (int var = 0; var < 4; var++)
                {
                    double Flux = 0.0;
                    for (int k = 0; k < stage + 1; k++)
                    {

                        Flux = Flux
                            + block.gauss_weight[num_gauss] *
                            (block.timecoefficient[stage][k][0] * yinterface[i * n + j].flux[stage][num_gauss][var]
                                + block.timecoefficient[stage][k][1] * yinterface[i * n + j].der_flux[stage][num_gauss][var]);

                    }
                    yinterface[i * n + j].Flux[stage][num_gauss][var] = Flux;
                    // calculate the final flux of the interface, in x in xfluxes, by the obtained the flux and its derivative (f, def, der2f) at guass points, in xfluxes, and the corresponding weight factors
                    // calculate by several stages according to the time marching method. same for yfluxes
                }
            }
        }
    }

    Update_with_gauss(xinterface, yinterface, fluid, block, stage);
}

void free_boundary(Fluid2d* fluid, Block2d block)
{
    int n = block.N + 2 * block.ghost;
    int ghost = block.ghost, N = block.N;
    // left boudnary
    for (int i = ghost - 1; i >= 0; i--)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                fluid[i * n + j].W[k] = fluid[(i + 1) * n + j].W[k];
            }
        }
    }
    // rigth boundary
    for (int i = ghost + N; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                fluid[i * n + j].W[k] = fluid[(i - 1) * n + j].W[k];
            }
        }
    }
    // down boundary
    for (int j = ghost - 1; j >= 0; j--)
    {
        for (int i = 0; i < n; i++)
        {
            for (int k = 0; k < 4; k++)
            {
                fluid[i * n + j].W[k] = fluid[i * n + j + 1].W[k];
            }
        }
    }
    // up boundary
    for (int j = ghost + N; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int k = 0; k < 4; k++)
            {
                fluid[i * n + j].W[k] = fluid[i * n + j - 1].W[k];
            }
        }
    }
}

double Dtx(double dtx, Block2d block, double density, double u, double v, double pressure)
{
    double tmp;
    tmp = sqrt(u * u + v * v) + sqrt(Gamma * pressure / density);
    if (tmp > block.CFL * block.dx / dtx)
    {
        dtx = block.CFL * block.dx / tmp;
    }
    return dtx;
}

double Get_CFL(Fluid2d* fluid, Block2d block)
{
    double dt = block.dx, density, u, v, pressure;
    int n = block.N + 2 * block.ghost, ghost = block.ghost;
    for (int i = ghost; i < n - ghost; i++)
    {
        for (int j = ghost; j < n - ghost; j++)
        {
            density = fluid[i * n + j].prim[0];
            u = fluid[i * n + j].prim[1];
            v = fluid[i * n + j].prim[2];
            pressure = fluid[i * n + j].prim[3];
            dt = Dtx(dt, block, density, u, v, pressure);
        }
    }
    if (block.tnow + dt > block.tstop)
    {
        dt = block.tstop - block.tnow + tol;
    }
    return dt;
}

void S1O2_2D(Block2d& block)
{
    block.stages = 1;
    block.timecoefficient[0][0][0] = 1.0;
    block.timecoefficient[0][0][1] = 1.0;
}

void S2O4_2D(Block2d& block)
{
    block.stages = 2;
    block.timecoefficient[0][0][0] = 0.5;
    block.timecoefficient[0][0][1] = 1.0 / 8.0;
    block.timecoefficient[1][0][0] = 1.0;
    block.timecoefficient[1][1][0] = 0.0;
    block.timecoefficient[1][0][1] = 1.0 / 6.0;
    block.timecoefficient[1][1][1] = 1.0 / 3.0;
}

void output2d(Interface2d* interface, Fluid2d* fluid, Block2d block)
{
    int N = block.N, ghost = block.ghost;
    int n = block.N + 2 * block.ghost;
    ofstream ResultFile;
    const char* filePath_R = "/Users/hongzhang/Desktop/data/R2d.txt";
    const char* filePath_U = "/Users/hongzhang/Desktop/data/U2d.txt";
    const char* filePath_V = "/Users/hongzhang/Desktop/data/V2d.txt";
    const char* filePath_p = "/Users/hongzhang/Desktop/data/p2d.txt";

    ResultFile.open(filePath_R);
    for (int i = ghost; i < N + ghost; i++)
    {
        for (int j = ghost; j < N + ghost; j++)
        {
            ResultFile << fluid[i * n + j].W[0] << endl;
        }
    }
    ResultFile.close();

    ResultFile.open(filePath_U);
    for (int i = ghost; i < N + ghost; i++)
    {
        for (int j = ghost; j < N + ghost; j++)
        {
            ResultFile << fluid[i * n + j].W[1] << endl;
        }
    }
    ResultFile.close();

    ResultFile.open(filePath_V);
    for (int i = ghost; i < N + ghost; i++)
    {
        for (int j = ghost; j < N + ghost; j++)
        {
            ResultFile << fluid[i * n + j].W[2] << endl;
        }
    }
    ResultFile.close();

    ResultFile.open(filePath_p);
    for (int i = ghost; i < N + ghost; i++)
    {
        for (int j = ghost; j < N + ghost; j++)
        {
            ResultFile << fluid[i * n + j].W[3] << endl;
        }
    }
    ResultFile.close();
}

void (*boundary_type_2d)(Fluid2d*, Block2d);

void (*initial_stage_2d)(Block2d&);

void (*cell_reconstruction_2d)(Interface2d*, Interface2d*, Fluid2d*, Block2d);

void (*g0reconstruction_2D_tangent)(Interface2d*, Interface2d*, Fluid2d*, Block2d);

void (*flux_function_2d)(Interface2d*, int, int, int, int, int, double);

void Riemann_problem_2d()
{
    Block2d block;
    block.N = 300;
    block.ghost = 2;
    block.CFL = 0.5;
    block.tnow = 0.0;
    int n = pow(block.N + 2 * block.ghost, 2);
    auto* fluid = new Fluid2d[n];
    auto* xinterface = new Interface2d[n];
    auto* yinterface = new Interface2d[n];

    initial_stage_2d = S1O2_2D;
    boundary_type_2d = free_boundary;
    cell_reconstruction_2d = VanLeer;
    flux_function_2d = GKS2D_non_smooth;
    g0reconstruction_2D_tangent = Center_all_collision_multi;

    // time coefficient
    (*initial_stage_2d)(block);
    // initial variables;
    
    block.dx = 1.0 / block.N;
    block.x = 0.7;
    block.y = 0.7;
    block.tstop = 0.6;
    double zone1[4]{ 0.138, 1.206, 1.206, 0.029 };
    double zone2[4]{ 0.5323, 0, 1.206, 0.3 };
    double zone3[4]{ 1.5, 0, 0, 1.5 };
    double zone4[4]{ 0.5323, 1.206, 0, 0.3 };
    /*
    fluid[0].dx = 1.0 / N;
    double x = 0.5, y = 0.5, tstop = 0.2;
    double zone1[] = { 0.8, 0, 0, 1.0 };
    double zone2[]{ 1.0, 0, 0.7276, 1.0 };
    double zone3[]{ 0.5313, 0, 0, 0.4 };
    double zone4[]{ 1.0, 0.7276, 0, 1.0 };
    */
    /*
    fluid[0].dx = 2.0 / N;
    double x = 1.0, y = 1.0, tstop = 0.4;
    double zone1[] = { 1.0, -0.75, 0.5, 1.0 };
    double zone2[]{ 3.0, -0.75, -0.5, 1.0 };
    double zone3[]{ 1.0, 0.75, -0.5, 1.0 };
    double zone4[]{ 2.0, 0.75, 0.5, 1.0 };
    */
    Initial_variables(fluid, zone1, zone2, zone3, zone4, block);
    // initial gauss points
    block.gausspoint = 1;
    SetGaussPoint(block);

    while (block.tnow < block.tstop)
    {
        // copy the variables in step n
        Copy_array(fluid, block);
        // get the time increment
        block.dt = Get_CFL(fluid, block);
        // time iteration
        for (int stage = 0; stage < block.stages; stage++)
        {
            // boundary condition
            (*boundary_type_2d)(fluid, block);
            
            // reconstruction within cell
            (*cell_reconstruction_2d)(xinterface, yinterface, fluid, block);

            // g0 reconstruction
            (*g0reconstruction_2D_tangent)(xinterface, yinterface, fluid, block);

            // Calculate flux
            Calculate_flux(xinterface, yinterface, fluid, block, stage, flux_function_2d);

            // update variables
            Update_variables(xinterface, yinterface, fluid, block, stage);
        }
        block.tnow += block.dt;
        block.step += 1;
        if (block.step % 10 == 0)
        {
            cout << "Progress degree: " << block.tnow / block.tstop * 100 << "%" << endl;
        }
    }
    output2d(xinterface, fluid, block);
}

void Accuracy_test_2d()
{

}
