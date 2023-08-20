/*
 * One dimension second order gas-kinetic scheme
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

#include "GKS1D.h"
#define K 4
#define gamma 1.4
#define pi acos(-1)
#define tol 1e-20

double Alpha1d(double lambda, double u)
{
    return erfc(sqrt(lambda) * u);
}

double Beta1d(double lambda, double u)
{
    return exp(-lambda * u * u) / sqrt(pi * lambda);
}

void GL1d(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
    psi[0] = a[0] * m.upxi[no_u][no_xi] + a[1] * m.upxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]);
    psi[1] = a[0] * m.upxi[no_u + 1][no_xi] + a[1] * m.upxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]);
    psi[2] = 0.5 * (a[0] * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]) +
        a[1] * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]) +
        a[2] * 0.5 * (m.upxi[no_u + 4][no_xi] + m.upxi[no_u][no_xi + 2] + 2 * m.upxi[no_u + 2][no_xi + 1]));
}

void GR1d(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
    psi[0] = a[0] * m.unxi[no_u][no_xi] + a[1] * m.unxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]);
    psi[1] = a[0] * m.unxi[no_u + 1][no_xi] + a[1] * m.unxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]);
    psi[2] = 0.5 * (a[0] * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]) +
        a[1] * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]) +
        a[2] * 0.5 * (m.unxi[no_u + 4][no_xi] + m.unxi[no_u][no_xi + 2] + 2 * m.unxi[no_u + 2][no_xi + 1]));

}

void G1d(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
    psi[0] = a[0] * m.uxi[no_u][no_xi] + a[1] * m.uxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]);
    psi[1] = a[0] * m.uxi[no_u + 1][no_xi] + a[1] * m.uxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]);
    psi[2] = 0.5 * (a[0] * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]) +
        a[1] * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]) +
        a[2] * 0.5 * (m.uxi[no_u + 4][no_xi] + m.uxi[no_u][no_xi + 2] + 2 * m.uxi[no_u + 2][no_xi + 1]));

}

double sign(double x)
{
    return (x < 0) ? -1 : (x > 0) ? 1 : 0;
}

double lambda(double R, double RU, double RE)
{
    return 1.25 / (RE / R - 0.5 * (RU / R) * (RU / R));
}

void Convar_to_ULambda_1d(double* primvar, double* convar)
{
    primvar[0] = convar[0];
    primvar[1] = convar[1] / convar[0];
    primvar[2] = lambda(convar[0], convar[1], convar[2]);
}

void Convar_to_Char(double character[], const double base[], const double convar[])
{
    double c = sqrt(1.4 * base[2] / base[0]);
    double alfa = 0.4 / (2.0 * c * c);
    double u = base[1];
    double s[3][3];
    s[0][0] = alfa * (0.5 * u * u + u * c / 0.4);
    s[0][1] = alfa * (-u - c / 0.4);
    s[0][2] = alfa;
    s[1][0] = alfa * (-u * u + 2.0 * c * c / 0.4);
    s[1][1] = alfa * 2.0 * u;
    s[1][2] = -2.0 * alfa;
    s[2][0] = alfa * (0.5 * u * u - u * c / 0.4);
    s[2][1] = alfa * (-u + c / 0.4);
    s[2][2] = alfa;

    for (int i = 0; i < 3; i++)
    {
        character[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            character[i] = character[i] + s[i][j] * convar[j];
        }
    }
}

void Char_to_Convar(double var[], const double base[], const double W[])
{
    double r = 1.4;
    double c = sqrt(r * base[2] / base[0]);

    double u = base[1];
    double	h = 0.5 * u * u + c * c / (r - 1.0);
    double s[3][3];
    s[0][0] = 1.0;
    s[0][1] = 1.0;
    s[0][2] = 1.0;
    s[1][0] = u - c;
    s[1][1] = u;
    s[1][2] = u + c;
    s[2][0] = h - u * c;
    s[2][1] = u * u / 2.0;
    s[2][2] = h + u * c;

    for (int i = 0; i < 3; i++)
    {
        var[i] = 0;
        for (int j = 0; j < 3; j++)
        {
            var[i] = var[i] + s[i][j] * W[j];
        }
    }
}

void VanLeer(Block1d block, Fluid* fluid, Interface* interface)
{
    /*
     * cell left correspond to the interface right
     * cell right correspond to the interface left
     * */
    double s, r;

    // cell reconstruction
    for (int i = block.ghost - 1; i <= block.nodex + block.ghost; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            s = (fluid[i + 1].W[j] - fluid[i].W[j]) / block.dx;
            r = (fluid[i].W[j] - fluid[i - 1].W[j]) / block.dx;
            fluid[i].DerL[j] = (sign(s) + sign(r)) * s * r / (abs(s) + abs(r) + tol);
            fluid[i].DerR[j] = fluid[i].DerL[j];
            fluid[i].WL[j] = fluid[i].W[j] - 0.5 * block.dx * fluid[i].DerL[j];
            fluid[i].WR[j] = fluid[i].W[j] + 0.5 * block.dx * fluid[i].DerL[j];
        }
        double flag_l = fluid[i].WL[2] - 0.5 * fluid[i].WL[1] * fluid[i].WL[1] / fluid[i].WL[0];
        double flag_r = fluid[i].WR[2] - 0.5 * fluid[i].WR[1] * fluid[i].WR[1] / fluid[i].WR[0];
        if (flag_l <= 0 || flag_r <= 0)
        {
            for (int j = 0; j <= 2; j++)
            {
                fluid[i].DerL[j] = 0;
                fluid[i].DerR[j] = 0;
                fluid[i].WL[j] = fluid[i].W[j];
                fluid[i].WR[j] = fluid[i].W[j];
            }
        }
    }
    // interface reconstruction
    for (int i = block.ghost; i <= block.nodex + block.ghost; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            interface[i].WL[j] = fluid[i - 1].WR[j];
            interface[i].WR[j] = fluid[i].WL[j];
            interface[i].DerL[j] = fluid[i - 1].DerR[j];
            interface[i].DerR[j] = fluid[i].DerL[j];
        }
        //cout << interface[i].WL[0] << endl;
    }
}

void WENO5_AO_left(double& var, double& der1, double w0, double wp1, double wp2, double wn1, double wn2, Fluid* fluid, double h)
{
    double dhi = 0.85;
    double dlo = 0.85;
    //-- - parameter of WENO-- -
    double beta[4], d[4], ww[4], alpha[4];
    double epsilonW = 1e-8;
    //-- - intermediate parameter-- -
    double p[4], px[4], pxx[4];
    double sum_alpha;

    //three small stencil
    d[0] = (1 - dhi) * (1 - dlo) / 2.0;
    d[1] = (1 - dhi) * dlo;
    d[2] = (1 - dhi) * (1 - dlo) / 2.0;
    //one big stencil
    d[3] = dhi;

    //cout << "here" << endl;
    beta[0] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w0), 2);
    beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
    beta[2] = 13.0 / 12.0 * pow((w0 - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w0 - 4 * wp1 + wp2), 2);

    beta[3] = (1.0 / 5040.0) * (231153 * w0 * w0 + 104963 * wn1 * wn1 + 6908 * wn2 * wn2 -
        38947 * wn2 * wp1 + 104963 * wp1 * wp1 +
        wn1 * (-51001 * wn2 + 179098 * wp1 - 38947 * wp2) -
        3 * w0 * (99692 * wn1 - 22641 * wn2 + 99692 * wp1 - 22641 * wp2) +
        8209 * wn2 * wp2 - 51001 * wp1 * wp2 + 6908 * wp2 * wp2);

    double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

    sum_alpha = 0.0;
    for (int i = 0; i < 4; i++)
    {
        double global_div = tau5 / (beta[i] + epsilonW);
        alpha[i] = d[i] * (1 + global_div * global_div);
        sum_alpha += alpha[i];
    }

    for (int k = 0; k < 4; k++)
    {
        ww[k] = alpha[k] / sum_alpha;
    }
    //-- - candidate polynomial-- -
    p[0] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w0;
    p[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w0 - 1.0 / 6.0 * wp1;
    p[2] = 11.0 / 6.0 * w0 - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
    p[3] = (1.0 / 60.0) * (47 * w0 + 27 * wn1 - 3 * wn2 - 13 * wp1 + 2 * wp2);

    px[0] = (w0 - wn1) / h;
    px[1] = (w0 - wn1) / h;
    px[2] = -((2 * w0 - 3 * wp1 + wp2) / h);
    px[3] = (15 * w0 - 15 * wn1 + wn2 - wp1) / (12 * h);

    pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
    pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
    pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
    pxx[3] = ((-8 * w0 + 2 * wn1 + wn2 + 6 * wp1 - wp2) / (4 * h * h));

    //-- - combination-- -
    var = 0.0;
    der1 = 0.0;
    double final_weight[4];
    final_weight[3] = ww[3] / d[3];
    for (int k = 0; k < 3; k++)
    {
        final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
    }

    for (int k = 0; k < 4; k++)
    {
        var += final_weight[k] * p[k];
        der1 += final_weight[k] * px[k];
    }
}

void WENO5_AO_right(double& var, double& der1, double w0, double wp1, double wp2, double wn1, double wn2, Fluid* fluid, double h)
{
    double dhi = 0.85;
    double dlo = 0.85;
    //-- - parameter of WENO-- -
    double beta[4], d[4], ww[4], alpha[4];
    double epsilonW = 1e-8;

    //-- - intermediate parameter-- -
    double p[4], px[4], pxx[4];
    double sum_alpha;

    //three small stencil
    d[0] = (1 - dhi) * (1 - dlo) / 2.0;
    d[1] = (1 - dhi) * dlo;
    d[2] = (1 - dhi) * (1 - dlo) / 2.0;
    //one big stencil
    d[3] = dhi;

    beta[0] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w0), 2);
    beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
    beta[2] = 13.0 / 12.0 * pow((w0 - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w0 - 4 * wp1 + wp2), 2);

    beta[3] = (1.0 / 5040.0) * (231153 * w0 * w0 + 104963 * wn1 * wn1 + 6908 * wn2 * wn2 -
        38947 * wn2 * wp1 + 104963 * wp1 * wp1 +
        wn1 * (-51001 * wn2 + 179098 * wp1 - 38947 * wp2) -
        3 * w0 * (99692 * wn1 - 22641 * wn2 + 99692 * wp1 - 22641 * wp2) +
        8209 * wn2 * wp2 - 51001 * wp1 * wp2 + 6908 * wp2 * wp2);

    double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

    sum_alpha = 0.0;
    for (int i = 0; i < 4; i++)
    {
        double global_div = tau5 / (beta[i] + epsilonW);
        alpha[i] = d[i] * (1 + global_div * global_div);
        sum_alpha += alpha[i];
    }

    for (int k = 0; k < 4; k++)
    {
        ww[k] = alpha[k] / sum_alpha;
    }
    //-- - candidate polynomial-- -

    p[0] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w0;
    p[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w0 + 1.0 / 3.0 * wp1;
    p[2] = 1.0 / 3.0 * w0 + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
    p[3] = (1.0 / 60.0) * (47 * w0 - 13 * wn1 + 2 * wn2 + 27 * wp1 - 3 * wp2);

    px[0] = (2 * w0 - 3 * wn1 + wn2) / h;
    px[1] = (-w0 + wp1) / h;
    px[2] = (-w0 + wp1) / h;
    px[3] = (-15 * w0 + wn1 + 15 * wp1 - wp2) / (12 * h);

    pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
    pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
    pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
    pxx[3] = (-8 * w0 + 6 * wn1 - wn2 + 2 * wp1 + wp2) / (4 * h * h);

    //-- - combination-- -
    var = 0.0;
    der1 = 0.0;
    double final_weight[4];
    final_weight[3] = ww[3] / d[3];
    for (int k = 0; k < 3; k++)
    {
        final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
    }

    for (int k = 0; k < 4; k++)
    {
        var += final_weight[k] * p[k];
        der1 += final_weight[k] * px[k];
    }
}

void WENO_AO_cal(int i, Fluid* fluid, Block1d block)
{
    double w[3], wp1[3], wp2[3], wn1[3], wn2[3], var[3]{}, der1[3]{};
    double base_left[3], base_right[3];
    double W[3], Wp1[3], Wp2[3], Wn1[3], Wn2[3];
    double ConW[3]{}, DerConW[3]{};
    for (int j = 0; j < 3; j++)
    {
        W[j] = fluid[i].W[j];
        Wp1[j] = fluid[i + 1].W[j];
        Wp2[j] = fluid[i + 2].W[j];
        Wn1[j] = fluid[i - 1].W[j];
        Wn2[j] = fluid[i - 2].W[j];

        w[j] = fluid[i].prim[j];
        wp1[j] = fluid[i + 1].prim[j];
        wp2[j] = fluid[i + 2].prim[j];
        wn1[j] = fluid[i - 1].prim[j];
        wn2[j] = fluid[i - 2].prim[j];

        base_left[j] = 0.5 * (wn1[j] + w[j]);
        base_right[j] = 0.5 * (w[j] + wp1[j]);
    }

    // correspond to conservative variables
    if (block.variable_type == "conservative")
    {
        // cell left side -- interface right side
        for (int j = 0; j < 3; j++)
        {
            WENO5_AO_left(fluid[i].WL[j], fluid[i].DerL[j], W[j], Wp1[j], Wp2[j], Wn1[j], Wn2[j], fluid, block.dx);
        }

        // cell right side -- interface left side
        for (int j = 0; j < 3; j++)
        {
            WENO5_AO_right(fluid[i].WR[j], fluid[i].DerR[j], W[j], Wp1[j], Wp2[j], Wn1[j], Wn2[j], fluid, block.dx);
        }
    }
    // correspond to characteristic variables
    else
    {
        // cell left side -- interface right side
        Convar_to_Char(wn2, base_left, Wn2);
        Convar_to_Char(wn1, base_left, Wn1);
        Convar_to_Char(w, base_left, W);
        Convar_to_Char(wp1, base_left, Wp1);
        Convar_to_Char(wp2, base_left, Wp2);

        for (int j = 0; j < 3; j++)
        {
            WENO5_AO_left(var[j], der1[j], w[j], wp1[j], wp2[j], wn1[j], wn2[j], fluid, block.dx);
        }
        Char_to_Convar(ConW, base_left, var);
        Char_to_Convar(DerConW, base_left, der1);
        for (int j = 0; j < 3; j++)
        {
            fluid[i].WL[j] = ConW[j];
            fluid[i].DerL[j] = DerConW[j];
        }
        // cell right side -- interface left side
        Convar_to_Char(wn2, base_right, Wn2);
        Convar_to_Char(wn1, base_right, Wn1);
        Convar_to_Char(w, base_right, W);
        Convar_to_Char(wp1, base_right, Wp1);
        Convar_to_Char(wp2, base_right, Wp2);

        for (int j = 0; j < 3; j++)
        {
            WENO5_AO_right(var[j], der1[j], w[j], wp1[j], wp2[j], wn1[j], wn2[j], fluid, block.dx);
        }
        Char_to_Convar(ConW, base_right, var);
        Char_to_Convar(DerConW, base_right, der1);
        for (int j = 0; j < 3; j++)
        {
            fluid[i].WR[j] = ConW[j];
            fluid[i].DerR[j] = DerConW[j];
        }
    }

    // check if order reduce
    // if lambda < 0, then the order reduce to the first order
    double flag_l = fluid[i].WL[2] - 0.5 * fluid[i].WL[1] * fluid[i].WL[1] / fluid[i].WL[0];
    double flag_r = fluid[i].WR[2] - 0.5 * fluid[i].WR[1] * fluid[i].WR[1] / fluid[i].WR[0];
    if (flag_l <= 0 || flag_r <= 0 || fluid[i].WL[0] <= 0 || fluid[i].WR[0] <= 0)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].DerL[j] = 0;
            fluid[i].DerR[j] = 0;
            fluid[i].WL[j] = fluid[i].W[j];
            fluid[i].WR[j] = fluid[i].W[j];
        }
    }
}

void WENO5_AO(Block1d block, Fluid* fluid, Interface* interface)
{
    // w0, w_positive1, w_positive2, w_negative1, w_negative2
    // for characteristic variables, the base correspond to the prim variables
    // the W correspond to the conservative variables and w the characteristic variables

#pragma omp parallel  for
    for (int i = block.ghost - 1; i <= block.nodex + block.ghost; i++)
    {
        WENO_AO_cal(i, fluid, block);
    }

    // interface reconstruction
#pragma omp parallel  for
    for (int i = block.ghost; i <= block.nodex + block.ghost; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            interface[i].WL[j] = fluid[i - 1].WR[j];
            interface[i].WR[j] = fluid[i].WL[j];
            interface[i].DerL[j] = fluid[i - 1].DerR[j];
            interface[i].DerR[j] = fluid[i].DerL[j];
        }
    }
}

void Polynomial_3rd_left(double& var, double& der1, double w0, double wp1, double wn1, double alpha, double h)
{
    // Zero-mean form
    // p(x) = w0 + b(x - x0) + c(x^2 - x1)
    // Integrate[x - x0]_{-\delta x}^0 = 0
    // Integrate[x^2 - x1]_{-\delta x}^0 = 0
    double x0 = -1.0 / 2.0, x1 = 1.0 / 3.0;
    double b = (wp1 - w0) * alpha;
    double c = (1.0 / 2.0 * wn1 - w0 + 1.0 / 2.0 * wp1) * alpha;
    
    var = w0 + b * (-1.0 - x0) + c * (1.0 - x1);
    der1 = (b - 2.0 * c) / h;
}

void Polynomial_3rd_right(double& var, double& der1, double w0, double wp1, double wn1, double alpha, double h)
{
    // Zero-mean form
    // p(x) = w0 + b(x - x0) + c(x^2 - x1)
    // Integrate[x - x0]_{-\delta x}^0 = 0
    // Integrate[x^2 - x1]_{-\delta x}^0 = 0
    double x0 = -1.0 / 2.0, x1 = 1.0 / 3.0;
    double b = (wp1 - w0) * alpha;
    double c = (1.0 / 2.0 * wn1 - w0 + 1.0 / 2.0 * wp1) * alpha;
    
    var = w0 + b * (-x0) + c * (-x1);
    der1 = b / h;
}

void Polynomial_3rd_cal(int i, Fluid* fluid, Block1d block)
{
    double w[3], wp1[3], wn1[3], var[3]{}, der1[3]{};
    double base_left[3], base_right[3];
    double W[3], Wp1[3], Wn1[3];
    double ConW[3]{}, DerConW[3]{};
    for (int j = 0; j < 3; j++)
    {
        W[j] = fluid[i].W[j];
        Wp1[j] = fluid[i + 1].W[j];
        Wn1[j] = fluid[i - 1].W[j];

        w[j] = fluid[i].prim[j];
        wp1[j] = fluid[i + 1].prim[j];
        wn1[j] = fluid[i - 1].prim[j];

        base_left[j] = 0.5 * (wn1[j] + w[j]);
        base_right[j] = 0.5 * (w[j] + wp1[j]);
    }

    if (block.variable_type == "conservative")
    {
        // cell left side -- interface right side
        for (int j = 0; j < 3; j++)
        {
            Polynomial_3rd_left(fluid[i].WL[j], fluid[i].DerL[j], W[j], Wp1[j], Wn1[j], fluid[i].alpha, block.dx);
        }
        
        // cell right side -- interface left side
        for (int j = 0; j < 3; j++)
        {
            Polynomial_3rd_right(fluid[i].WR[j], fluid[i].DerR[j], W[j], Wp1[j], Wn1[j], fluid[i].alpha, block.dx);
        }
    }
    if (block.variable_type == "characteristic")
    {
        // cell left side -- interface right side
        Convar_to_Char(wn1, base_left, Wn1);
        Convar_to_Char(w, base_left, W);
        Convar_to_Char(wp1, base_left, Wp1);

        for (int j = 0; j < 3; j++)
        {
            Polynomial_3rd_left(var[j], der1[j], w[j], wp1[j], wn1[j], fluid[i].alpha, block.dx);
        }
        Char_to_Convar(ConW, base_left, var);
        Char_to_Convar(DerConW, base_left, der1);
        for (int j = 0; j < 3; j++)
        {
            fluid[i].WL[j] = ConW[j];
            fluid[i].DerL[j] = DerConW[j];
        }
        // cell right side -- interface left side
        Convar_to_Char(wn1, base_right, Wn1);
        Convar_to_Char(w, base_right, W);
        Convar_to_Char(wp1, base_right, Wp1);

        for (int j = 0; j < 3; j++)
        {
            Polynomial_3rd_right(var[j], der1[j], w[j], wp1[j], wn1[j], fluid[i].alpha, block.dx);
        }
        Char_to_Convar(ConW, base_right, var);
        Char_to_Convar(DerConW, base_right, der1);
        for (int j = 0; j < 3; j++)
        {
            fluid[i].WR[j] = ConW[j];
            fluid[i].DerR[j] = DerConW[j];
        }
    }

    // check if order reduce
    // if lambda < 0, then the order reduce to the first order
    double flag_l = fluid[i].WL[2] - 0.5 * fluid[i].WL[1] * fluid[i].WL[1] / fluid[i].WL[0];
    double flag_r = fluid[i].WR[2] - 0.5 * fluid[i].WR[1] * fluid[i].WR[1] / fluid[i].WR[0];
    if (flag_l <= 0 || flag_r <= 0 || fluid[i].WL[0] <= 0 || fluid[i].WR[0] <= 0)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].DerL[j] = 0;
            fluid[i].DerR[j] = 0;
            fluid[i].WL[j] = fluid[i].W[j];
            fluid[i].WR[j] = fluid[i].W[j];
        }
    }
}

void Polynomial_3rd(Block1d block, Fluid* fluid, Interface* interface)
{
    // w0, w_positive1, w_positive2, w_negative1, w_negative2
    // for characteristic variables, the base correspond to the prim variables
    // the W correspond to the conservative variables and w the characteristic variables

#pragma omp parallel  for
    for (int i = block.ghost - 1; i <= block.nodex + block.ghost; i++)
    {
        Polynomial_3rd_cal(i, fluid, block);
    }

    // interface reconstruction
#pragma omp parallel  for
    for (int i = block.ghost; i <= block.nodex + block.ghost; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            interface[i].WL[j] = fluid[i - 1].WR[j];
            interface[i].WR[j] = fluid[i].WL[j];
            interface[i].DerL[j] = fluid[i - 1].DerR[j];
            interface[i].DerR[j] = fluid[i].DerL[j];
        }
    }
}

void MMDF1d_calculate(MMDF1d& m, double* prim)
{
    double u = prim[1], lambda = prim[2];
    m.uwhole[0] = 1;
    m.uwhole[1] = u;
    m.uplus[0] = 0.5 * Alpha1d(lambda, -u);
    m.uminus[0] = 0.5 * Alpha1d(lambda, u);
    m.uplus[1] = u * m.uplus[0] + 0.5 * Beta1d(lambda, u);
    m.uminus[1] = u * m.uminus[0] - 0.5 * Beta1d(lambda, u);

    for (int i = 2; i <= 9; i++)
    {
        m.uwhole[i] = u * m.uwhole[i - 1] + 0.5 * (i - 1) / lambda * m.uwhole[i - 2];
    }
    for (int i = 2; i <= 9; i++)
    {
        m.uplus[i] = u * m.uplus[i - 1] + 0.5 * (i - 1) / lambda * m.uplus[i - 2];
        m.uminus[i] = u * m.uminus[i - 1] + 0.5 * (i - 1) / lambda * m.uminus[i - 2];
    }
    m.xi2 = 0.5 * K / lambda;
    m.xi4 = 0.25 * (K * K + 2 * K) / (lambda * lambda);

    for (int i = 0; i < 10; i++)
    {
        for (int k = 0; k < 4; k++)
        {
            if ((i + 2 * k) <= 9)
            {
                if (k == 0)
                {
                    m.upxi[i][k] = m.uplus[i];
                    m.unxi[i][k] = m.uminus[i];
                }
                if (k == 1)
                {
                    m.upxi[i][k] = m.uplus[i] * m.xi2;
                    m.unxi[i][k] = m.uminus[i] * m.xi2;
                }
                if (k == 2)
                {
                    m.upxi[i][k] = m.uplus[i] * m.xi4;
                    m.unxi[i][k] = m.uminus[i] * m.xi4;
                }
            }
        }
    }
    for (int i = 0; i < 10; i++)
    {
        for (int k = 0; k < 4; k++)
        {
            if ((i + 2 * k) <= 9)
            {
                if (k == 0)
                {
                    m.uxi[i][k] = m.uwhole[i];
                }
                if (k == 1)
                {
                    m.uxi[i][k] = m.uwhole[i] * m.xi2;
                }
                if (k == 2)
                {
                    m.uxi[i][k] = m.uwhole[i] * m.xi4;
                }
            }
        }
    }
}

void Microslope(double* a, double der[3], double prim[3])
{
    double R4, R2;
    R4 = der[2] / prim[0] - 0.5 * (prim[1] * prim[1] + 0.5 * (K + 1) / prim[2]) * der[0] / prim[0];
    R2 = (der[1] - prim[1] * der[0]) / prim[0];
    a[2] = 4 * prim[2] * prim[2] / (K + 1) * (2 * R4 - 2 * prim[1] * R2);
    a[1] = 2 * prim[2] * R2 - prim[1] * a[2];
    a[0] = der[0] / prim[0] - prim[1] * a[1] - 0.5 * a[2] * (prim[1] * prim[1] + 0.5 * (K + 1) / prim[2]);
}

void Center_3rd(Interface& interface, Fluid* fluid, Block1d block)
{
    double w[3]{}, wp[3]{};
    double dx = block.dx;
    for (int i = 0; i < 3; i++)
    {
        w[i] = fluid[0].W[i];
        wp[i] = fluid[1].W[i];
    }

    double convar_left[3]{}, convar_right[3]{};
    for (int i = 0; i < 3; i++)
    {
        convar_left[i] = interface.WL[i];
        convar_right[i] = interface.WR[i];
    }

    double prim_left[3], prim_right[3];
    Convar_to_ULambda_1d(prim_left, convar_left);
    Convar_to_ULambda_1d(prim_right, convar_right);

    MMDF1d ml, mr;
    MMDF1d_calculate(ml, prim_left);
    MMDF1d_calculate(mr, prim_right);

    double unit[3]{ 1.0, 0.0, 0.0 };

    double gl[3], gr[3];
    GL1d(0, 0, gl, unit, ml); // gl, means Wl(u>0), by input uint
    GR1d(0, 0, gr, unit, mr); // gr, means Wr(u<0), by input uint
    for (int i = 0; i < 3; i++)
    {
        interface.WC[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
    }
    for (int i = 0; i < 3; i++)
    {
        interface.DerC[i] = (wp[i] - w[i]) / dx;

    }
}

void Center_collision(Interface& interface, Fluid* fluid, Block1d block)
{
    double convar_left[3], convar_right[3];
    double prim_left[3], prim_right[3]; //rho, U, lambda
    for (int i = 0; i < 3; i++)
    {
        convar_left[i] = interface.WL[i];
        convar_right[i] = interface.WR[i];
    }
    Convar_to_ULambda_1d(prim_left, convar_left);
    Convar_to_ULambda_1d(prim_right, convar_right);

    MMDF1d ml, mr;
    MMDF1d_calculate(ml, prim_left);
    MMDF1d_calculate(mr, prim_right);

    double unit[3]{ 1.0, 0.0, 0.0 };

    double gl[3], gr[3];
    GL1d(0, 0, gl, unit, ml); // gl, means Wl(u>0), by input uint
    GR1d(0, 0, gr, unit, mr); // gr, means Wr(u<0), by input uint
    double axl[3], axr[3];
    Microslope(axl, interface.DerL, prim_left); // axl, means a coefficient indicating slope
    Microslope(axr, interface.DerR, prim_right); // axr, means a coefficient indicating slope
    double ax0l[3], ax0r[3];
    GL1d(0, 0, ax0l, axl, ml);  // ax0l, means Wlx(u>0), by input axl
    GR1d(0, 0, ax0r, axr, mr); // ax0r, means Wrx(u<0), by input axr
    for (int i = 0; i < 3; i++)
    {
        interface.WC[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
        interface.DerC[i] = convar_left[0] * ax0l[i] + convar_right[0] * ax0r[i];
    }
}

void Reconstruction_forg0(Interface* interface, Fluid* fluid, Block1d block, void (*g0reconstruction)(Interface&, Fluid*, Block1d))
{

    for (int i = block.ghost; i < block.nx - block.ghost + 1; i++)
    {
        (*g0reconstruction)(interface[i], &fluid[i-1], block);
    }
}

double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt, Block1d block)
{
    
    double C = block.c2 * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right);
    return block.c1 * dt + dt * C;
}

void flux_function(Interface& interface, Block1d block, double dt, int stage)
{
    double Flux[2][3];
    //change conservative variables to rho u lambda
    double convar_left[3], convar_right[3], convar0[3];
    for (int i = 0; i < 3; i++)
    {
        convar_left[i] = interface.WL[i];
        convar_right[i] = interface.WR[i];
        convar0[i] = interface.WC[i];
    }

    double prim_left[3], prim_right[3], prim0[3];
    Convar_to_ULambda_1d(prim_left, convar_left);
    Convar_to_ULambda_1d(prim_right, convar_right);
    Convar_to_ULambda_1d(prim0, convar0);

    double tau, tau_num;
    tau = 0;
    tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[2], prim_right[2], prim0[2], dt, block);
    double eta = exp(-dt / tau_num);
    double t[10];
    // non equ part time coefficient for gks_2nd algorithm (f0)
    t[0] = tau_num * (1 - eta); // this refers glu, gru part
    t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
    t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
    // then, equ part time coefficient for gks 2nd (g0)
    t[3] = tau_num * eta + dt - tau_num; //this refers g0u part
    t[4] = tau_num * (tau_num - eta * (dt + tau_num) - tau * (eta - 1)) - dt * tau; //this refers a0uu part
    t[5] = 0.5 * dt * dt - tau * tau_num * (eta - 1) - tau * dt; //this refers A0u part

    MMDF1d ml, mr, m0;
    
    MMDF1d_calculate(ml, prim_left);
    MMDF1d_calculate(mr, prim_right);
    MMDF1d_calculate(m0, prim0);

    double unit[3] = { 1, 0.0, 0.0 };

    double glu[3], gru[3], g0u[3];
    GL1d(1, 0, glu, unit, ml);
    GR1d(1, 0, gru, unit, mr);
    G1d(1, 0, g0u, unit, m0);

    //only one part, the kfvs1st part
    for (int i = 0; i < 3; i++)
    {
        Flux[0][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
    }

    //the equ g0u part, the gks1st part
    for (int i = 0; i < 3; i++)
    {
        Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
    }

    //for kfvs2nd part
    double der1left[3], der1right[3];
    for (int i = 0; i < 3; i++)
    {
        der1left[i] = interface.DerL[i];
        der1right[i] = interface.DerR[i];
    }

    double alx[3];
    Microslope(alx, der1left, prim_left);

    double alxuul[3];
    GL1d(2, 0, alxuul, alx, ml);

    double arx[3];
    Microslope(arx, der1right, prim_right);
    double arxuur[3];
    GR1d(2, 0, arxuur, arx, mr);
    for (int i = 0; i < 3; i++)
    {	// t1 part
        Flux[0][i] = Flux[0][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
    }

    // then we still need t[2], t[4] t[5] part for gks 2nd
    //for t[2] Aru,Alu part
    double alxu[3];
    double arxu[3];

    //take <u> moment for al, ar
    G1d(1, 0, alxu, alx, ml);
    G1d(1, 0, arxu, arx, mr);

    double Al[3], Ar[3];
    double der_AL[3], der_AR[3];

    //using compatability condition to get the time derivative
    for (int i = 0; i < 3; i++)
    {
        der_AL[i] = -prim_left[0] * (alxu[i]);
        der_AR[i] = -prim_right[0] * (arxu[i]);
    }
    // solve the coefficient martix b=ma
    Microslope(Al, der_AL, prim_left);
    Microslope(Ar, der_AR, prim_right);

    //to obtain the Alu and Aru
    double Alul[3];
    double Arur[3];
    GL1d(1, 0, Alul, Al, ml);
    GR1d(1, 0, Arur, Ar, mr);

    for (int i = 0; i < 3; i++)
    {	// t2 part
        Flux[0][i] = Flux[0][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
    }

    // for t[4] a0xuu part

    double a0x[3];
    double der1[3];

    for (int i = 0; i < 3; i++)
    {
        der1[i] = interface.DerC[i];
    }

    //solve the microslope
    Microslope(a0x, der1, prim0);
    //a0x <u> moment
    double a0xu[3];
    G1d(1, 0, a0xu, a0x, m0);
    //a0x <u^2> moment
    double a0xuu[3];
    G1d(2, 0, a0xuu, a0x, m0);

    for (int i = 0; i < 3; i++)
    {	// t4 part
        Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i]);
    }

    // for t[5] A0u part
    double derA0[3];

    for (int i = 0; i < 3; i++)
    {
        derA0[i] = -prim0[0] * (a0xu[i]);
    }
    double A0[3];
    Microslope(A0, derA0, prim0);
    double A0u[3];
    G1d(1, 0, A0u, A0, m0);
    for (int i = 0; i < 3; i++)
    {	// t5 part
        Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
    }

    // half time part
    double dt2 = 0.5 * dt;
    eta = exp(-dt2 / tau_num);
    // non equ part time coefficient for gks_2nd algorithm
    t[0] = tau_num * (1 - eta); // this refers glu, gru part
    t[1] = tau_num * (eta * (dt2 + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
    t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
    // then, equ part time coefficient for gks 2nd
    t[3] = tau_num * eta + dt2 - tau_num; //this refers g0u part
    t[4] = tau_num * (tau_num - eta * (dt2 + tau_num) - tau * (eta - 1)) - dt2 * tau; //this refers a0uu part
    t[5] = 0.5 * dt2 * dt2 - tau * tau_num * (eta - 1) - tau * dt2; //this refers A0u part

    for (int i = 0; i < 3; i++)
    {
        // t0 part
        Flux[1][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
        // t1 part
        Flux[1][i] = Flux[1][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
        // t2 part
        Flux[1][i] = Flux[1][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
        // t3 part
        Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
        // t4 part
        Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i]);
        // t5 part
        Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
    }

    for (int i = 0; i < 3; i++)
    {
        interface.fluxes[stage][i] = (4.0 * Flux[1][i] - Flux[0][i]);
        interface.der1fluxes[stage][i] = 4.0 * (Flux[0][i] - 2.0 * Flux[1][i]);
    }
}

void Calculate_Flux(Interface* interface, Block1d block, int stage)
{
    for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
    {
        flux_function(interface[i], block, block.dt, stage);
    }
}

void Update(Fluid* fluid, Interface* interface, Block1d block, int stage)
{
//#pragma omp parallel  for
    for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            double Flux = 0.0;
            for (int k = 0; k < stage + 1; ++k)
            {
                Flux = Flux
                    + block.timecoefficient[stage][k][0] * interface[i].fluxes[k][j]
                    + block.timecoefficient[stage][k][1] * interface[i].der1fluxes[k][j];
            }
            interface[i].Fluxes[stage][j] = Flux;
        }
    }
//#pragma omp parallel  for
    for (int i = block.ghost; i < block.nodex + block.ghost; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            fluid[i].W[j] = fluid[i].W_old[j] + 1.0 / block.dx * (interface[i].Fluxes[stage][j] - interface[i+1].Fluxes[stage][j]);
        }
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }

}

void free_boundary(Fluid* fluid, Block1d block)
{
    // free boundary left
    for (int i = block.ghost - 1; i >= 0; i--)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W[j] = fluid[i + 1].W[j];
        }
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }

    // free boundary right
    for (int i = block.nodex + block.ghost; i < block.nx; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W[j] = fluid[i - 1].W[j];
        }
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }
}

void free_boundary_left(Fluid* fluid, Block1d block)
{
    // free boundary left
    for (int i = block.ghost - 1; i >= 0; i--)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W[j] = fluid[i + 1].W[j];
        }
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }
}

void periodic_boundary(Fluid* fluid, Block1d block)
{
    // left boundary
    for (int i = 0; i < block.ghost; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W[j] = fluid[i + block.nodex].W[j];
        }
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }

    // right boundary
    for (int i = block.nodex + block.ghost; i < block.nx; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W[j] = fluid[i - block.nodex].W[j];
        }
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }
}

void reflect_boundary(Fluid* fluid, Block1d block)
{
    // left boundary
    for (int i = block.ghost - 1; i >= 0; i--)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W[j] = fluid[2 * block.ghost - 1 - i].W[j];
        }
        fluid[i].W[1] = -fluid[i].W[1];
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }

    // right boundary
    for (int i = block.nodex + block.ghost; i < block.nx; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W[j] = fluid[2 * (block.nodex + block.ghost) - 1 - i].W[j];
        }
        fluid[i].W[1] = -fluid[i].W[1];
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }
}

double Dtx(Block1d block, double dt, double R, double RU, double RE)
{
    double prim[3];
    prim[0] = R;
    prim[1] = RU / R;
    prim[2] = 0.4 * (RE - 0.5 * RU * RU / R);

    double tmp = abs(prim[1]) + sqrt(1.4 * prim[2] / prim[0]);
    // Euler flow
    if (tmp > block.CFL * block.dx / dt)
    {
        dt = block.CFL * block.dx / tmp;
    }
    return dt;
}

double Get_CFL(Fluid* fluid, Block1d block)
{
    double dt, R, RU, RE;
    dt = block.dx;
    for (int i = block.ghost; i <= block.nodex + block.ghost - 1; i++)
    {
        R = fluid[i].W[0];
        RU = fluid[i].W[1];
        RE = fluid[i].W[2];
        dt = Dtx(block, dt, R, RU, RE);
    }
    if (block.tnow + dt > block.tstop) {
        dt = block.tstop - block.tnow + 1e-20;
    }

    return dt;
}

void Copy_Array(Fluid* fluid, Block1d block)
{
    for (int i = 0; i < block.nx; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            fluid[i].W_old[j] = fluid[i].W[j];
        }
    }
}

void S1O2(Block1d& block)
{
    block.stages = 1;
    block.timecoefficient[0][0][0] = 1.0;
    block.timecoefficient[0][0][1] = 0.5;
}

void S2O4(Block1d& block)
{
    block.stages = 2;
    block.timecoefficient[0][0][0] = 0.5;
    block.timecoefficient[0][0][1] = 1.0 / 8.0;
    block.timecoefficient[1][0][0] = 1.0;
    block.timecoefficient[1][1][0] = 0.0;
    block.timecoefficient[1][0][1] = 1.0 / 6.0;
    block.timecoefficient[1][1][1] = 1.0 / 3.0;
}

void output1d(Fluid* fluid, Block1d block)
{
    ofstream ResultFile;
    const char* filePath_R = "/Users/hongzhang/Desktop/data/R.txt";
    const char* filePath_U = "/Users/hongzhang/Desktop/data/U.txt";
    const char* filePath_p = "/Users/hongzhang/Desktop/data/p.txt";
    const char* filePath_alpha = "/Users/hongzhang/Desktop/data/alpha.txt";

    ResultFile.open(filePath_R);
    for (int i = block.ghost; i < block.nodex + block.ghost; i++)
    {
        ResultFile << fluid[i].prim[0] << endl;
    }
    ResultFile.close();

    ResultFile.open(filePath_U);
    for (int i = block.ghost; i < block.nodex + block.ghost; i++)
    {
        ResultFile << fluid[i].prim[1] << endl;
    }
    ResultFile.close();

    ResultFile.open(filePath_p);
    for (int i = block.ghost; i < block.nodex + block.ghost; i++)
    {
        ResultFile << fluid[i].prim[2] << endl;
    }
    ResultFile.close();

    ResultFile.open(filePath_alpha);
    for (int i = block.ghost; i < block.nodex + block.ghost; i++)
    {
        ResultFile << fluid[i].alpha << endl;
    }
    ResultFile.close();
}

void (*boundary_type)(Fluid*, Block1d);

void (*reconstruction_type)(Block1d block, Fluid*, Interface*);

void (*g0reconstruction)(Interface&, Fluid*, Block1d);

void (*timecoe_list)(Block1d&); 

void Accuracy_test_1d()
{
    Block1d block;
    block.nodex = 80;
    block.ghost = 3;
    block.CFL = 0.5;
    block.dx = 2.0 / block.nodex;
    block.tstop = 2.0;
    block.tnow = 0.0;
    block.c1 = 0.0;
    block.c2 = 0.0;
    block.variable_type = "conservative";
    block.nx = block.nodex + 2 * block.ghost;
    auto* R_exact = new double[block.nodex];
    auto* fluid = new Fluid[block.nx];
    auto* interface = new Interface[block.nx];

    // multi-stages multi-derivatives
    boundary_type = periodic_boundary;
    reconstruction_type = Polynomial_3rd;
    g0reconstruction = Center_3rd;
    timecoe_list = S2O4;
    // time coefficient
    (*timecoe_list)(block);

    // Initial variables----periodic condition
    for (int i = block.ghost; i <= block.nodex + block.ghost - 1; i++)
    {
        fluid[i].W[0] = 1 - 0.2 / pi / block.dx * (cos(pi * (i - block.ghost + 1) * block.dx)
            - cos(pi * (i - block.ghost) * block.dx));
        fluid[i].W[1] = fluid[i].W[0] * 1;
        fluid[i].W[2] = 2.5 * 1 + 0.5 * fluid[i].W[0] * 1 * 1;
        R_exact[i - block.ghost] = fluid[i].W[0];
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }

    while (block.tnow < block.tstop)
    {
        // copy the flow variables
        Copy_Array(fluid, block);

        // get the time increment
        block.dt = Get_CFL(fluid, block);
        for (int stage = 0; stage < block.stages; stage++)
        {
            // boundary condition
            (*boundary_type)(fluid, block);

            // reconstruction within cells
            (*reconstruction_type)(block, fluid, interface);

            // center reconstruction
            Reconstruction_forg0(interface, fluid, block, *g0reconstruction);

            // calculate the flux
            Calculate_Flux(interface, block, stage);

            // update variables
            Update(fluid, interface, block, stage);
        }
        //Update_alpha(interface, fluid, block);
        block.tnow += block.dt;
    }

    double L1 = 0, L2 = 0, LInf = 0;
    for (int i = block.ghost; i < block.ghost + block.nodex; i++)
    {
        L1 += abs(fluid[i].W[0] - R_exact[i - block.ghost]) / block.nodex;
        L2 += pow(fluid[i].W[0] - R_exact[i - block.ghost], 2) / block.nodex;
        LInf = (LInf > abs(fluid[i].W[0] - R_exact[i - block.ghost])) ? LInf : abs(fluid[i].W[0] - R_exact[i - block.ghost]);
    }
    cout << L1 << endl;
}

void Riemann_problem_1d()
{
    /*
     * 1. sod shock tube case: tstop = 0.14, c1 = 0.05, c2 = 1, tau = 0, free boundary;
     * [rho, rhoU, rhoE] = [1, 0, 2.5] if 0.0 <= x < 0.5
     * [rho, rhoU, rhoE] = [0.125, 0, 0.25] if 0.5 <= x <= 1.0
     * 2. sjogreen test case: tstop = 0.15, c1 = 0.05, c2 = 1, tau = 0; free boundary;
     * [rho, rhoU, rhoE] = [1, -2, 3] if 0.0 <= x < 0.5
     * [rho, rhoU, rhoE] = [1, 2, 3] if 0.5 <= x <= 1.0
     * 3. blast wave case: tstop = 0.038, c1 = 0.05, c2 = 1, tau = 0, reflect boundary;
     * [rho, rhoU, rhoE] = [1, 0, 2500] if 0.0 <= x < 0.1
     * [rho, rhoU, rhoE] = [1, 0, 0.025] if 0.1 <= x < 0.9
     * [rho, rhoU, rhoE] = [1, 0, 250] if 0.9 <= x <= 1.0
     * 4. Noh case: tstop = 0.5, c1 = 0.05, c2 = 1, tau = 0, free boundary
     * [rho, rhoU, rhoE] = [1, 1, p0] if 0.0 <= x < 0.5
     * [rho, rhoU, rhoE] = [1, -1, p0] if 0.5 <= x <= 1.0
     * 5. Lax Harten case: tstop = 0.1, c1 = 0.05, c2 = 1, tau = 0, free boundary
     * [rho, rhoU, rhoE] = [0.445, 0.31061, 8.93] if 0.0 <= x < 0.5
     * [rho, rhoU, rhoE] = [0.5, 0, 1.4275] if 0.5 <= x <= 1.0
     * 6. Shu Osher case: tstop = 1.8, c1 = 0.05, c2 = 1, tau = 0, left free boundary, right exact boundary
     * [rho, U, p] = [3.857134, 2.629369, 10.33333] if 0.0 <= x < 1
     * [rho, U, p] = [1-(0.2/(5*dx))*(cos(5*(i-2)*dx)-cos(5*(i-3)*dx)), 0, 1]
     * if 1 <= x <= 10.0 + glost_cell * dx ==> exact boundary
     * flag = 0 correspond to S1O1
     * flag = 1 correspond to S1O2 & S2O4
     * variable_type = 0 correspond to characteristic variable
     * variable_type = 1 correspond to conservative variable
     * #pragma omp parallel  for
     * */
    Block1d block;
    block.nodex = 100;
    block.ghost = 3;
    block.CFL = 0.5;
    block.dx = 1.0 / block.nodex;
    block.tstop = 0.14;
    block.tnow = 0.0;
    block.c1 = 0.01;
    block.c2 = 5.0;
    block.variable_type = "characteristic";
    block.nx = block.nodex + 2 * block.ghost;

    auto* fluid = new Fluid[block.nx];
    auto* interface = new Interface[block.nx];

    // multi-stages multi-derivatives
    boundary_type = free_boundary;
    reconstruction_type = Polynomial_3rd;
    g0reconstruction = Center_collision;
    timecoe_list = S1O2;
    // time coefficient
    (*timecoe_list)(block);

    // initial condition
#pragma omp parallel for
    for (int i = block.ghost; i <= block.nodex + block.ghost - 1; i++)
    {
        if ((double)(i - block.ghost + 1) * block.dx <= 0.5)
        {
            fluid[i].W[0] = 1.0; fluid[i].W[1] = -2.0; fluid[i].W[2] = 2.5 * 0.4 + 2.0;
        }
        else
        {
            fluid[i].W[0] = 1.0; fluid[i].W[1] = 2.0; fluid[i].W[2] = 2.5 * 0.4 + 2.0;
        }
        fluid[i].prim[0] = fluid[i].W[0];
        fluid[i].prim[1] = fluid[i].W[1] / fluid[i].W[0];
        fluid[i].prim[2] = 0.4 * (fluid[i].W[2] - 0.5 * fluid[i].W[0] * pow(fluid[i].prim[1], 2));
    }
    while (block.tnow < block.tstop)
    {
        // copy the flow variables
        Copy_Array(fluid, block);

        // get the time increment
        block.dt = Get_CFL(fluid, block);
        
        // boundary condition
        (*boundary_type)(fluid, block);
        
        for (int stage = 0; stage < block.stages; stage++)
        {
            // boundary condition
            (*boundary_type)(fluid, block);

            // reconstruction within cells
            (*reconstruction_type)(block, fluid, interface);
            
            // center reconstruction
            Reconstruction_forg0(interface, fluid, block, *g0reconstruction);

            // calculate the flux
            Calculate_Flux(interface, block, stage);

            // update variables
            Update(fluid, interface, block, stage);
        }
        block.tnow += block.dt;
        block.step += 1;
        if (block.step % 100 == 0)
        {
            cout << "Progress degree: " << block.tnow / block.tstop * 100.0 << "%" << endl;
        }
    }
    output1d(fluid, block);
}
