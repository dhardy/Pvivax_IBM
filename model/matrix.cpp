/////////////////////////////////////////////////////////////////////////////
///  Individual-based Plasmodium vivax transmission model.                ///
///  Matrix functions                                                     ///
///  Dr Michael White, Imperial College London, m.white08@imperial.ac.uk  ///
/////////////////////////////////////////////////////////////////////////////

#include "matrix.h"
#include <cmath>

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//       //                                                                                           //
//  3.1  //  LU decomposition of a matrix                                                             // 
//       //  Based on ludcmp.cpp from Numerical Recipes in C++                                        //
//       //                                                                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                    //
//  Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise       //
//  permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;    //
//  indx[1..n] is an output vector that records the row permutation effected by the partial           //
//  pivoting; d is output as ±1 depending on whether the number of row interchanges was even          //
//  or odd, respectively. This routine is used in combination with lubksb to solve linear equations   //
//  or invert a matrix                                                                                //
//                                                                                                    // 
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ludcmp(vector<vector<double>> &a, int n_dim, vector<int> &indx, double &d)
{
    const double TINY = 1.0e-20;
    int i, imax, j, k;
    double big, dum, sum, temp;


    vector<double> vv(n_dim);
    d = 1.0;
    for (i = 0; i<n_dim; i++) {
        big = 0.0;
        for (j = 0; j<n_dim; j++)
            if ((temp = fabs(a[i][j])) > big) big = temp;
        if (big == 0.0) throw("Singular matrix in routine ludcmp");
        vv[i] = 1.0 / big;
    }
    for (j = 0; j<n_dim; j++) {
        for (i = 0; i<j; i++) {
            sum = a[i][j];
            for (k = 0; k<i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i<n_dim; i++) {
            sum = a[i][j];
            for (k = 0; k<j; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 0; k<n_dim; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        if (j != n_dim - 1) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i<n_dim; i++) a[i][j] *= dum;
        }
    }
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//       //                                                      // 
//  3.2  //  Matrix back substitution                            //
//       //  Based on lubksb.cpp from Numerical Recipes in C++   //
//       //                                                      //
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

void lubksb(vector<vector<double>> &a, int n_dim, vector<int> &indx, vector<double> &b)
{
    int i, ii = 0, ip, j;
    double sum;


    for (i = 0; i<n_dim; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii != 0)
            for (j = ii - 1; j<i; j++) sum -= a[i][j] * b[j];
        else if (sum != 0.0)
            ii = i + 1;
        b[i] = sum;
    }
    for (i = n_dim - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j<n_dim; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//       //                                                      //
//  3.3  //  Matrix inversion and calculation of determinant.    //
//       //  Based on ludcmp.cpp from Numerical Recipes in C++   //
//       //                                                      //
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

void matrix_inv(vector<vector<double>> &a, int n, vector<vector<double>> &a_inv)
{
    vector<int> a_index(n);
    vector<double> col(n);
    double d;

    ludcmp(a, n, a_index, d);

    for (int j = 0; j<n; j++)
    {
        for (int i = 0; i<n; i++)
        {
            col[i] = 0.0;
        }
        col[j] = 1.0;

        lubksb(a, n, a_index, col);

        for (int i = 0; i<n; i++)
        {
            a_inv[i][j] = col[i];
        }
    }

}


///////////////////////////////////////////////////
///////////////////////////////////////////////////
//       //                                      // 
//  3.4  //  Matrix multiplication               //
//       //  Calculates inv(MM*)xx               //
//       //                                      //
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void inv_MM_bb(vector<vector<double>> &MM, vector<double> &bb, vector<double> &xx, int n_dim)
{
    ///////////////////////////////////////////////
    // 3.4.1. calculate inv(MM)

    vector<vector<double>> MM_inv;
    MM_inv.resize(n_dim);
    for (int k = 0; k < n_dim; k++)
    {
        MM_inv[k].resize(n_dim);
    }

    matrix_inv(MM, n_dim, MM_inv);


    ///////////////////////////////////////////////
    // 3.4.2. calculate xx = MM_inv*bb

    for (int i = 0; i<n_dim; i++)
    {
        xx[i] = 0.0;

        for (int j = 0; j<n_dim; j++)
        {
            xx[i] = xx[i] + MM_inv[i][j] * bb[j];
        }
    }

}


///////////////////////////////////////////////////
///////////////////////////////////////////////////
//       //                                      // 
//  3.5  //  Equilibrium matrix                  //
//       //                                      //
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void MM_ij(int i, int j, Params& theta, double r_age[], vector<vector<double>>& MM,
    vector<vector<double>> lam_eq, vector<vector<vector<double>>> phi_LM_eq,
    vector<vector<vector<double>>> phi_D_eq, vector<vector<vector<double>>> r_PCR_eq)
{
    //////////////////////////////////////////////
    // 4.5.1. Initialise matrix with zeros

    for (int i1 = 0; i1<(N_H_comp * (K_max + 1)); i1++)
    {
        for (int j1 = 0; j1<(N_H_comp * (K_max + 1)); j1++)
        {
            MM[i1][j1] = 0.0;
        }
    }


    //////////////////////////////////////////////
    // 4.5.2. Fill out non-zero elements

    for (int k1 = 0; k1<(K_max + 1); k1++)
    {
        for (int k2 = 0; k2<(K_max + 1); k2++)
        {
            MM[0 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = - lam_eq[i][j] * theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];
            MM[0 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + r_PCR_eq[i][j][k2] * theta.D_MAT[k1][k2];
            MM[0 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = + theta.r_P*theta.D_MAT[k1][k2];

            MM[1 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * (1.0 - phi_LM_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_LM_eq[i][j][k2])*theta.K_MAT[k1][k2];
            MM[1 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = - lam_eq[i][j] * theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2] - r_PCR_eq[i][j][k2] * theta.D_MAT[k1][k2]
                                                             + lam_eq[i][j] * (1.0 - phi_LM_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_LM_eq[i][j][k2])*theta.K_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];
            MM[1 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + theta.r_LM*theta.D_MAT[k1][k2];

            MM[2 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta.K_MAT[k1][k2];
            MM[2 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2] * (1.0 - phi_D_eq[i][j][k2])*theta.K_MAT[k1][k2];
            MM[2 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = - lam_eq[i][j] * theta.D_MAT[k1][k2] - theta.ff*theta.K_MAT[k1][k2] - theta.r_LM*theta.D_MAT[k1][k2]
                                                             + lam_eq[i][j] * (1.0 - phi_D_eq[i][j][k2])*theta.OD_MAT[k1][k2] + theta.ff*(1.0 - phi_D_eq[i][j][k2])*theta.K_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];
            MM[2 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = + theta.r_D*theta.D_MAT[k1][k2];

            MM[3 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta.treat_cov*theta.treat_eff)*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta.treat_cov*theta.treat_eff)*theta.K_MAT[k1][k2];
            MM[3 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta.treat_cov*theta.treat_eff)*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * (1.0 - theta.treat_cov*theta.treat_eff)*theta.K_MAT[k1][k2];
            MM[3 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_D_eq[i][j][k2] * (1.0 - theta.treat_cov*theta.treat_eff)*theta.OD_MAT[k1][k2] + theta.ff*phi_D_eq[i][j][k2] * (1.0 - theta.treat_cov*theta.treat_eff)*theta.K_MAT[k1][k2];
            MM[3 * (K_max + 1) + k1][3 * (K_max + 1) + k2] = - lam_eq[i][j] * theta.D_MAT[k1][k2] - theta.r_D*theta.D_MAT[k1][k2] + lam_eq[i][j] * theta.OD_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];

            MM[4 * (K_max + 1) + k1][0 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta.treat_cov*theta.treat_eff*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta.treat_cov*theta.treat_eff*theta.K_MAT[k1][k2];
            MM[4 * (K_max + 1) + k1][1 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta.treat_cov*theta.treat_eff*theta.OD_MAT[k1][k2] + theta.ff*phi_LM_eq[i][j][k2] * phi_D_eq[i][j][k2] * theta.treat_cov*theta.treat_eff*theta.K_MAT[k1][k2];
            MM[4 * (K_max + 1) + k1][2 * (K_max + 1) + k2] = + lam_eq[i][j] * phi_D_eq[i][j][k2] * theta.treat_cov*theta.treat_eff*theta.OD_MAT[k1][k2] + theta.ff*phi_D_eq[i][j][k2] * theta.treat_cov*theta.treat_eff*theta.K_MAT[k1][k2];
            MM[4 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = - lam_eq[i][j] * theta.D_MAT[k1][k2] - theta.r_T*theta.D_MAT[k1][k2] + lam_eq[i][j] * theta.OD_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];

            MM[5 * (K_max + 1) + k1][4 * (K_max + 1) + k2] = + theta.r_T*theta.D_MAT[k1][k2];
            MM[5 * (K_max + 1) + k1][5 * (K_max + 1) + k2] = - lam_eq[i][j] * theta.D_MAT[k1][k2] - theta.r_P*theta.D_MAT[k1][k2] + lam_eq[i][j] * theta.OD_MAT[k1][k2]
                                                             + theta.gamma_L*theta.L_MAT[k1][k2] - theta.mu_H*theta.D_MAT[k1][k2] - r_age[i] * theta.D_MAT[k1][k2];
        }
    }

}
