/////////////////////////////////////////////////////////////////////////////
///  Individual-based Plasmodium vivax transmission model.                ///
///  Matrix functions                                                     ///
///  Dr Michael White, Imperial College London, m.white08@imperial.ac.uk  ///
/////////////////////////////////////////////////////////////////////////////

// Include guard
#ifndef PVIVAX_MODEL_MATRIX
#define PVIVAX_MODEL_MATRIX

#include "Params.h"

void ludcmp(vector<vector<double>>& a, int n_dim, vector<int>& indx, double& d);
void lubksb(vector<vector<double>>& a, int n_dim, vector<int>& indx, vector<double>& b);
void matrix_inv(vector<vector<double>>& a, int n, vector<vector<double>>& a_inv);
void inv_MM_bb(vector<vector<double>>& MM, vector<double>& bb, vector<double>& xx, int n_dim);
void MM_ij(int i, int j, Params& theta, double r_age[], vector<vector<double>>& MM,
           vector<vector<double>>& lam_eq, vector<vector<vector<double>>>& phi_LM_eq,
           vector<vector<vector<double>>>& phi_D_eq, vector<vector<vector<double>>>& r_PCR_eq);

#endif
