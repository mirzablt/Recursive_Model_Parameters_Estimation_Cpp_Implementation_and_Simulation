#ifndef ERLSLIBRARY_H
#define ERLSLIBRARY_H



// Deklaracija funkcija definisanih u ERLSLibrary.cpp

double * allocate_array (int n);

void clear_matrix (double ** matrix, int m);

double ** allocate_matrix (int m, int n);

double ** allocate_rough_matrix (int m, const int *n);


double delay (double * array, double u, int d);

struct Variable regressor_grl (double * psi, double u, double y,
                              double w_hat, double eps_hat, double v_hat,
                              int na, int nb, int nf, int nc, int nd, int model);

struct Variable regressor_grls (double * psi, double * u, double y,
                                double w_hat, double eps_hat, double v_hat,
                                int na, const int * nb, int nu, int nf, int nc, int nd, int model);

void regressor_assembly(struct Variable& psi, struct Variable& theta, double y,
                        double* w_hat, double* eps_hat, double* v_hat,
                        int na, const int* nb, int nu, int nf, int nc, int nd, int model);

struct Variable regressor_flow (double ** regressor_matrix, struct Variable psi_k1, int N);

double y0_flow (double * y0, double y, int N);


void B_params (PSAU_BUS * Params, struct Variable theta, const int *nb, int nu, int na);


#endif
