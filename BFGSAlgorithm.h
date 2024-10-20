#ifndef BFGSALGORITHM_H
#define BFGSALGORITHM_H



//Deklaracija funkcija definisanih u BFGSAlgorithm.cpp


struct Variable regressor_window (double ** regressor_matrix, struct Variable psi_k1, int N);

struct Variable y0_window (double * y0, double y, int N);


double criterion (struct Variable regressor_mat, struct Variable y0, struct Variable theta);

struct Variable gradient (double(*f)(struct Variable, struct Variable, struct Variable),
                          struct Variable dh, struct Variable dpc, struct Variable var);


struct Variable line (struct Variable point, struct Variable direction, double alpha);

double derivative_on_line (double( *f)(struct Variable, struct Variable, struct Variable),
                           struct Variable dh, struct Variable dpc, struct Variable point,
                           struct Variable direction, double t);

double second_derivative (double( *f)(struct Variable, struct Variable, struct Variable),
                          struct Variable dh, struct Variable dpc, struct Variable point,
                          struct Variable direction, double t);


struct Variable newton_two_points (double( *f)(struct Variable, struct Variable, struct Variable),
                                   struct Variable dh, struct Variable dpc,
                                   struct Variable point, struct Variable direction);


struct Variable BFGS ( double( *f)(struct Variable, struct Variable, struct Variable),
                       struct Variable dh, struct Variable dp_htm,
                       struct Variable start_point, double tolerance);






 #endif
