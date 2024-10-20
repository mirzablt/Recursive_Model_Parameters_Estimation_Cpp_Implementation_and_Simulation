#ifndef MATRIXLIBRARY_H
#define MATRIXLIBRARY_H


//Deklaracija funkcija definisanih u matrixlibrary.cpp

struct Variable initialize_matrix (int m, int n);

struct Variable initialize_vector_column (int n);

struct Variable initialize_vector_row (int n);


struct Variable transpose (struct Variable A);

struct Variable add (struct Variable A, struct Variable B);

struct Variable subtract (struct Variable A, struct Variable B);

struct Variable c (struct Variable A, double scalar);

struct Variable eye (int order);

struct Variable multiply (struct Variable A, struct Variable B);


struct Variable submatrix ( struct Variable A, int p, int q);

double determinant (struct Variable A, int n);

struct Variable adjung (struct Variable A);

struct Variable inv (struct Variable A);


double vector_norm (struct Variable vec);

double inner (struct Variable A, struct Variable B, int lb, int ub);

double max_ (double * x, int N);

#endif
