#include <math.h>
#include "MatrixLibrary.h"

#define pi 3.141592653


const double h = 10e-5;

using namespace std;


// Zbog mogucnosti proslijedjivanja cijelog vektora/matrice funkciji deklarisana je opsta struktura 'Variable'.
struct Variable{
	            vector< vector<double> >X;
               };


// s - skalarni dio kvaterniona, v - vektorski dio kvaterniona.               
struct Quaternion {
            	double s;
            	double v[3];
                  };  


// Inicijalizacija matrice dimenzija (m x n). 
                  
struct Variable initialize_matrix (int m, int n){
    struct Variable matrix;
    
	for (int i = 0; i < m; i++){
		vector<double> row;
		for (int j = 0; j < n; j++){
			row.push_back(0.0);	
		}
		matrix.X.push_back(row);
	}
	return matrix;
}


// Inicijalizacija vektor-kolone dimenzija (n x 1).

struct Variable initialize_vector_column (int n){
   	struct Variable matrix;

	for (int i = 0; i < n; i++){
		vector<double> column;
		column.push_back(0.0);	
		matrix.X.push_back(column);
	}
	return matrix;
}


// Inicijalizacija vektora dimenzija (1 x n).

struct Variable initialize_vector_row (int n){
    struct Variable matrix;
	vector<double> row;

	for (int i = 0; i < n; i++){
		row.push_back(0.0);	
	}
	matrix.X.push_back(row);
	return matrix;
}


// Transponovanje matrice ili vektora.

struct Variable transpose (struct Variable A){
    struct Variable B;
	int i, j;
    
	if (A.X.size() > 1 && A.X[0].size() == 1){
        
		B = initialize_vector_row(A.X.size());
		for (i=0; i < A.X.size(); i++){
		    B.X[0][i] = A.X[i][0];
		 }
        
	} else if (A.X.size() > 1 && A.X[0].size() > 1){
        
		B = initialize_matrix(A.X.size(), A.X[0].size());
		for (i = 0; i < A.X.size(); i++){
			for (j = 0; j < A.X[0].size(); j++){
				B.X[j][i] = A.X[i][j];
			}
		}
        
	} else if (A.X.size() == 1 && A.X[0].size() > 1){
        
		B = initialize_vector_column(A.X[0].size());
		for (i = 0; i < A.X[0].size(); i++){
			B.X[i][0] = B.X[0][i];
		}
        
	}  else return A;
    
	return B;
}


/* Suma matrica A i B. Zbog prirode koda gdje je ova funkcija pozvana,
   nije ispitivana jednakost dimenzija matrica A i B. */

struct Variable add (struct Variable A, struct Variable B){
    struct Variable C = initialize_matrix(A.X.size(), A.X[0].size());
    
	for(int i = 0; i < A.X.size(); i++){
		for(int j = 0; j < A.X[0].size(); j++){
			C.X[i][j] = A.X[i][j] + B.X[i][j];
		}
	}
	return C;	
}


// Razlika matrica A i B.

struct Variable subtract (struct Variable A, struct Variable B){
    struct Variable C = initialize_matrix(A.X.size(), A.X[0].size());
    
	for(int i = 0; i < A.X.size(); i++){
		for(int j=0; j < A.X[0].size(); j++){
			C.X[i][j] = A.X[i][j] - B.X[i][j];
		}
	}
	return C;	
}


// Mnozenje vektora/matrice A skalarom.

struct Variable c (struct Variable A, double scalar){
    
	struct Variable B = initialize_matrix(A.X.size(), A.X[0].size());
    
	for (int i = 0; i < A.X.size(); i++){
		for (int j = 0; j < A.X[0].size(); j++){
			B.X[i][j] = scalar * A.X[i][j];
		}
	}
	return B;
}


// Jedinicna matrica.

struct Variable eye (int order){ 
    
  struct Variable A = initialize_matrix(order, order);
  
  for (int i = 0; i < order; i++){
  	   for (int j = 0; j < order; j++){
  		   if(i == j){
               A.X[i][j] = 1.0;
			}
           else {
               A.X[i][j] = 0.0;
			}
	    }
    }
	return A; 
}   


// Proizvod matrica A i B.

struct Variable multiply (struct Variable A, struct Variable B){
    
	 int i, j, k;
     
	 int m = A.X.size();
	 int n = B.X.size();
	 int p = B.X[0].size();
     
	struct Variable C = initialize_matrix(m, p);
      
	double s;
	for (i = 0; i < m; i++){
        
		for (k = 0; k < p; k++){
            
			s = 0;
			for (j = 0; j < n; j++){
			s += (A.X[i][j]) * (B.X[j][k]);
		}
		C.X[i][k] = s;
		}
	}
	return C;
}


// Submatrica matrice A koja se dobije izostavljanjem p-te vrste i q-te kolone od A.

struct Variable submatrix ( struct Variable A, int p, int q){
    
	int i = 0, j = 0;
	int n = A.X.size();
    
	struct Variable B = initialize_matrix(n-1, n-1);
    
	for(int row = 0; row < n; row++){
		for(int col = 0; col < n; col++){
            
			if(row != p && col != q){
				B.X[i][j++] = A.X[row][col];
                
				if(j == n - 1) {
                    j = 0; 
                    i++;
                }
			}
		}
	}
	return B;
}


// Rekurzivna fukcija koja racuna determinantu matrice A.

double determinant (struct Variable A, int n){
    
	double D = 0;
    
	if(n == 1)
	   return A.X[0][0]; 
    
	int sign = 1;
	for(int f = 0; f < n; f++) {
		D += sign * A.X[0][f] * determinant( submatrix(A, 0, f), n - 1);
		sign = -sign;
	}
	return D;
}


// Racunanje adjungirane  matrice od A.

struct Variable adjung (struct Variable A){
    
	int i, j;
	int n = A.X.size();
    
	struct Variable B = initialize_matrix(n, n);
    
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			B.X[j][i] = pow(-1, i + j) * determinant( submatrix(A, i, j), n - 1);
		}
	}
	return B;		
}


// Racunanje inverzne matrice od A.

struct Variable inv (struct Variable A){
    
	int n = A.X.size();
    
	struct Variable B;
	struct Variable C = initialize_matrix(n, n);
	B = adjung(A);
    
	if(determinant(A, n) != 0){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				C.X[i][j] = B.X[i][j] / determinant(A, n);
			}
		}
		return C;
   }
}


// Euklidska norma vektora 'vec'. 

double vector_norm (struct Variable vec){
    
	double norm = 0.0;
    
	for (int i = 0; i < vec.X.size(); i++)
	     norm += vec.X[i][0] * vec.X[i][0];
	return sqrt(norm);     
}


// Unutarnji-skalarni proizvod dva vektora.

double inner (struct Variable A, struct Variable B, int lb, int ub){

 	double s = 0;
 	for (int i = lb; i < ub; i++){
 		s = s + A.X[i][0] * B.X[i][0];
	}
	 return s;
 }


// Racunanje maksimalnog clana niza.

 double max_ (double * x, int N){
     
	double mx = x[0];

	for (int i = 1; i < N; i++){
		if (x[i] > mx)
		    mx = x[i];
	}
	return mx;
}
