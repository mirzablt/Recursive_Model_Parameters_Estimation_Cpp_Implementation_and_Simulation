#include "ERLSLibrary.h"


// Alokacija vektora demenzija 1 x n. 

double * allocate_array (int n){
    
	double * pok = (double *) calloc(n, sizeof(double));

	if (pok == NULL){
		return NULL;
	}	
	return pok;
}


// Oslobadjanje alociranog memorijskog prostora.

void clear_matrix (double ** matrix, int m){

	if(matrix == NULL) return;

	for (int i = 0; i < m; i++){
		free(matrix[i]);
	}
	free(matrix);	
}


// Alokacija matrice dimenzija m x n.

double ** allocate_matrix (int m, int n){
    
 	double ** matrix = (double** )calloc(m, sizeof(double *));

 	if (matrix == NULL) return NULL;

 	for (int i = 0; i < m; i++){
 		matrix[i] = (double * )calloc(n, sizeof(double));
 		if (matrix[i] == NULL){
		   clear_matrix(matrix,i);
		   return NULL;
	    } 	
    }
    return matrix;
}


//Alokacija matrice cije vrste imaju rezlicit broj clanova.

double ** allocate_rough_matrix (int m, const int *n){
    
 	double ** matrix = (double** )calloc(m, sizeof(double *));
    
 	if (matrix == NULL) return NULL;
    
 	for (int i = 0; i < m; i++){
 		matrix[i] = (double * )calloc(n[i] + 2, sizeof(double));
 		if (matrix[i] == NULL){
		   clear_matrix(matrix,i);
		   return NULL;
	    } 	
    }
    return matrix;
}


// Kasnjenje niza odabiraka za d odabiraka.

double delay (double * array, double u, int d){
    
	double out;
	out = array[d];
	array[0] = u;
    
	for (int i = d; i >= 1; i--){
         array[i] = array[i-1];
    }
    return out;
}


/* Funkcija regressor_grls  vrsi punjenje regresora definisanog prema (2.60), 
   dinamicki alociranog, sa pokazivacem psi. Parametri funkcije su: 
     u - matrica u ciji redovi odredjuju kasnjenja pojedinih mjerenih ulaza MISO procesa;
     y - mjereni izlaz procesa; w_hat, eps_hat;    
     v_hat – rezidali za opci linearni model (2.51) odredjeni relacijama (2.57);
     na, *nb, nf, nc, nd, nu – brojevi odredjeni prema (2.51); 
     model – struktura modela: 1 = ARX, 2 = ARMAX, 3 = BJ, 4 = OE. 
   Izlaz funkcije je vektor kolina   deklarisan kao struktura. */

struct Variable regressor_grls (double * psi, double * u, double y,
				double w_hat, double eps_hat, double v_hat,
				int na, const int * nb, int nu, int nf, int nc, int nd, int model){

	struct Variable out;
    
	int snb = 0;
	for (int i = 0; i < nu; i++) {
		snb += nb[i]; 
	}
	switch (model){
		case 1:    //ARX model
			out = initialize_vector_column(na + snb);
            
			for (int i = na - 1; i >= 1; i--){
    	         psi[i] = psi[i - 1];
            }

            snb = 0;
            for (int j = 0; j < nu; j++){
	            for (int i = nb[j] - 1; i >= 1; i--){
    	             psi[na + snb + i] = psi[na + snb + i-1];    
    	             snb += nb[j];
    	        }        
            }
            psi[0]            = -y;

            snb = 0;
            for (int i = 0; i < nu; i++){
       	         psi[na + snb]  =  u[i];
       	         snb += nb[i];
       	     }
       	    for (int i = 0; i < na + snb; i++){
                 out.X[i][0] =  psi[i];
            }     
		break;
		case 2:     // ARMAX model
			out = initialize_vector_column(na + snb + nc);
            
			for (int i = na - 1; i >= 1; i--){
    	         psi[i] = psi[i - 1];
            }

            snb = 0;
            for (int j = 0; j < nu; j++){
	            for (int i = nb[j] - 1; i >= 1; i--){
    	             psi[na + snb + i] = psi[na + snb + i - 1];    
    	             snb += nb[j];
    	        }        
            }
            for (int i = nc - 1; i >= 1; i--){
    	         psi[na + snb + i] = psi[na + snb + i - 1];
            }
            psi[0]     = -y;

            snb = 0;
            for (int i = 0; i < nu; i++){
       	         psi[na + snb]  =  u[i];
       	         snb += nb[i];
       	     }
       	    psi[na + snb] =  eps_hat;
       	    for (int i = 0; i < na + snb + nc; i++){
                 out.X[i][0] = psi[i];
            }     
        break;
        case 3:     // BJ model 
        	out = initialize_vector_column(snb + nf + nc + nd);
            
            snb = 0;
            for (int j = 0; j < nu; j++){
	            for (int i = nb[j] - 1; i >= 1; i--){
    	             psi[snb + i] = psi[snb + i - 1];    
    	             snb += nb[j];
    	        }        
            }
	        for (int i = nf - 1; i >= 1; i--){
    	         psi[snb + i] = psi[snb + i - 1];    
            }
            for (int i = nc - 1; i >= 1; i--){
    	         psi[snb + nf + i] = psi[snb + nf + i - 1];
            }
            for (int i = nd - 1; i >= 1; i--){
    	         psi[snb + nf + nc + i] = psi[snb + nf + nc + i - 1];
            }
            snb = 0;
            for (int i = 0; i < nu; i++){
       	         psi[snb]  =  u[i];
       	         snb += nb[i];
       	     }

       	    psi[snb]           = -w_hat;
       	    psi[snb + nf]      =  eps_hat;
       	    psi[snb + nf + nc] = -v_hat;

       	    for (int i = 0; i < snb + nf + nc + nd; i++){
                 out.X[i][0] = psi[i];
            }     
        break;
        case 4:     // OE model
        	out = initialize_vector_column(snb + nf);

            snb = 0;
            for (int j = 0; j < nu; j++){
	            for (int i = nb[j] - 1; i >= 1; i--){
    	             psi[snb + i] = psi[snb + i - 1];    
    	             snb += nb[j];
    	        }        
            }
	        for (int i = nf-1; i >= 1; i--){
    	         psi[snb + i] = psi[snb + i - 1];    
            }

            snb = 0;
            for (int i = 0; i < nu; i++){
       	         psi[snb]  =  u[i];
       	         snb += nb[i];
       	    }
       	    psi[snb]       = -w_hat;
       	    for (int i = 0; i < snb + nf; i++){
                 out.X[i][0] = psi[i];
            }     
        break;
	}
    return out;
}


/* Funkcija regressor_assembly racuna komponente w_hat, eps_hat, v_hat regresora
   psi prema relacijama (2.57). */
 
 void regressor_assembly(struct Variable& psi, struct Variable& theta, double y,
                        double* w_hat, double* eps_hat, double* v_hat,
                        int na, const int* nb, int nu, int nf, int nc, int nd, int model) {
    
    if (model == 1) {
        
        *w_hat   = 0;
        *eps_hat = 0;
        *v_hat   = 0;

    } else if (model == 2) {
        
        *w_hat   = 0;
        *eps_hat = y - inner(psi, theta, 0, na);
        *v_hat   = 0;

    } else if (model == 3) {
        
        *w_hat   = inner(psi, theta, 0, *nb) + inner(psi, theta, *nb, *nb + nf);
        *v_hat   = -*w_hat;
        *eps_hat = *v_hat - inner(psi, theta, *nb + nf + nc, *nb + nf + nc + nd) -
                            inner(psi, theta, *nb + nf, *nb + nf + nc);

    } else if (model == 4) {
        
        *w_hat   = inner(psi, theta, 0, *nb + nf);
        *eps_hat = 0;
        *v_hat   = 0;
    }
}
 
 
/* Funkcija regressor_flow pohranjuje regresore psi_k1 iz prethodnih od k-N+1 do k
   trenutaka odabiranja u matricu sa pokazivacem regressor_matrix. Sa obzirom na 
   regresore, ova  matrica cini FIFO listu, cime je njen izlaz za N odabiraka 
   zakasnjeli  regresor. Ova funkcija se koristi u rekurzivnom LS algoritmu sa 
   prozorom (sa konacnom historijom).  */

struct Variable regressor_flow (double ** regressor_matrix, struct Variable psi_k1, int N) {
	struct Variable out = initialize_vector_column(psi_k1.X.size());

 	for (int i = N - 1; i >= 1; i--){
 		for (int j = 0; j < psi_k1.X.size(); j++){
 			regressor_matrix[i][j] = regressor_matrix[i - 1][j];
		 }
	}
    
 	for (int j = 0; j < psi_k1.X.size(); j++){
 		regressor_matrix[0][j] = psi_k1.X[j][0];
	}
    for (int j = 0; j < psi_k1.X.size(); j++){
    	 out.X[j][0] = regressor_matrix[N - 1][j];	 
    }
	return out;	  
 } 


/* Funkcija y0_flow implementira shift registar sa pokazivacem y0, koji vrsi
   kasnjenje signala y za N odabiraka. */

double y0_flow (double * y0, double y, int N){
    
	double out;
    
 	for (int i = N - 1; i >= 1; i--){
 		y0[i] = y0[i - 1];
	}
    
	y0[0] = y;
	out   = y0[N - 1];
    
	return out;
 }


/* Za proces sa vise ulaza, estimirani parametri su struktuirani u matricu. Izlaz
   bloka S-funkcije u kojoj su implementirani algoritmi estimacije parametara 
   sistema je sabirnica PSAU_BUS, koju cine A, B, F, C, D estimirani parametri.
   Funkcija B_params vrsi struktuiranje B estimiranih parametara u matricu B. */

void B_params (PSAU_BUS * Params, struct Variable theta, const int *nb, int nu, int na){
	for (int i = 0; i < nu; i++) Params->B[i] = 0;
	int k = nu;
	 int j = 0;   
	for ( ; ; ){
		int snb = 0;
		for (int i = 0; i < nu; i++){
			if (j < nb[i]){ 
			    Params->B[k] = theta.X[na + j + snb][0];
			    k++;
			} else {
				Params->B[k] = 0;
				k++;
			}
			snb += nb[i];
		}
		j++;
		int mx = nb[0];
		for (int l = 0; l<nu; l++){
			if (nb[l] > mx) mx = nb[l];
		}
		if (j >= mx) break;
	}
}



 
 
 


