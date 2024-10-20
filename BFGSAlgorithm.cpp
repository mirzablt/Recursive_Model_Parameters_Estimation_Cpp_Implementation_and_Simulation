
#include "BFGSAlgorithm.h"



/* Izlaz  funkcije regressor_window je matrica sa N vrsta i brojem kolona jednakom
   dimenziji vektora parametara (2.60), cije vste cine regresori od k-N+1 do k-tog
   (sadasnjeg) trenutka odabiranja.*/

struct Variable regressor_window (double ** regressor_matrix, struct Variable psi_k1, int N) {
    
	struct Variable out = initialize_matrix(N, psi_k1.X.size());

 	for (int i = N-1; i >= 1; i--){
 		for (int j = 0; j < psi_k1.X.size(); j++){
 			regressor_matrix[i][j] = regressor_matrix[i - 1][j];
		 }
	}
    
 	for (int j = 0; j < psi_k1.X.size(); j++){
 		regressor_matrix[0][j] = psi_k1.X[j][0];
	}
    
	for (int i = 0; i < N; i++){
        for (int j = 0; j < psi_k1.X.size(); j++){
    	     out.X[i][j] = regressor_matrix[i][j];	 
        }
    }
	return out;	  
 }  


/* Izlaz funkcije y0_window jeste vektor ciji su elementi jednaki izmjerenim 
   izlazima procesa od trenutka k-N+1 do k-tog (tekuceg) trenutka iteracije. */

struct Variable y0_window (double * y0, double y, int N){
    
	struct Variable out = initialize_vector_column(N);

 	for (int i = N-1; i >= 1; i--){
 		y0[i] = y0[i-1];
	}
	y0[0] = y;
    
	for (int i = 0; i < N; i++){
		out.X[i][0] = y0[i];
	}
	return out;
 }


/* U funkciji criterion je implementirana fukcija cilja (A.9) koja izrazava mjeru
   odstupanja predikcije modela od izmjerene vrijednosti izlaza sistema y0. Ovdje
   se uzimaju vrijednosti izlaza od trenutka k-N+1 do k-tog (tekuceg) trenutka.
   Vrijednosti izlaza modela dobiju se mnozenjem vektora parametara theta i regresora. */

double criterion (struct Variable regressor_mat, struct Variable y0, struct Variable theta){
     
	struct Variable s = subtract(y0, multiply(regressor_mat, theta));

	double d = 0;
	for (int i = 0; i < regressor_mat.X.size(); i++){
		d += pow(s.X[i][0], 2);
	}
	return d;
}
 
 
 /* Naredne linije koda su implementirane u prethodnom projektu u kome je impementirana
    inverzna kinematika robota, pa se u komentarima pozivalo na reference u .pdf-u koji
    je prilozen u repozitoriju tog projekta. */
 
 /* Funkcija gradient racuna gradijent funkcije f po promjenjivoj var. Kako racunamo
    gradijent funkcije cilja, to ovom metodu moramo proslijediti i parametre analogno
    kao u slucaju metoda u kojem je implementirana funk. cilja. */
 
struct Variable gradient (double(*f)(struct Variable, struct Variable, struct Variable),
                          struct Variable dh, struct Variable dpc, struct Variable var){
	
    struct Variable var_1, var_2;
    struct Variable partial_derivatives;
    partial_derivatives = initialize_vector_column(var.X.size());

	for(int i = 0; i < var.X.size(); i++){
		         var_1 = var;
	             var_2 = var;
		var_1.X[i][0] += h;
	    var_2.X[i][0] -= h;
		partial_derivatives.X[i][0] = ( f(dh, dpc, var_1) - f(dh, dpc, var_2) ) / (2.0 * h);
	}
	return partial_derivatives;	
}


/* Funcija line odredjuje tacke linije po kojoj vrsimo jednodimenzionalno pretrazivanje
   po parametru alpha. Ona je odredjena tackom point i vektorom pravca direction. */

struct Variable line (struct Variable point, struct Variable direction, double alpha){
    
	struct Variable linea=initialize_vector_column(point.X.size());

	for (int i=0; i<point.X.size(); i++){
		linea.X[i][0] = point.X[i][0] + alpha * direction.X[i][0];
	}
	return linea;
}	


/* Izvod funkcije f na pravoj odredjenom tackom point i pravcem direction, 
   za parametar t. */

double derivative_on_line (double( *f)(struct Variable, struct Variable, struct Variable),
						   struct Variable dh, struct Variable dpc, struct Variable point,
                           struct Variable direction, double t){

	double df = f(dh, dpc, line(point, direction, t + h)) -
                f(dh, dpc, line(point, direction, t - h));

	return df / (2.0 * h);
}


/* Drugi izvod funkcije f na pravoj odredjenom tackom point i pravcem direction, 
   za parametar t. */

double second_derivative (double( *f)(struct Variable, struct Variable, struct Variable),
						  struct Variable dh, struct Variable dpc, struct Variable point,
                          struct Variable direction, double t){

	double ddf = derivative_on_line (f, dh, dpc, point, direction, t + h) -
                 derivative_on_line (f, dh, dpc, point, direction, t-h);

	return ddf / (2.0 * h); 
} 


/* Jednodimenzionalno pretrazivanje za odredjivanje optimalnog koraka prema (2.16) 
   odnosno (2.61) odredjeno sa (2.24) po pravoj  kroz tacku point sa pravcem direction. 
   Kriterij zaustavljanja (eps) je dostizanje minimalne promjene problemske varijable.
   Pocetna vrijednost x1k je uzeta proizvoljno. */

struct Variable newton_two_points (double( *f)(struct Variable, struct Variable, struct Variable),
                                   struct Variable dh, struct Variable dpc,
								   struct Variable point, struct Variable direction){

	const double EPSILON =0.1;
    int i;
	double xk1;
	double eps;
    
	struct Variable XK = point;
	struct Variable optimal_step = initialize_vector_column(1);
    
	double xk = 0.0; 
	double x1k = 0.02; 

	do{
		double dxk = xk - x1k;
		double d_derivative_on_line_k = derivative_on_line( f, dh, dpc, XK, direction, xk) -
									    derivative_on_line( f, dh, dpc, XK, direction, x1k);

		xk1 = x1k + (
			  dxk / (
			  1 - ((derivative_on_line(f, dh, dpc, XK, direction, xk) /
			       derivative_on_line(f, dh, dpc, XK, direction, x1k)) *
				   (d_derivative_on_line_k /(dxk * second_derivative( f, dh, dpc, XK, direction, xk))))));

		eps = fabs(xk1-xk);
		x1k = xk;
		xk = xk1;

	}while (eps > EPSILON); 

	optimal_step.X[0][0] = xk1;
	return optimal_step;
}


/* Implementacija BFGS algoritma. Algoritam odredjuje vektor parametara za kojeg
   funkcija cilja f ima minimalnu vrijednost. U idealnom slucaju ta vrijednost je nula,
   sto odgovara slucaju da se je predikcija modela jednaka izmjerenoj vrijednosti.
   U funciji cilja predikcija modela (od k-N+1 do k-tog trenutka) se dobije na osnovu
   vektora parametara theta i regresora (parametar dh).
   A izmjerene vrijednosti izlaza se prislijedjuju preko parametra dp_htm.
   Xk je vektor problemskih varijabli (i predstavlja vektor parametara sistema theta). */

struct Variable BFGS ( double( *f)(struct Variable, struct Variable, struct Variable),
                       struct Variable dh, struct Variable dp_htm, struct Variable start_point, double tolerance){
	
    struct Variable Xk, Xk1, DXk;
	       Variable Gk, Gk1, DGk;
	       Variable step, Rk;
           
	       Variable Hk, Hk1;
	       Variable Mk, Nk;

	double gamma, beta;	
	Xk = start_point;
	Hk = eye(start_point.X.size());

	while (vector_norm(gradient(f, dh, dp_htm, Xk)) > tolerance){      // Kriterij zaustavljanja dat sa (2.59).
        
		Gk    = gradient(f, dh, dp_htm, Xk);
		Rk    = c(multiply(Hk, Gk), -1);                               // Pravac pretrazivanja dat sa (2.60).
        
		step  = newton_two_points (f, dh, dp_htm, Xk, Rk);             // Jednodimenzionalno pretrazivanje (2.61).
        
		Xk1   = add(Xk, multiply(Rk, step));                           // Racunanje nove aproksimacije (2.62).
		Gk1   = gradient(f, dh, dp_htm, Xk1);
        
		DXk   = subtract(Xk1, Xk);                                     // Razlika data sa(2.63a).
		DGk   = subtract(Gk1, Gk);                                     // Razlika data sa(2.63b).
        
		gamma = (multiply( transpose(DXk), DGk)).X[0][0];
		beta  = (multiply( multiply( transpose(DGk), Hk), DGk)).X[0][0];
        
		Mk    = c(multiply(DXk, transpose(DXk)), (gamma+beta)/pow(gamma,2));
		Nk    = c( add( multiply( multiply(Hk, DGk), transpose(DXk)), 
                   multiply( multiply(DXk, transpose(DGk)), Hk)), -1.0/gamma);
        
		Hk1   = add( add(Hk, Mk), Nk);                                 // Racunanje aproksimacije inverznog Hessiana H(k+1) prema (2.56).
        
		Hk    = Hk1;
		Xk    = Xk1;	
	}
	return Xk;
}




