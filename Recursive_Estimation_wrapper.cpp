

#include "simstruc.h"

#include "Recursive_Estimation_bus.h"

#include <math.h>
#include <vector>
#include <stdlib.h>

#include "MatrixLibrary.cpp"
#include "ERLSLibrary.cpp"
#include "BFGSAlgorithm.cpp"

#define u_width 2
#define y_width 1



void Recursive_Estimation_Start_wrapper(void **pW,
			const int32_T *na, const int_T p_width0,
			const int32_T *nb, const int_T p_width1,
			const int32_T *nf, const int_T p_width2,
			const int32_T *nc, const int_T p_width3,
			const int32_T *nd, const int_T p_width4,
			const int32_T *d, const int_T p_width5,
			const real_T *lambda, const int_T p_width6,
			const real_T *alpha, const int_T p_width7,
			const int32_T *model, const int_T p_width8,
			const int32_T *method, const int_T p_width9,
			const int32_T *win_length, const int_T p_width10,
			const real_T *tolerance, const int_T p_width11,
			const real_T *sample_time, const int_T p_width12,
			const int32_T *history, const int_T p_width13,
			const real_T *pnc, const int_T p_width14,
			SimStruct *S)
{

/* Kod u narednim linijama je implementiran u bloku S-Function Builder u mapi Start.
   Sve linije koda koje prethode ovom komentaru (kao i header Recursive_Estimation_bus.h) 
   su automatski generirane pritiskom tastera Build u bloku S-Function Builder. 
   U mapi Libraries u bloku S-Function Builder su navedene implemetirane biblioteke:
        -ERLSLibrary, 
        -BFGSAlgorithm,
        -MatrixLibrary. 
   
   Parametri bloka S-funkcije u kojoj su implementirani algoritmi estimacije parametara 
   modela sistema, su:
     -na, nb, nc, nd, nf - red polonoma A, B, C, D, F koji karakterisu izabrani model;
     -d                  - kasnjenje izlaza modela sistema;
     -lambda             - faktor iscezavanja. Parametar koji karakerise ERLS metod sa 
                           eksponencijalnim faktorom iscezavanja;
     -alpha              - dijagonalni elementi pocetne vrijednosti matrice kovarijanse
                           procesnog suma. Ovaj param. karakterise ERLS metod;
     -model              - struktura modela sistema. 1-ARX, 2-ARMAX, 3-BJ, 4-OE;
     -method             - izbor metode (agoritma) estimacije: 1-ERLS sa faktorom iscezavanja,
                          2-BFGS, 3-Kalmanov filter;
     -win_length         - sirina vremenskog okvira. Ovaj parametar karakterise LS algoritam 
                           sa konacnim vrem. okv. (konacnom historijom) i  funkciju kriterija 
                           kod  primjene BFGS algoritma u estimaciji parametara modela;    
     -tolerance          - tacnost do koje BFGS algoritam racuna optimum funkcije krterija;
     -sample_time        - perioda uzorkovanja;
     -pnc                - vrijednosti dijagonalnih elemenata matrice kovarijanse procesnog
                           suma (kod metode KF).
   Ovi parametri su deklarisani u mapi: Data Properties/Parameters bloka S-funkcije.
 
   Ulazi bloka su:
     -u0     - vektor ulaza sistema (procesa) cije parametre modela estimiramo;
     -y0     - izlaz sistema cije parametre modela estimiramo;
     -Enable - ulaz kojim aktiviramo/deaktiviramo  estmaciju parametara sistema;
   i deklarisani u mapi Data Properties/Input ports bloka S-funkcije.
 
   Izlazi bloka su:
     -Params     - vektor estimiranih parametara modela sistema;
     -Error      - greksa predikcije modela;
     -Excitation - mjera pobudjenosti sistema cije parametre modela estimiramo. */


int DIM;             // Dimenzija vektora regresora
int nu = p_width1;   //Broj ulaza sistema.

// Ovisno od izbora modela sistema, imamo razlicite dimenzije DIM vektora regresora.
int snb = 0;
for (int i = 0; i < nu; i++){
    snb += nb[i];
}
if (*model == 1){                   // ARX model
   DIM = *na + snb;
} 
else if (*model == 2){              // ARMAX 
   DIM = *na + snb + *nc;
} 
else if (*model == 3){              // BJ 
   DIM = snb + *nf + *nc + *nd;
} 
else if (*model == 4){              // OE 
    DIM = snb + *nf;
} 

/* pW je niz pokazivaca i njihov broj se odabire u S-Function Bulderu, u mapi Initialization
   i u polju Number of PWorks. Ovaj niz pokazivača se koristi za pohranu vrijednosti varijabli
   izmedju različitih koraka simulacije. */

// Alokacija radnog vektora regresora.
double * psi = allocate_array(DIM);
pW[0] = psi;

// Alokacija radnog vektora (shift registra) ulaza u sistem kojeg identificiramo.
double ** inputs_delays = allocate_rough_matrix(nu, d);
pW[1] = inputs_delays;

// Alokacija radnog vektora (shift registra) izlaznog signala iz sistema kojeg identificiramo.
pW[2] = allocate_array(2);

// Alokacija matrice kovarijanse P.
double ** p = allocate_matrix(DIM, DIM);

// Inicijalizacija matrice P prema (2.17)
for (int i = 0; i < DIM; i++){
    p[i][i] = *alpha * *alpha;
}    
pW[3] = p;

// Alokacija radnog vektora  estimiranih parametara sistema.
double * theta = allocate_array(DIM);
pW[4] = theta;

if (*method == 3){
    for (int i = 0; i < DIM; i++){
        theta[i] = -0.1;
    }  
}

// Alokacija prostora za cuvanje vrijednosti residuala eps_hat prethodne iteracije.
pW[5] = allocate_array(1);

// Alokacija prostora za cuvanje vrijednosti residuala w_hat prethodne iteracije.
pW[6] = allocate_array(1);

// Alokacija prostora za cuvanje vrijednosti residuala v_hat prethodne iteracije.
pW[7] = allocate_array(1);


/* U slucaju izbora odgovarajuceg  modela, metode i duzine intervala, imamo alokaciju 
   matrice regresora regressor_mat. U slucaju izbora metoda sa konacnom historijom 
   win_lenght predstavlja sirinu intervala od sadasnjeg trenutka odabiranja k 
   do k - win_lemght + 1 proslog trenutka */ 

if ((*method==1 && *history==2 && (*model==1 || *model==2 || *model==4)) || 
    (*method==2 && (*model==1 || *model==2 || *model==3 || *model==4)))    {
    
    double ** regressor_mat = allocate_matrix(*win_length, DIM);
    
    for (int i=0; i<DIM; i++){
        regressor_mat[i][DIM-1-i] = 1 / (*alpha);
    } 
    
    pW[8] = regressor_mat;
    double * y0_vec = allocate_array(*win_length);
    pW[9] = y0_vec;
}

}


 
void Recursive_Estimation_Outputs_wrapper(const real_T *u0,
			const real_T *y0,
			const real_T *Enable,
			PSAU_BUS *Params,
			real_T *Error,
			real_T *Excitation,
			void **pW,
			const int32_T *na, const int_T p_width0,
			const int32_T *nb, const int_T p_width1,
			const int32_T *nf, const int_T p_width2,
			const int32_T *nc, const int_T p_width3,
			const int32_T *nd, const int_T p_width4,
			const int32_T *d, const int_T p_width5,
			const real_T *lambda, const int_T p_width6,
			const real_T *alpha, const int_T p_width7,
			const int32_T *model, const int_T p_width8,
			const int32_T *method, const int_T p_width9,
			const int32_T *win_length, const int_T p_width10,
			const real_T *tolerance, const int_T p_width11,
			const real_T *sample_time, const int_T p_width12,
			const int32_T *history, const int_T p_width13,
			const real_T *pnc, const int_T p_width14,
			SimStruct *S)
{
	
// Kod u narednim linijama je implementiran u bloku S-Function Builder u mapi Outputs.


    double y;
    double fade, pn_cov;
    int    DIM;
    int nu = p_width1;
    double u[nu];

    int snb = 0;
    for (int i=0; i<nu; i++){
        snb += nb[i];
    }

    // Ovisno od izbora modela sistema, imamo razlicite dimenzije DIM vektora regresora
    if (*model==1){                    // ARX
       DIM = *na + snb;
    } 
    else if (*model==2){             // ARMAX
       DIM = *na + snb + *nc;
    }
    else if (*model==3){             // BJ
       DIM = snb + *nf + *nc + *nd;
    }
    else if (*model==4){             // OE
        DIM = snb + *nf;
    } 

    struct Variable psi_k1;
    struct Variable residual_k1 = initialize_vector_column(1);
    struct Variable theta_hat_k = initialize_vector_column(DIM);
    struct Variable theta_hat_k1;
    struct Variable P_k         = initialize_matrix(DIM, DIM);
    struct Variable P_k1;
    struct Variable Pk_psik1, psik1tr_Pk_psik1;
    struct Variable gamma_k;
    struct Variable psi_kN;
    double y0kN;


    /* pW je niz pokazivaca i njihov broj se odabire u S-Function Bulderu, u mapi Initialization
       i u polju Number of PWorks. Ovaj niz pokazivača se koristi za pohranu vrijednosti varijabli
       izmedju različitih koraka simulacije. */

    // Deklaracije pokazivaca za prethodno alocirane vektore i matrice.
    double *  psi        = (double *) pW[0];
    double ** u_buf      = (double **)pW[1];
    double *  y_buf      = (double *) pW[2];
    double ** P_buf      = (double **)pW[3];
    double *  theta_buf  = (double *) pW[4];
    double *  eps_hat    = (double *) pW[5];
    double *  w_hat      = (double *) pW[6];
    double *  v_hat      = (double *) pW[7];


    /* Pohrana vrijednosti i-tog ulaza procesa u grbavu matricu, cija i-ta kolona  ima d[i] elemenata,
       sa ciljem kasnjenja i-tog ulaza za vrijednost d[i]. */
    for (int i = 0; i < nu; i++){
         u[i] = delay(u_buf[i], u0[i], d[i]+1);
    }

    // Kasnjenje izlaza procesa za jedan odabirak.
    y       = delay(y_buf, *y0, 1);

    /* Punjenje regresora definisanog sa (2.60) preko pokazivaca‘psi’. w_hat, eps_hat, v_hat 
       koji su odredjeni sa (2.57). */
    psi_k1  = regressor_grls( psi, u, y, *w_hat, *eps_hat, *v_hat, *na, nb, nu, *nf, *nc, *nd, *model);


    /* U slucaju izbora odgovarajucih metoda sa konacnom historijom duzime win_lenght, 
        regresore pohranjujemo u matricu regressor_matrix,koja cini FIFO listu, cime je
        njen izlaz (posljednja vrsta matrice) za N odabiraka zakasnjeli  regresor. */

    if (*method==1){ 
        if (*history==2 && (*model==1 || *model==2 || *model==4)){

            double ** regressor_matrix  = (double **)pW[8];
            double *  y0_vector         = (double *) pW[9];

            psi_kN               = regressor_flow( regressor_matrix, psi_k1, *win_length);
            y0kN                 = y0_flow( y0_vector, *y0, *win_length); 
        } 
    }


    // Uzimanje pohranjenog vektora parametara theta_buf (iz prethodnog koraka simulacije). 
    for (int i = 0; i < DIM; i++){
        theta_hat_k.X[i][0] = theta_buf[i];
    }


    //  Slucaj izbora ERLS i KF  metoda.
    if (*method==1 || *method==3){

            /* Ova dva metoda estimacije parametara u obzir uzimaju matricu kovarijanse P_k. 
               Ovdje se uzima pohranjena vrijednost P_buf (iz prethodnog koraka simulacije). */
            for (int i = 0; i < DIM; i++){
                for (int j = 0; j < DIM; j++){
                    P_k.X[i][j] = P_buf[i][j];
                }
            } 

        /* U slucaju da imamo izbor metode sa konacnom historijom, faktor iscezavanja je jednak 1,
           u protivnom je faktor iscezavanja fade jednak vrijednosti koja se zadaje 
           kao parametar bloka (lambda). */
        *history==2 ? fade = 1      : fade = *lambda; 

        /* U slucaju da je za metod estimacije uzet KF, kovarijansa procesnog suma pn_cov je 
           jednaka vrijednosti koja se zadaje kao parametar bloka (pnc). */
        *method==3  ? pn_cov = *pnc : pn_cov = 0;      // process noise covariance


        /* Implementacija rekurzivnog algoritma najmanjih kvadrata sa ekponencijalnim 
           faktorom iščezavanja (relacije (2.34)). */
        struct Variable psi_k1_tr = transpose(psi_k1);
        residual_k1.X[0][0] = *y0 - inner(psi_k1, theta_hat_k, 0, DIM);
        Pk_psik1            = multiply(P_k, psi_k1);
        psik1tr_Pk_psik1    = multiply(psi_k1_tr, Pk_psik1);

        gamma_k             = c( Pk_psik1, (1/(fade + psik1tr_Pk_psik1.X[0][0])));


        /* Ako je blok  estimacije parametara omogucen preko ulaza Enable, tada cemo vektor
           parametara theta_hat_k1 racunati na osnovu (2.34a). */  
           if (*Enable > 0) {
              theta_hat_k1 = add(theta_hat_k, multiply(gamma_k, residual_k1));
           } else {
              theta_hat_k1 = theta_hat_k;
           }

        // Azuriranje matrice kovarijanse P i vektora parametara theta.   
        P_k1                = c( subtract( P_k, multiply( gamma_k, multiply( transpose(psi_k1), P_k) ) ), (1/ fade));
        P_k                 = add(P_k1, c(eye(DIM), pn_cov));
        theta_hat_k         = theta_hat_k1;

        /* Sastavljanje psi_k1 regresora i cuvanje njegove vrijednosti za sljedeci 
           korak simulacije. Regresor se sastavlja od elemenata w_hat, eps_hat, v_hat
           odredjenih relacijama (2.57) (u ovisnosti od izbora modela). */
        regressor_assembly(psi_k1, theta_hat_k1, *y0, w_hat, eps_hat, v_hat, *na, nb, nu, *nf, *nc, *nd, *model);


        // Racunanje greske predikcije i mjere pobudjenosti za ERLS i KF algoritme.
        *Error = residual_k1.X[0][0];
        *Excitation = 1 / (fade + psik1tr_Pk_psik1.X[0][0]);


        // U slucaju izbora metode (ERLS ili KF) sa konacnom historijom (i odgovarajuceg modela).
        if (*history==2 && (*model==1 || *model==2 || *model==4)){
            *method==3  ? pn_cov = *pnc : pn_cov = 0;

            //Implementacija ERLS algoritma sa konacnom historijom (relacije (2.37)).
            struct Variable Pk_psikN, psikNtr_Pk_psikN;
            struct Variable psi_kN_tr  = transpose(psi_kN);
            residual_k1.X[0][0]        = y0kN - multiply(psi_kN_tr, theta_hat_k).X[0][0];
            Pk_psikN                   = multiply(P_k, psi_kN);
            psikNtr_Pk_psikN           = multiply(psi_kN_tr, Pk_psikN);
            gamma_k                    = c( Pk_psikN, (1/(1 - psikNtr_Pk_psikN.X[0][0])));

            /* Ako je blok  estimacije parametara omogucen, tada cemo vektor parametara theta_hat_k1 
               racunati na osnovu (2.36c). */ 
            if (*Enable > 0) {
               theta_hat_k1 = add(theta_hat_k, multiply(gamma_k, residual_k1));
            } else {
               theta_hat_k1 = theta_hat_k;
            }

            // Azuriranje matrice kovarijanse P i vektora parametara theta.   
            P_k1                       = add( P_k, multiply( gamma_k, multiply( psi_kN_tr, P_k) ) );
            P_k                        = add(P_k1, c(eye(DIM), pn_cov));
            theta_hat_k                = theta_hat_k1;
        }

        // Pohrana matrice kovarijanse P_k u P_buf koji se cuva za naredni korak simulacije.
         for (int i = 0; i < DIM; i++){
              for (int j = 0;j < DIM; j++){
                   P_buf[i][j] = P_k.X[i][j];
              }
         }

    // U slucaju izbora BFGS metode:
    } else if (*method==2){

        struct Variable regressor_win, y0_win;
        double ** regressor_matrix  = (double **)pW[8];
        double *  y0_vector         = (double *) pW[9];

        /* U funkciji criterion je implementirana fukcija cilja (A.9) koja izrazava mjeru
           odstupanja predikcije modela od izmjerene vrijednosti izlaza sistema y0. 
           BFGS optimizacijski algoritam racuna vektor parametara theta_hat_k za koji 
           funkcija cilja ima minimum. U racunanju vrijednosti funkcije cilja uzimaju se 
           regresori i vrijednosti izlaza od trenutka k-N+1 do k-tog (tekuceg) trenutka.
           Ovi regresori su struktuirani u matricu regressor_win, a izlazi u y0_win.*/
        regressor_win   = regressor_window( regressor_matrix, psi_k1, *win_length);
        y0_win          = y0_window( y0_vector, y0[0], *win_length);

        /* Ako je blok estimacije parametara omogucen, tada cemo primjenom BFGS algoritma
           vektor parametara theta_hat_k1 racunati kao vrijednost za koju odstupanje
           predikcije modela od izmjerene vrijednosti ima minimum. */
           if (*Enable>0){
           theta_hat_k1    = BFGS( criterion, regressor_win, y0_win, theta_hat_k, *tolerance);
        } else{
           theta_hat_k1 = theta_hat_k;
        }  

         /* Sastavljanje psi_k1 regresora i cuvanje njegove vrijednosti za sljedeci 
            korak simulacije. Regresor se sastavlja od elemenata w_hat, eps_hat, v_hat
            odredjenih relacijama (2.57) (u ovisnosti od izbora modela). */
        regressor_assembly(psi_k1, theta_hat_k1, *y0, w_hat, eps_hat, v_hat, *na, nb, nu, *nf, *nc, *nd, *model);

        // Azuriranje vektora parametara theta.
        theta_hat_k     = theta_hat_k1;

        // Racunanje greske predikcije parametara za BGFS algoritam.
        struct Variable ss  = subtract(y0_win, multiply(regressor_win, theta_hat_k));
        double dd = 0;

        for (int i = 0; i < regressor_win.X.size(); i++){
             dd += ss.X[i][0];
        }
        *Error   = dd;
    } // <-- Izbor  metode: 1, 2, 3

     /* Ako je blok  estimacije parametara omogucen, tada cemo za sljedeci korak 
        simulacije sacuvati vektor estimiranih parametara iz tekuceg koraka. */
     if (*Enable > 0){
         for (int i = 0; i < DIM; i++){
              theta_buf[i] = theta_hat_k.X[i][0];
         } 
     }


     /* Izlaz bloka S-funkcije implementirane u ovom radu je sabirnica PSAU_BUS, 
        koju cine A, B, F, C, D estimirani parametri, u ovisnosti od izbora modela.
        Ovdje se vrsi struktuiranje estimiranih parametara u sabirnicu Params.  */

     if (*model==1){

         Params->A[0] = 1;
         for (int i=0; i<*na; i++){
            Params->A[i+1] = theta_buf[i];
         }
         B_params( Params, theta_hat_k, nb, nu, *na);

     } else if (*model==2){

         Params->A[0] = 1;
         for (int i=0; i<*na; i++){
            Params->A[i+1] = theta_buf[i];
         }   
         B_params( Params, theta_hat_k, nb, nu, *na);
         Params->C[0] = 1;
         for (int i=*na+snb; i<*na+snb+*nc; i++){
            Params->C[i+1-*na-snb] = theta_buf[i];
         }

     } else if (*model==3){  

         B_params( Params, theta_hat_k, nb, nu, 0);
         Params->F[0] = 1;
         for (int i=snb; i<snb+*nf; i++){
            Params->F[i+1-snb] = theta_buf[i];
         } 
         Params->C[0] = 1;
         for (int i=snb+*nf; i<snb+*nf+*nc; i++){
            Params->C[i+1-snb-*nf] = theta_buf[i];
         }
         Params->D[0] = 1;
         for (int i=snb+*nf+*nc; i<snb+*nf+*nc+*nd; i++){
            Params->D[i+1-snb-*nf-*nc] = theta_buf[i];
         }

     } else if (*model==4){

         B_params( Params, theta_hat_k, nb, nu, 0);
         Params->F[0] = 1;
         for (int i=snb; i<snb+*nf; i++){
            Params->B[i+1-snb] = theta_buf[i];
         }
     }

     Params->Ts = *sample_time;
}


/* Naredna funkcija je automatski generirana pritiskom tastera Build u bloku S-Function Builder. */
void Recursive_Estimation_Terminate_wrapper(void **pW,
			const int32_T *na, const int_T p_width0,
			const int32_T *nb, const int_T p_width1,
			const int32_T *nf, const int_T p_width2,
			const int32_T *nc, const int_T p_width3,
			const int32_T *nd, const int_T p_width4,
			const int32_T *d, const int_T p_width5,
			const real_T *lambda, const int_T p_width6,
			const real_T *alpha, const int_T p_width7,
			const int32_T *model, const int_T p_width8,
			const int32_T *method, const int_T p_width9,
			const int32_T *win_length, const int_T p_width10,
			const real_T *tolerance, const int_T p_width11,
			const real_T *sample_time, const int_T p_width12,
			const int32_T *history, const int_T p_width13,
			const real_T *pnc, const int_T p_width14,
			SimStruct *S)
{
free(pW[0]);
free(pW[1]);
free(pW[2]);
free(pW[3]);
free(pW[4]);
free(pW[5]);
free(pW[6]);
free(pW[7]);
free(pW[8]);
free(pW[9]);
}

