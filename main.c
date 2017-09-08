/********************************************

21/08/2017
Written by  Conor O'Toole
			University College Dublin

Supported by the ESA Summer of Code in Space 2017

********************************************/

#include "SWSH.h"

int main(int argc, char *argv[]){
	
	clock_t begin = clock();
	
	int l, m, s;	
	int lmin, nmin, nmax;
	int i, j;
	int Ninv;

  double a_omega;
	
	l = 2;
	m = 0;
	s = -2;
	a_omega = 5.5;
	
	Ninv = l - fabs(m);
	
	if(fabs(s) >= fabs(m)) lmin = fabs(s);
	else lmin = fabs(m);
	nmax = ceil(fabs(1.5 * a_omega - (a_omega*a_omega)/250.)) + 5;
	if(nmax <= l-lmin) nmin = nmax;
	else nmin = l - lmin;
	
	int dim = nmin + nmax + 2;
	
	double E_1, E_2, E_3, E_4;
	double E_vec[dim], E_vec_2[dim];
	
	E_1 = SWSH_Eigenvalue_Eigenvector_Spectral_gsl(l, m, s, a_omega, E_vec);	//Most accurate method currently
	
	
	//E_2 = SWSH_Eigenvalue_Leaver(200, Ninv, m, s, a_omega, E_3 - 1./pow(a_omega, 10.), E_3 + 1./pow(a_omega, 10.), Secant_method);// - 2.*m*a_omega + a_omega*a_omega;
	//E_3 = SWSH_Eigenvalue_Leaver(200, Ninv, m, s, a_omega, E_3 - 1./pow(a_omega, 10.), E_3 + 1./pow(a_omega, 10.), Newton_Raphson);
	
	//printf("Leaver: %.20lf \n", E_2);
	//printf("Leaver: %.20lf \n", E_3);
	
	E_4 = SWSH_Eigenvalue_Eigenvector_Spectral_custom(l, m, s, a_omega, E_vec_2);
	
	
	//Printing to allow for checking calculated values
	printf("Spectral (gsl): %.20lf \n\n", E_1);
	printf("Spectral (custom): %.20lf \n", E_4);
	printf("Eigenvalue Difference: %.20lf \n \n", fabs(E_1-E_4));
	printf("GSL                        Custom                        Difference \n");
	for (i = 0; i < dim; i++){
		printf("%.20lf %.20lf %.20lf \n", E_vec[i], E_vec_2[i], fabs(fabs(E_vec[i]) - fabs(E_vec_2[i])));
	}
	printf("\n");
	
	clock_t end = clock();
	double duration = (double)(end - begin)/CLOCKS_PER_SEC;
	printf("Duration = %lf s \n", duration);
	fflush(stdout);
	
	return 0;
}