/********************************************

21/08/2017
Written by  Conor O'Toole
			University College Dublin
Supervisor:	Barry Wardell
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
	double theta;
	int res;
	FILE *spheroidal_file;
	FILE *eigenvalues_file;
	char output_filename[50];
	char *eigenvalues_filename = "C_Test_Results/Spheroidal_Eigenvalues.txt";
	res = 100;

	//l = 6;
	//m = 4;
	//s = -1;
	a_omega = 0.5;
	//theta = M_PI/4.;
	
	eigenvalues_file = fopen(eigenvalues_filename, "w");
	
	double S;
	double Lambda;
	int dim;
	
	for(s = -2; s < 3; s++){
		for(l = abs(s); l < 21; l++){
			for(m = -l; m < l+1; m++){
				
				if(abs(s) >= abs(m)) lmin = abs(s);
				else lmin = abs(m);
				nmax = ceil(fabs(1.5 * a_omega - (a_omega*a_omega)/250.)) + 5;
				if(nmax <= l-lmin) nmin = nmax;
				else nmin = l - lmin;
	
				dim = nmin + nmax + 2;
				double eigenvector[dim];
			
				Lambda = SWSH_Eigenvalue_Eigenvector_Spectral_gsl(l, m, s, a_omega, eigenvector);
					
				fprintf(eigenvalues_file, "%.20lf \n", Lambda);
			}
		}
	}
	
	fclose(eigenvalues_file);
	free(eigenvalues_file);
	eigenvalues_file = NULL;
	
	clock_t end = clock();
	double duration = (double)(end - begin)/CLOCKS_PER_SEC;
	printf("Duration = %lf s \n", duration);
	fflush(stdout);
	
	return 0;
}
