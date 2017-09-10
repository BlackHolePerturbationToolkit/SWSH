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
	double theta, phi;
	FILE *Spherical_test, *Spherical_test_theta, *Spherical_test_phi, *Spherical_test_r, *Spherical_test_i;
	const char* spherical_file = "Spherical_Harmonic.txt";
	const char* spherical_file_theta = "Spherical_Harmonic_Theta.txt";
	const char* spherical_file_phi = "Spherical_Harmonic_Phi.txt";
	const char* spherical_file_r = "Spherical_Harmonic_R.txt";
	const char* spherical_file_i = "Spherical_Harmonic_I.txt";
	
	Spherical_test = fopen(spherical_file, "w");
	Spherical_test_theta = fopen(spherical_file_theta, "w");
	Spherical_test_phi = fopen(spherical_file_phi, "w");
	Spherical_test_r = fopen(spherical_file_r, "w");
	Spherical_test_i = fopen(spherical_file_i, "w");
	
	l = 5;
	m = 0;
	s = -1;
	a_omega = 5.5;
	
	double Re_Y, Im_Y, c_phi, s_phi;
	double Re_Y_2, Im_Y_2;
	
	theta = 0.;
	do{
		phi = 0.;
		
		c_theta = cos(theta);
		
		if(c_theta > 1.) c_theta = 1.;
			if(c_theta < -1.) c_theta = -1.;
		
		do{
			c_phi = cos(((double)m)*phi);
			s_phi = sin(((double)m)*phi);
		
			if(c_phi > 1.) c_phi = 1.;
			if(c_phi < -1.) c_phi = -1.;
			if(s_phi > 1.) s_phi = 1.;
			if(s_phi < -1.) s_phi = -1.;
			
			Re_Y = c_phi * SWSpherical_Harmonic_Wigner(l, m, s, theta);
			Im_Y = s_phi * SWSpherical_Harmonic_Wigner(l, m, s, theta);

			if(s == 0){
				Re_Y_2 = c_phi * gsl_sf_legendre_sphPlm(l, m, cos(theta)); 
				Im_Y_2 = s_phi * gsl_sf_legendre_sphPlm(l, m, cos(theta));
			
				if(fabs(Re_Y - Re_Y_2) != 0.) printf("%.20lf \n", fabs(Re_Y - Re_Y_2));
			}
			
			fprintf(Spherical_test, "%.15lf %.15lf %.15lf %.15lf \n", theta, phi, Re_Y, Im_Y);
			fprintf(Spherical_test_theta, "%.15lf \n", theta);
			fprintf(Spherical_test_phi, "%.15lf \n", phi);


			fprintf(Spherical_test_r, "%.15lf \n", Re_Y);
			fprintf(Spherical_test_i, "%.15lf \n", Im_Y);
			
			phi += M_PI/100.;
		} while(phi <= M_PI);
		
		theta += (2.*M_PI)/100.;
	} while(theta < 2.*M_PI);
	
	fclose(Spherical_test);
	fclose(Spherical_test_theta);
	fclose(Spherical_test_phi);
	fclose(Spherical_test_r);
	fclose(Spherical_test_i);
	free(Spherical_test);
	free(Spherical_test_theta);
	free(Spherical_test_phi);
	free(Spherical_test_r);
	free(Spherical_test_i);
	
	clock_t end = clock();
	double duration = (double)(end - begin)/CLOCKS_PER_SEC;
	printf("Duration = %lf s \n", duration);
	fflush(stdout);
	
	return 0;
}
