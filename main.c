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
	FILE *Spherical_test;
	const char* spherical_file = "Spherical_Harmonic_R.txt";
	
	Spherical_test = fopen(spherical_file, "w");
	
	l = 5;
	m = 0;
	s = 0;
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
			Re_Y_2 = c_phi * gsl_sf_legendre_sphPlm(l, m, c_theta);
			Im_Y_2 = s_phi * gsl_sf_legendre_sphPlm(l, m, c_theta);
			
			if(fabs(Re_Y - Re_Y_2) >= 1.0e-10) printf("%.20lf \n", fabs(Re_Y - Re_Y_2));
			if(fabs(Im_Y - Im_Y_2) >= 1.0e-10) printf("%.20lf \n", fabs(Im_Y - Im_Y_2));
			//printf("%lf, %lf, %.20lf \n", Re_Y, Re_Y_2, fabs(Re_Y - Re_Y_2));
			//printf("%lf, %lf, %.20lf \n", Im_Y, Im_Y_2, fabs(Im_Y - Im_Y_2));
			
			fprintf(Spherical_test_r, "%.15lf \n", Re_Y);
			fprintf(Spherical_test_i, "%.15lf \n", Im_Y);
			
			phi += M_PI/100.;
		} while(phi <= M_PI);
		
		theta += (2.*M_PI)/100.;
	} while(theta < 2.*M_PI);
	
	
	fclose(Spherical_test);
	free(Spherical_test);
	
	clock_t end = clock();
	double duration = (double)(end - begin)/CLOCKS_PER_SEC;
	printf("Duration = %lf s \n", duration);
	fflush(stdout);
	
	return 0;
}
