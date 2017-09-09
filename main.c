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
	const char* spherical_file = "Spherical_Harmonic.txt";
	
	Spherical_test = fopen(spherical_file, "w");
	
	l = 20;
	m = 4;
	s = 0;
	a_omega = 5.5;
	
	double Re_Y, Im_Y, c_phi, s_phi;
	double Re_Y_2;
	
	//double Y2, diff, p;
	
	//diff = d(l, m, cos(0.5));
	//Y2 = s2_Ylm(l, m, s, 0.5);
	
	//printf("%lf, %lf \n", diff, Y2);
	
	theta = 0.;
	do{
		phi = 0.;
		c_phi = cos(phi);
		s_phi = sin(phi);
		
		if(c_phi > 1.) c_phi = 1.;
		if(c_phi < -1.) c_phi = -1.;
		if(s_phi > 1.) s_phi = 1.;
		if(s_phi < -1.) s_phi = -1.;
		
		do{
			//Re_Y = c_phi * SWSpherical_Harmonic(l, m, s, theta);
			//Im_Y = s_phi * SWSpherical_Harmonic(l, m, s, theta);
			Re_Y_2 = c_phi * gsl_sf_legendre_sphPlm(l, m, cos(theta));
			
			//printf("%lf, %lf, %.20lf \n", Re_Y, Re_Y_2, fabs(Re_Y - Re_Y_2));
			
			//fprintf(Spherical_test, "%.15lf %.15lf %.15lf %.15lf \n", theta, phi, Re_Y, Im_Y);
			
			phi += M_PI/100.;
		} while(phi < M_PI);
		
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
