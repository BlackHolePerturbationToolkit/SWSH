/********************************************

13/07/2017
Written by  Conor O'Toole
			University College Dublin
Supervisor:	Barry Wardell
			University College Dublin
			
Supported by the ESA Summer of Code in Space 2017

********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_legendre.h>

/**************************************/
/*           Leaver's Method          */
/**************************************/
/*

alpha_n
beta_n
gamma_n						Three polynomials necessary for Leaver's recursion relation

Cont_Frac_Nth_inversion		For a given l, m, the corresponding eigenvalue will be the most stable
							root of the Nth inversion, where N = l -|m|
							
Newton_Raphson
Secant_method				Root finding functions for the continued fractions  (do not support solving general problems)

SWSH_Eigenvalue_Leaver		Collation of all above functions to find eigenvalue corresponding to given l, m, s, gamma (denoted a_omega)
							Allows fo specification of root finding method, with Newton_Raphson and Secant_method being the only options currently.

*/
double alpha_n(int n, int m, int s, double a_omega);
double beta_n(double A, int n, int m, int s, double a_omega);
double gamma_n(int n, int m, int s, double a_omega);

double Cont_Frac_Nth_inversion(double a, int Nmax, int Ninv, int m, int s, double a_omega);

double Newton_Raphson(int nmax, int ninv, int m, int s, double a_omega, double A0, double A1);	//A0, A1 passed but only A1 used, to make it ewasier to call a specific solver in SWSH_Eigenvalue_Leaver
double Secant_method(int nmax, int ninv, int m, int s, double a_omega, double A0, double A1);

/**************************************/
/*           Spectral Method          */
/**************************************/
/*

clm2m_2
clm1m_2
clm_2
clp1m_2
clp2m_2
clm1m_1
clm_1
clp1m_1					Calculates coefficients of spectral matrix. Can be simplified to three functions (to be done later)

matrix_mult				Simple matrix multiplication functions

SpectralMatrix			Calculate the spectral matrix for a given l, m, s, gamma

Givens_Rotation			Calculate rotation matrix which will eliminate an off-diagonal element of the spectral matrix

Reduction5_3			Reduction of spectral matrix from pentadiagonal to ridiagonal

QR_algorithm			QR algorithm to calculate eigenvalues and eigenvectors of spectral matrix.
						Currently, issues with reduction algorithm, so QR algorithm is applied directly to the spectral matrix.
						Efficiency can be improved by applying to reduced matrix.
						
SWSH_Eigenvalue_Eigenvector_Spectral_custom		Function to collate the above and find all eigenvalues and eigenvectors of spectral matrix

SWSH_Eigenvalue_Eigenvector_Spectral_gsl		GSL version of the above, which uses standard GSL functions to find all eigenvalues and eigenvectors of spectral matrix
												Currently the preferred method.

*/
double clm2m_2(double l, double m, double s);
double clm1m_2(double l, double m, double s);
double clm_2(double l, double m, double s);
double clp1m_2(double l, double m, double s);
double clp2m_2(double l, double m, double s);
double clm1m_1(double l, double m, double s);
double clm_1(double l, double m, double s);
double clp1m_1(double l, double m, double s);

void matrix_mult(int dim, double mat1[dim][dim], double mat2[dim][dim], double mat_out[dim][dim]);
void SpectralMatrix(int dim, double mat[dim][dim], int l, int m, int s, double a_omega, int nmin);
void Givens_Rotation(int dim, double mat[dim][dim], int i, int j, double c, double s);
void Reduction5_3(int dim, double input[dim][dim], double output[dim][dim]);
void QR_algorithm(int dim, double input[dim][dim], double eigenvalues[dim], int spin, double orthog[dim][dim]);

double SWSH_Eigenvalue_Leaver(int nmax, int ninv, int m, int s, double a_omega, double A0, double A1,
								double (*solver)(int Nmax, int Ninv, int M, int S, double A_Omega, double a0, double a1));
								
double SWSH_Eigenvalue_Eigenvector_Spectral_custom(int l, int m, int s, double a_omega, double eigenvector[]);
double SWSH_Eigenvalue_Eigenvector_Spectral_gsl(int l, int m, int s, double a_omega, double eigenvector[]);

/******************************************************/
/*  Spin-Weighted Spherical Harmonics Legendre method */
/******************************************************/

double A(int l, int m);
double P(int l, int m, double x);
double d(int l, int m, double x);
double D(int l, int m, double x);
double s0_Ylm(int l, int m, double theta);
double s1_Ylm(int l, int m, int s, double theta);
double s2_Ylm(int l, int m, int s, double theta);
double SWSpherical_Harmonic_Hughes(int l, int m, int s, double theta);

/****************************************************/
/* Spin-Weighted Spherical Harmonics Wigner method  */
/****************************************************/

double SWSpherical_Harmonic_Wigner(int l, int m, int s, double theta);

/****************************************************/
/*        Spin-Weighted Spheroidal Harmonics        */
/****************************************************/

double SWSpheroidal_Harmonic(int l, int m, int s, double a_omega, double theta);
