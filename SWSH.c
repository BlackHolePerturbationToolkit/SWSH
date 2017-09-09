/********************************************

13/07/2017
Written by  Conor O'Toole
			University College Dublin

Supported by the ESA Summer of Code in Space 2017

********************************************/

#include "SWSH.h"

/**************************************/
/*           Leaver's Method          */
/**************************************/

double alpha_n(int n, int m, int s, double a_omega){
	double alpha;
	double kplus, kminus;
	
	kplus = (double)abs(m + s);
	kminus = (double)abs(m - s);
	
	if(n < 0) return 0.; 
	
	alpha = -4.*a_omega*((double)n + kplus + 1.)*((double)n + kminus + 1.)*((double)n + (kplus + kminus)/2. + 1. + (double)s);
	alpha /= (2.*(double)n + kplus + kminus + 2.)*(2.*(double)n + kplus + kminus + 3.);
	
	
	return alpha;
}

double beta_n(double A, int n, int m, int s, double a_omega){
	double beta;
	double kplus, kminus;
	
	kplus = (double)abs(m + s);
	kminus = (double)abs(m - s);
	
	beta = A + (double)s*((double)s + 1.) + a_omega*a_omega - ((double)n + (kplus+kminus)/2.)*((double)n + (kplus + kminus)/2. + 1.);
	if(s != 0) beta += (8.*(double)m*(double)s*(double)s*a_omega)/((2.*(double)n + kplus + kminus)*(2.*(double)n + kplus + kminus + 2.));
	
	return beta;
}

double gamma_n(int n, int m, int s, double a_omega){
	double gamma;
	double kplus, kminus;
	
	kplus = (double)abs(m + s);
	kminus = (double)abs(m - s);
	
	gamma = 4.*a_omega*(double)n*((double)n + kplus + kminus)*((double)n + (kplus+kminus)/2. - (double)s);
	gamma /= (2.*(double)n + kplus + kminus -1.)*(2.*(double)n + kplus + kminus);
	
	return gamma;
}

double Cont_Frac_Nth_inversion(double a, int Nmax, int Ninv, int m, int s, double a_omega){

	int i;
	double CF1, CF2;
	
	CF1 = 0.;
	for(i = 0; i < Ninv; i++){
		CF1 = gamma_n(i+1, m, s, a_omega)/(beta_n(a, i, m, s, a_omega) - alpha_n(i-1, m, s, a_omega) * CF1);
	}
	CF1 *= alpha_n(Ninv-1, m, s, a_omega);
	
	CF2 = 1.;
	for( i = Nmax; i > Ninv; i--){
		CF2 = beta_n(a, i-1, m, s, a_omega) - (alpha_n(i-1, m, s, a_omega)*gamma_n(i, m, s, a_omega))/CF2;
	}
	
	return CF2 - CF1;		//Returns the difference of the n-th inversion and the rest of the continued fraction
							//To solve for the eigenvalue, set to zero.
}

double Newton_Raphson(int nmax, int ninv, int m, int s, double a_omega, double A0, double A1){
						
	double delta = 1.0e-6;		//The value used to numerically calculate the derivative of the continued fraction
	double prec = 1.0e-10;		//Will iterate until the continued fraction returns a value less than this.
	double CFprime;
	double A;
	
	A = A0;
	do{
		
		CFprime = Cont_Frac_Nth_inversion(A + 0.5*delta, nmax, ninv, m, s, a_omega) - Cont_Frac_Nth_inversion(A - 0.5*delta, nmax, ninv, m, s, a_omega);
		CFprime /= delta;
		
		A -= Cont_Frac_Nth_inversion(A, nmax, ninv, m, s, a_omega)/CFprime; 
		
	}while(fabs(Cont_Frac_Nth_inversion(A, nmax, ninv, m, s, a_omega)) > prec);
	
	return A;
}

double Secant_method(int nmax, int ninv, int m, int s, double a_omega, double A0, double A1){
							
	double a0, a1;
	double A;
	double prec = 1.0e-10;		//Will iterate until the continued fraction returns a value less than this.
	
	a0 = A0;			//Initial guesses either side of actual root
	a1 = A1;
	
	do{
		A = a0*Cont_Frac_Nth_inversion(a1, nmax, ninv, m, s, a_omega) - a1*Cont_Frac_Nth_inversion(a0, nmax, ninv, m, s, a_omega);
		A /= (Cont_Frac_Nth_inversion(a1, nmax, ninv, m, s, a_omega) - Cont_Frac_Nth_inversion(a0, nmax, ninv, m, s, a_omega));
		
		a0 = a1;
		a1 = A;
		
	}while(fabs(Cont_Frac_Nth_inversion(A, nmax, ninv, m, s, a_omega)) > prec);
	
	return A;
}


double SWSH_Eigenvalue_Leaver(int nmax, int ninv, int m, int s, double a_omega, double A0, double A1,
								double (*solver)(int Nmax, int Ninv, int M, int S, double A_Omega, double a0, double a1)){
	//A0 and A1 are the initial guesses. Both solvers take each as an input, to simplify the choice of solver. May improve later
	//Ideally, they will be either side of the actual root.
	
	double A;
	
	A = solver(nmax, ninv, m, s, a_omega, A0, A1);
	
	return A;
}

/**************************************/
/*           Spectral Method          */
/**************************************/

double clm2m_2(double l, double m, double s){
	double c;
	
	if( l < 2. || l + m < 2. || l - m < 2. || l-s < 2. || l + s < 2. ){
		return 0.;
	}
	else{
		c = sqrt((2.*l - 3.)/(2.*l + 1.));
		c *= sqrt((l - m)*(l + m)*(l - m - 1.)*(l + m - 1.));
		c *= sqrt((l - s)*(l + s)*(l - s - 1.)*(l + s - 1.));
		c /= (l*(2.*l - 1.)*(l - 1.)*(2.*l - 3.));
		
		return c;
	}
}

double clm1m_2(double l, double m, double s){
	double c;
	
	if( l < 3./2. || l + m < 1. || l - m < 1. || l - s < 1. || l + s < 1. ){
		return 0.;
	}
	else{
		c = -2. * m * s * sqrt((2.*l + 1.)/(2.*l - 1.));
		c *= sqrt((l - m)*(l + m));
		c *= sqrt((l - s)*(l + s));
		c /= (l*(2.*l*l*l + l*l -2.* l - 1.));
		
		return c;
	}
}

double clm_2(double l, double m, double s){
	double c;
	
	if( l < 1. || l + m < 0. || l - m < 0. || l-s < 0. || l + s < 0. ){
		return 1./3.;
	}
	else{
		c = 2.*(l*l + l -3.*m*m)*(l*l + l - 3.*s*s);
		c /= (l*(2.*l + 3.)*(l + 1.)*(2.*l - 1.));
		c += 1.;
		c /= 3.;
		
		return c;
	}
}

double clp1m_2(double l, double m, double s){
	double c;
	
	if( l < 1./2. || l + m < 0. || l - m < 0. || l-s < 0. || l + s < 0. ){
		return 0.;
	}
	else{
		c = -2. * m * s * sqrt((2.*l + 3.)/(2.*l + 1.));
		c *= sqrt((l - m + 1.)*(l + m + 1.));
		c *= sqrt((l - s + 1.)*(l + s + 1.));
		c /= (l*(l + 2.)*(2.*l + 3.)*(l + 1.));
		
		return c;
	}
}

double clp2m_2(double l, double m, double s){
	double c;
	
	if( l < 0. || l + m < 0. || l - m < 0. || l - s < 0. || l + s < 0. ){
		return 0.;
	}
	else{
		c = sqrt((2.*l + 5.)/(2.*l + 1.));
		c *= sqrt((l + m + 2.)*(l + m + 1.)*(l - m + 2.)*(l - m + 1.));
		c *= sqrt((l + s + 2.)*(l + s + 1.)*(l - s + 2.)*(l - s + 1.));
		c /= ((2.*l + 5.)*(l + 2.)*(2.*l + 3.)*(l + 1.));
		
		return c;
	}
}

double clm1m_1(double l, double m, double s){
	double c;
	
	if( l < 1. || l + m < 1. || l - m < 1. || l - s < 1. || l + s < 1. ){
		return 0.;
	}
	else{
		c = 1.*sqrt((2.*l - 1.)/(2.*l + 1.));
		c *= sqrt((l - m)*(l + m));
		c *= sqrt((l - s)*(l + s));
		c /= (l*(2.*l - 1.));
		
		return c;
	}
}

double clm_1(double l, double m, double s){
	double c;
	
	if( l < 0.5 || l + m < 0. || l - m < 0. || l-s < 0. || l + s < 0. ){
		return 0.;
	}
	else{
		c = -1. * m * s;
		c /= (l*(l + 1.));
		
		return c;
	}
}

double clp1m_1(double l, double m, double s){
	double c;
	
	if( l < 0. || l + m < 0. || l - m < 0. || l - s < 0. || l + s < 0. ){
		return 0.;
	}
	else{
		c = 1.*sqrt((2.*l + 3.)/(2.*l + 1.));
		c *= sqrt((l - m + 1.)*(l + m + 1.));
		c *= sqrt((l - s + 1.)*(l + s + 1.));
		c /= ((2.*l + 3.)*(l + 1.));
		
		return c;
	}
}

void matrix_mult(int dim, double mat1[dim][dim], double mat2[dim][dim], double mat_out[dim][dim]){
	int i, j, k;
	
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			mat_out[i][j] = 0.;
			
			for(k = 0; k < dim; k++){
				mat_out[i][j] += mat1[i][k]*mat2[k][j];
			}
		}
	}
}

void SpectralMatrix(int dim, double mat[dim][dim], int l, int m, int s, double a_omega, int nmin){
	int i, j;
	
	for( i = 0; i < dim; i ++ ){
		for( j = 0; j < dim; j++){
			if( i == j){
				mat[i][j] = a_omega*a_omega*clm_2((double)(i + l - nmin), (double)m, (double)s) - 2.*a_omega*s*clm_1((double)(i + l - nmin), (double)m, (double)s) - (double)(i + l - nmin)*((double)(i + l - nmin) + 1.);
			}
			else if( j - i == 1){	//ie. one column further right, so one above the diagonal
				mat[i][j] = a_omega*a_omega*clp1m_2((double)(i + l - nmin), (double)m, (double)s) - 2.*a_omega*s*clp1m_1((double)(i + l - nmin), (double)m, (double)s);
			}
			else if( j - i == -1){	//ie. one column further left, so one below the diagonal
				mat[i][j] = a_omega*a_omega*clm1m_2((double)(i + l - nmin), (double)m, (double)s) - 2.*a_omega*s*clm1m_1((double)(i + l - nmin), (double)m, (double)s);
			}
			else if( j - i == 2){	//ie. two columns further right, so two above the diagonal
				mat[i][j] = a_omega*a_omega*clp2m_2((double)(i + l - nmin), (double)m, (double)s);
			}
			else if( j - i == -2){	//ie. two columns further left, so one below the diagonal
				mat[i][j] = a_omega*a_omega*clm2m_2((double)(i + l - nmin), (double)m, (double)s);
			}
			else{
				mat[i][j] = 0.;
			}
		}
	}
}

void Givens_Rotation(int dim, double mat[dim][dim], int i, int j, double c, double s){
	int m, n;
	
	for(m = 0; m < dim; m++){
		for(n = 0; n < dim; n++){
			if(m == i && n == i) mat[m][n] = c;
			else if(m == j && n == j) mat[m][n] = c;
			//else if(m == i && n == j) mat[m][n] = -1.*s;	
			//else if(m == j && n == i) mat[m][n] = s;
			else if(m == i && n == j) mat[m][n] = s;	//These set the Given's rotation identical to the definition in Numerical Recipes Ch11.
			else if(m == j && n == i) mat[m][n] = -1.*s;
			else if(m == n && m != i && m != j) mat[m][n] = 1.;
			else mat[m][n] = 0.;
		}
	}
}

void Reduction5_3(int dim, double input[dim][dim], double output[dim][dim]){
	int i, j, m, n;
	double theta, t, sign;
	double r, c, s;
	double mid[dim][dim], mid2[dim][dim], rot[dim][dim], rot_T[dim][dim];
	double temp1, temp2, temp3, temp4;
	
	for(m = 0; m < dim; m++){
		for(n = 0; n < dim; n++){
			mid[m][n] = input[m][n];
		}
	}
	
	//This algorithm is defined to be identical to that in Numerical Recipes Ch 11.
	//There are issues with this algorithm. It is not advised to use this function for now.
	for(j = 0; j < dim-2; j++){	
		for(i = j+2; i < dim; i++){
			if(mid[i][j] != 0.){
				
				theta = (mid[j][j] - mid[i][i])/(2.*mid[i][j]);
				if(theta < 0.) sign = -1.;
				else sign = 1.;
				
				t = sign/(fabs(theta) + sqrt(theta*theta + 1.));
				c = 1./sqrt(t*t +1.);
				s = t*c;					//Non-trivial elements of spectral matrix
				
				Givens_Rotation(dim, rot, j+1, i, c, s);
				
				for(m = 0; m < dim; m++){
					for(n = 0; n < dim; n++){
						rot_T[m][n] = rot[n][m];
					}
				}
				
				matrix_mult(dim, rot, mid, mid2);	//Eliminate sub-diagonal element
				matrix_mult(dim, mid2, rot_T, mid);	//Eliminate super-diagonal element
				
				for(m = 0; m < dim; m++){
					for(n = 0; n < dim; n++){
						if(fabs(mid[m][n]) < 1.0e-10 ) mid[m][n] = 0.; 	//Simple check to allow for rounding error
					}
				}
			}
		}
	}
	
	double temp_col[dim];
	
	for(m = 0; m < dim; m++){
		for(n = 0; n < dim; n++){
			output[m][n] = mid[m][n];
		}
	}
	
}

void QR_algorithm(int dim, double input[dim][dim], double eigenvalues[dim], int spin, double orthog[dim][dim]){
	int i, j, m, n;
	double r, c, s;
	double Q[dim][dim], Q_T[dim][dim], G[dim][dim], G_T[dim][dim];
	double mid[dim][dim];
	double input_k[dim][dim];	// the result of the k-th step of the QR algorithm input_k+1 = Q_T*input_k*Q
	double convergence_check;
	double precision = 1.0e-10;
	
	for(m = 0; m < dim; m++){
		for(n = 0; n < dim; n++){
			input_k[m][n] = input[m][n];	//This will be input_0
		}
	}
	
	for(m = 0; m < dim; m++){
			for(n = 0; n < dim; n++){
				if(m == n) Q_T[m][n] = 1.;
				else Q_T[m][n] = 0.;	//initialize Q to identity, so that we can use *= in the following algorithm
			}
		}
	do{
		convergence_check = 0.;
		for( i = dim - 1; i > 0; i--){
			convergence_check += fabs(input_k[i][i-1]);	//all elements below the diagonal must be below precision
		}
		
		//Issues with this yet, must be corrected before it can be used.
		//Will likely be corrected after the reduction algorithm
		for( j = 0; j < dim - 1; j++){
			for( i = dim-1; i > j ; i--){
				if(input_k[i][j] != 0.){
				
					r = sqrt(input_k[i-1][j]*input_k[i-1][j] + input_k[i][j]*input_k[i][j]);
	
					c = input_k[i-1][j]/r;
					s = -1.*input_k[i][j]/r;
	
					Givens_Rotation(dim, G, i-1, i, c, s);
					
					for(m = 0; m < dim; m++){
						for(n = 0; n < dim; n++){
							G_T[m][n] = G[n][m];
						}
					}
				
					matrix_mult(dim, G, input_k, mid);
					matrix_mult(dim, mid, G_T, input_k);
					
					for(m = 0; m < dim; m++){
						for(n = 0; n < dim; n++){
							if(fabs(input_k[m][n]) < 1.0e-10 ) input_k[m][n] = 0.; 
						}
					}
				
					matrix_mult(dim, G, Q_T, mid);	//at each step will carry out Q*=G, where G is the rotation which eliminates input_k[i][i-1]
			
					for(m = 0; m < dim; m++){
						for(n = 0; n < dim; n++){
							Q_T[m][n] = mid[m][n];
						}
					}
				}
			}
		}
		
		for(m = 0; m < dim; m++){
			for(n = 0; n < dim; n++){
				Q[m][n] = Q_T[n][m];	//calculate transpose of Q
			}
		}
	} while(convergence_check > (double)(dim-1)*precision);
	
	//Printing functions for checking issues with this code.
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			printf("%.4lf ", input_k[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			printf("%.4lf ", Q[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
	for(i = 0; i<dim; i++){
		for(j = 0; j < dim; j++){
			orthog[i][j] = Q[i][j];		//Matrix with eigenvectors in columns
		}
	}
	
	for(i = 0; i < dim; i++){
		eigenvalues[i] = input_k[i][i];	//Returns eigenvalues.
	}
	
	
}

double SWSH_Eigenvalue_Eigenvector_Spectral_custom(int l, int m, int s, double a_omega, double eigenvector[]){
	//This function is not considered accurate currently. Issues with reduction and QR algorithm still need to be corrected
	int lmin, nmin, nmax;
	int i, j, k;
	double temp;
	
	if(abs(s) >= abs(m)) lmin = abs(s);
	else lmin = abs(m);
	nmax = ceil(fabs(1.5 * a_omega - (a_omega*a_omega)/250.)) + 5;
	if(nmax <= l-lmin) nmin = nmax;
	else nmin = l - lmin;
	
	int dim = nmin + nmax + 2;		//Calculation of matrix dimension.
									//May write in seperate function at a later point
	
	double mat[dim][dim], reduced_mat[dim][dim];
	double eigenvalues[dim];
	double orthog[dim][dim];
	double temp_vec[dim];				//May change all declarations of matrices to dynamic memory allocation to improve efficiency later.
										//Would have a significant effect for larger values of l and gamma, where the spectral matrix could be quite large.
	
	SpectralMatrix(dim, mat, l, m, s, a_omega, nmin);	//Generate spectral matrix
	
	//Printing for checking issues.
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++ ){
			printf("%.4lf ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
	Reduction5_3(dim, mat, reduced_mat);	//Reduction to tridiagonal
	//Printing for checking issues.
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++ ){
			printf("%.4lf ", reduced_mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
	QR_algorithm(dim, reduced_mat, eigenvalues, s, orthog);	//Apply QR algorithm
	
	for( i = 0; i < dim; i++ ){
		eigenvalues[i] *= -1.;
	}
	
	for ( i = 0; i < dim; i++ ){
		for(j = 0; j < dim; j++ ){
			if (eigenvalues[i] > eigenvalues[j]){
				
				temp = eigenvalues[i];
				
				eigenvalues[i] = eigenvalues[j];
				eigenvalues[j] = temp;
				
				for(k = 0; k < dim; k++){
					temp_vec[k] = orthog[k][i];
					
					orthog[k][i] = orthog[k][j];
					orthog[k][j] = temp_vec[k];
				}				//Sort eigenvalues and eigenvectors
			}					//The eigenvalue corresponding to l, m will be the l-th biggest
		}
	}
	
	for(k = 0; k < dim; k++){
		eigenvector[k] = orthog[k][dim-nmin-1];	
	}
	
	return eigenvalues[dim-nmin-1] - (double)(s*(s+1));		//This returns the eigenvalue Alm, which is the definition used here.
}

double SWSH_Eigenvalue_Eigenvector_Spectral_gsl(int l, int m, int s, double a_omega, double eigenvector[]){
	int lmin, nmin, nmax;
	double eig;
	int i, j;
	
	if(abs(s) >= abs(m)) lmin = abs(s);
	else lmin = abs(m);
	nmax = ceil(fabs(1.5 * a_omega - (a_omega*a_omega)/250.)) + 5;
	if(nmax <= l-lmin) nmin = nmax;
	else nmin = l - lmin;
	
	int dim = nmin + nmax + 2;			//dimension calculation
	
	double matrix[dim][dim], reduced_mat[dim][dim];
	
	SpectralMatrix(dim, matrix, l, m, s, a_omega, nmin);
	
	gsl_matrix_view mat = gsl_matrix_view_array (*matrix, dim, dim);	//Type casting to necessary GSL types to use standard GSl functions
	
	gsl_vector *eval = gsl_vector_alloc (dim);
	gsl_matrix *evec = gsl_matrix_alloc (dim, dim);
	
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dim);
	
	gsl_eigen_symmv(&mat.matrix, eval, evec, w);		//Eigensystem

	gsl_eigen_symmv_free(w);

	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC); 	//As we will be multiplying by -1., use ASC, not DESC
																//Eigenvectors will also be sorted

	eig = gsl_vector_get(eval, dim-nmin-1);
	gsl_vector_view eigvec = gsl_matrix_column (evec, dim-nmin-1);	//Obtain eigenvalue/eigenvector corresponding to given l, m
	
	for( i = 0; i < dim; i++){
		eigenvector[i] = gsl_vector_get(&eigvec.vector, i);		//Assign to eigenvector input argument.
	}
	
	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	
	return -1.*eig - (double)(s*(s+1));
	
}

/**************************************/
/* Spin-Weighted Spherical Harmonics  */
/**************************************/

double A(int l, int m){
	
	double mid1, mid2;
	int i;
	
	mid1 = 1.;
	mid2 = 1.;
	
	for(i = l-m; i > 0; i--){
		mid1 *= (double)i;
	}
	
	for(i = l+m; i > 0; i--){
		mid2 *= (double)i;
	}
	
	return sqrt(((2.* l + 1) * mid1)/( 4. * M_PI * mid2));
}

double P(int l, int m, double x){
	
	int i, j, k;
	double mid;
	double Legendre[l][l];		//dynamic memory allocation on each of the first indices would be far more efficient, but simplicity is desired for now
	
	for(i = 0; i <= l; i++){
		for(j = 0; j <= GSL_MIN_INT(i, m); j++){
			
			mid = 1.;
			
			for(k = 2*j-1; k > 0; k -= 2){
				mid *= (double)k;
			}
			
			if(i == j) Legendre[i][j] = pow(-1., (double)j) * mid * pow(1.- x*x, (double)j/2.);
			else if(i == j+1 ) Legendre[i][j] = x*(2*(double)j + 1.)*Legendre[i-1][j];
			else Legendre[i][j] = (1./(double)(i-j))*( x*(2.*(double)i -1.)*Legendre[i-1][j] - (double)(i+j-1) * Legendre[i-2][j] ); 
		}
	}
	
	return Legendre[l][m];
}

double d(int l, int m, double x){
	
	int i, j, k;
	double mid;
	double d_Legendre[l][l];
	
	for(i = 0; i <= l; i++){
		for(j = 0; j <= GSL_MIN_INT(i, m); j++){
			
			mid = 1.;
			
			for(k = 2*j-1; k > 0; k -= 2){
				mid *= (double)k;
			}
			
			if(i == j) d_Legendre[i][j] = -j * x * pow(-1., (double)j) * mid * pow(1.- x*x, (double)(j-2)/2.);
			else if(i == j+1 ) d_Legendre[i][j] = (2*(double)j + 1.)*(P(i-1, j, x) + x*d_Legendre[i-1][j]);
			else d_Legendre[i][j] = (1./(double)(i-j))*( (2.*(double)i -1.)*(P(i-1, j, x) + x*d_Legendre[i-1][j]) - (double)(i+j-1) * d_Legendre[i-2][j] ); 
		}
	}
	
	return d_Legendre[l][m];
}

double D(int l, int m, double x){
	
	int i, j, k;
	double mid;
	double D_Legendre[l][l];
	
	for(i = 0; i <= l; i++){
		for(j = 0; j <= GSL_MIN_INT(i, m); j++){
			
			mid = 1.;
			
			for(k = 2*j-1; k > 0; k -= 2){
				mid *= (double)k;
			}
			
			if(i == j) D_Legendre[i][j] = j * pow(-1., (double)j) * mid * pow(1.- x*x, (double)(j-4)/2.) * ( ((double)j - 2.) * x*x - (1. - x*x) );
			else if(i == j+1 ) D_Legendre[i][j] = (2*(double)j + 1.)*( 2.*d(i-1, j, x) + x*D_Legendre[i-1][j]);
			else D_Legendre[i][j] = (1./(double)(i-j))*( (2.*(double)i -1.)*( 2.* d(i-1, j, x) + x*D_Legendre[i-1][j]) - (double)(i+j-1) * D_Legendre[i-2][j] ); 
		}
	}
	
	return D_Legendre[l][m];
}

double s0_Ylm(int l, int m, double theta){
	
	double c_theta;
	
	c_theta = cos(theta);
	
	if(c_theta > 1.) c_theta = 1.;
	if(c_theta < -1.) c_theta = -1.;
	
	return A(l, m)*P(l, m, c_theta);
	
}

double s1_Ylm(int l, int m, int s, double theta){
	
	double sigma;
	double c_theta, s_theta;
	double mid1, mid2;
	
	if(s > 0) sigma = 1.;
	else sigma = -1.;
	
	c_theta = cos(theta);
	s_theta = sin(theta);
	
	if(c_theta > 1.) c_theta = 1.;
	if(c_theta < -1.) c_theta = -1.;
	
	if(s_theta > 1.) s_theta = 1.;
	if(s_theta < -1.) s_theta = -1.;
	
	if(fabs(s_theta) < 1.0e-6) return 0.;
	else{
		mid1 = s_theta * d(l, m, c_theta);
		mid2 = ((sigma*(double)m)/s_theta)*P(l, m, c_theta);
		
		return ((sigma*A(l, m))/(sqrt((double)(l*(l+1))))) * (mid1 + mid2);
	}
}

double s2_Ylm(int l, int m, int s, double theta){
	
	double sigma;
	double c_theta, s_theta;
	double mid1, mid2, mid3;
	
	if(s > 0) sigma = 1.;
	else sigma = -1.;
	
	c_theta = cos(theta);
	s_theta = sin(theta);
	
	if(c_theta > 1.) c_theta = 1.;
	if(c_theta < -1.) c_theta = -1.;
	
	if(s_theta > 1.) s_theta = 1.;
	if(s_theta < -1.) s_theta = -1.;
	
	if(fabs(s_theta) < 1.0e-6) return 0.;
	else{
		mid1 = s_theta*s_theta * D(l, m, c_theta);
		mid2 = sigma * 2. * m * d(l, m, c_theta);
		mid3 = ( ( m*m + sigma * 2. * m * c_theta )/( s_theta*s_theta ) ) * P(l, m, c_theta);
	
		return ( A(l, m)/( sqrt((double)((l-1)*l*(l+1)*(l+2)))) ) * (mid1 + mid2 + mid3);
	}
}

double SWSpherical_Harmonic(int l, int m, int s, double theta){
	
	double s_Y_lm;
	
	if(s == 0) s_Y_lm = s0_Ylm(l, m, theta);
	else if(s == -1 || s == 1) s_Y_lm = s1_Ylm(l, m, s, theta);
	else if(s == -2 || s == 2) s_Y_lm = s2_Ylm(l, m, s, theta);
	else {
		printf("WARNING@: The spin-weighted spheroidal harmonics for s > 2 have not been implemented in this code. \n");
		s_Y_lm = 0.;
	}
	
	return s_Y_lm;
}
