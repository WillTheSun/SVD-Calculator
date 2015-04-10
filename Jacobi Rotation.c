/* William Sun
Saturday, October 18, 2014
Applying Jacobi method to find SVD of symetric square matrix
Using the standard Jacobi rotation method explained by
http://www.netlib.org/lapack/lawnspdf/lawn15.pdf */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* ind(i,j) takes 2d matrix indices and returns the index in a 1d representation */
int ind(int i, int j, int n){
	return i*n + j;
}

/* transpose returns the transpose of the matrix Q as the matrix transpose */
void transpose(double Q[], double transpose[], int n){
	int i, j;
	for (i = 0; i < n; i ++){
		for (j = 0; j < n; j++){
			transpose[ind(j,i,n)] = Q[ind(i,j,n)];
		}
	}
}

/* Create the identity matrix */
void identity(double *x, int n){
	int i, j;
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			if (i == j)
				x[ind(i,i,n)]= 1;
			else
				x[ind(i,j,n)]= 0;
		}
	}
}

/*Multiplies two square matrices, x and y, together and returns z, 
the result*/
void multiply(double x[], double y[], double z[], int n){
	double sum = 0;
	int i, j, k;
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			for (k=0; k<n; k++){
				sum = sum + x[ind(i,k,n)] * y[ind(k,j,n)];
			}
			z[ind(i,j,n)] = sum;
			sum = 0;
		}
	}
}

/* Applies the Jacobi method onto the matrix g
Also updates the left and right singular values of the matrix g each time */
void jacobi(double *a, int n, double *s, double *u, double *v){
	double al, b, c, l, t, cs, sn, tmp, sign;
	/*i and j are the indices of the point we've chosen to zero out*/
	int i, j, p, q, k;
	identity(u,n);
	identity(v,n);

	/* Apply 10 times*/
	for (p = 0; p < 10; p++){
		/* for all pairs i<j */
		for (i=0; i<n; i++){
			for (j=i+1; j<n; j++){	
				al = b = c = l = t = cs = sn = tmp = sign = 0.0;
				/*Find the 2x2 submatrix*/
				for (k=0; k<n; k++){
					al += a[ind(k,i,n)] * a[ind(k,i,n)];
					b += a[ind(k,j,n)] * a[ind(k,j,n)];
					c += a[ind(k,i,n)] * a[ind(k,j,n)];
				}

				/*compute Jacobi rotation*/
				l = (b - al)/(2.0 * c);
				sign = 1.0;
				if (l < 0.0)
					sign = -1.0;
				t = sign / ((sign*l) + sqrt(1.0 + l*l));
				cs = 1.0/sqrt(1.0 + t*t);
				sn = cs *t;

				/* change columns i and j only*/
				for (k=0; k<n; k++){
					tmp = a[ind(k,i,n)];
					a[ind(k,i,n)] = cs*tmp - sn*a[ind(k,j,n)];
					a[ind(k,j,n)] = sn*tmp + cs*a[ind(k,j,n)];
				}

				/* update the right singular vectors */
				for (k=0; k<n; k++){
					tmp = v[ind(k,i,n)];
					v[ind(k,i,n)] = cs*tmp - sn*v[ind(k,j,n)];
					v[ind(k,j,n)] = sn*tmp + cs*v[ind(k,j,n)];
				}

			}
		}
	}

	/* find the singular values by adding up squares of each column, 
	then taking square root of each column */
	for (j=0; j<n; j++){
		for (i=0; i<n; i++){
			s[j] += a[ind(i,j,n)] * a[ind(i,j,n)];
		}	
		tmp = s[j];
		s[j] = sqrt(tmp);
	}
	
	/* sort the singular values largest to smallest, and the right matrix accordingly */
	for (p=0; p<(n-1); p++){
		for (j=0; j<n-p-1; j++){
			if (s[j] < s[j+1]){
				tmp = s[j];
				s[j] = s[j+1];
				s[j+1] = tmp;

				/* rearrange columns of u and v accordingly*/
				for (i=0; i<n; i++){
					tmp = v[ind(i,j,n)];
					v[ind(i,j,n)] = v[ind(i,j+1,n)];
					v[ind(i,j+1,n)] = tmp;
					tmp = a[ind(i,j,n)];
					a[ind(i,j,n)] = a[ind(i,j+1,n)];
					a[ind(i,j+1,n)] = tmp;
				}
			}
		}
	}

	/* a is A*V, so in order to get U, we divide by S in each column*/
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			a[ind(i,j,n)] = a[ind(i,j,n)]/s[j];
		}
	}

	/* Set U to A, since we have been making modifications to A */
	u=a;
}

int main()
{
	int n=4;
	double j[16]={0};
	double s[4] = {0};
	double u[16]={0};
	double v[16]={0};

	j[ind(0,0,n)] = 1.0;
	j[ind(0,1,n)] = 2.0;
	j[ind(1,0,n)] = 0.0;
	j[ind(1,1,n)] = 3.0;
	j[ind(1,2,n)] = 4.0;
	j[ind(1,3,n)] = 5.0;
	j[ind(2,0,n)] = 2.0;
	j[ind(2,1,n)] = 0.0;
	j[ind(2,2,n)] = 5.0;
	j[ind(2,3,n)] = 0.0;
	j[ind(3,0,n)] = 0.0;
	j[ind(3,1,n)] = 5.0;
	j[ind(3,2,n)] = 6.0;
	j[ind(3,3,n)] = 0.0;

	jacobi(j, n, s, u, v);
}
