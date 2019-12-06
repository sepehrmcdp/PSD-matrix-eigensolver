/*******************************************************************************
Modified version of Eigenvalue solvers, tred2 and tqli, from "Numerical Recipes in C" (Cambridge
Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery
In this version matrices are linearized 1D arrays. Modifications by Sepehr Jalalian.
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))



void eig_dec_full(double* mat,double* eigs,double* eivs,double* e,int n){// mat is kx,e is temp
/*******************************************************************************
Output eigenvalues in eigs[n] and eigenvectors in eivs[n*n]. need e[n] as temp memory for calculations.
*******************************************************************************/

    //for(int i=0;i<n*n;i++){// copy
    //    eivs[i]=mat[i];
    //}
    for(int i=0;i<n;i++){
        for(int j=0; j<n;j++){
            eivs[i*n+j]=(mat[i*n+j]+mat[j*n+i])/2;
        }
    }
    tred2(eivs,n,eigs,e);
    tqli(eigs-1,e-1,n,eivs);




}


void eig_dec_lite(double* mat,double* eigs,double* eivs,double* e,int n){// mat is kx,e is temp
/*******************************************************************************
Output eigenvalues in eigs[n] without calculating eigenvectors. note that eivs[n*n] works as temp matrix for
calculations and should still be given as input. Only difference with full version is shorter runtime.
need e[n] as temp memory for calculations.
*******************************************************************************/


    //for(int i=0;i<n*n;i++){// copy
    //    eivs[i]=mat[i];
    //}
    for(int i=0;i<n;i++){
        for(int j=0; j<n;j++){
            eivs[i*n+j]=(mat[i*n+j]+mat[j*n+i])/2;
        }
    }
    tred2_lite(eivs,n,eigs,e);
    tqli_lite(eigs-1,e-1,n,eivs);




}

void tred2(double *a, int n, double d[], double e[])
/*******************************************************************************
Householder reduction of a real, symmetric matrix a[1..n][1..n].
On output, a is replaced by the orthogonal matrix Q effecting the
transformation. d[1..n] returns the diagonal elements of the tridiagonal matrix,
and e[1..n] the off-diagonal elements, with e[1]=0. Several statements, as noted
in comments, can be omitted if only eigenvalues are to be found, in which case a
contains no useful information on output. Otherwise they are to be included.
*******************************************************************************/
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n-1;i>=1;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<=l;k++)
				scale += fabs(a[i*n+k]);
			//printf("reach0");
			if (scale == 0.0) /* Skip transformation. */
				e[i]=a[i*n+l];
			else {
				for (k=0;k<=l;k++) {
					a[i*n+k] /= scale; /* Use scaled a's for transformation. */
					h += a[i*n+k]*a[i*n+k]; /* Form sigma in h. */
				}
				f=a[i*n+l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g; /* Now h is equation (11.2.4). */
				a[i*n+l]=f-g; /* Store u in the ith row of a. */
				f=0.0;
				for (j=0;j<=l;j++) {
					/* Next statement can be omitted if eigenvectors not wanted */
					a[j*n+i]=a[i*n+j]/h; /* Store u/H in ith column of a. */
					g=0.0; /* Form an element of A.u in g. */
					for (k=0;k<=j;k++)
						g += a[j*n+k]*a[i*n+k];
					for (k=j+1;k<=l;k++)
						g += a[k*n+j]*a[i*n+k];
					e[j]=g/h; /* Form element of p in temporarily unused element of e. */
					f += e[j]*a[i*n+j];
				}
				hh=f/(h+h); /* Form K, equation (11.2.11). */
				for (j=0;j<=l;j++) { /* Form q and store in e overwriting p. */
					f=a[i*n+j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<=j;k++) /* Reduce a, equation (11.2.13). */
						a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
				}
			}
		} else
			e[i]=a[i*n+l];
		d[i]=h;
		}
		//printf("reach1");
		/* Next statement can be omitted if eigenvectors not wanted */
		d[0]=0.0;
		e[0]=0.0;
		/* Contents of this loop can be omitted if eigenvectors not
		   wanted except for statement d[i]=a[i][i]; */
		for (i=0;i<n;i++) { /* Begin accumulation of transformation matrices. */
			l=i-1;
		if (d[i]) { /* This block skipped when i=1. */
			for (j=0;j<=l;j++) {
				g=0.0;
				for (k=0;k<=l;k++) /* Use u and u/H stored in a to form P.Q. */
					g += a[i*n+k]*a[k*n+j];
				for (k=0;k<=l;k++)
					a[k*n+j] -= g*a[k*n+i];
			}
		}
		d[i]=a[i*n+i]; /* This statement remains. */
		a[i*n+i]=1.0; /* Reset row and column of a to identity matrix for next iteration. */
		for (j=0;j<=l;j++) a[j*n+i]=a[i*n+j]=0.0;
		//printf("reach2");
	}
}

void tqli(double d[], double e[], int n, double *z)
/*******************************************************************************
QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix
previously reduced by tred2 sec. 11.2. On input, d[1..n] contains the diagonal
elements of the tridiagonal matrix. On output, it returns the eigenvalues. The
vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix, with
e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues,
several lines may be omitted, as noted in the comments. If the eigenvectors of
a tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as the
identity matrix. If the eigenvectors of a matrix that has been reduced by tred2
are required, then z is input as the matrix output by tred2. In either case,
the kth column of z returns the normalized eigenvector corresponding to d[k].
*******************************************************************************/
{
	double pythag(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i]; /* Convenient to renumber the elements of e. */
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) { /* Look for a single small subdiagonal element to split the matrix. */
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) printf("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]); /* Form shift. */
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); /* This is dm - ks. */
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { /* A plane rotation as in the original QL, followed by Givens */
					f=s*e[i];          /* rotations to restore tridiagonal form.                     */
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) { /* Recover from underflow. */
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted */
					for (k=1;k<=n;k++) { /* Form eigenvectors. */
						f=z[(k-1)*n+i+1-1];
						z[(k-1)*n+i+1-1]=s*z[(k-1)*n+i-1]+c*f;
						z[(k-1)*n+i-1]=c*z[(k-1)*n+i-1]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
	printf("one done");
}
void tred2_lite(double *a, int n, double d[], double e[])
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n-1;i>=1;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<=l;k++)
				scale += fabs(a[i*n+k]);
			//printf("reach0");
			if (scale == 0.0) /* Skip transformation. */
				e[i]=a[i*n+l];
			else {
				for (k=0;k<=l;k++) {
					a[i*n+k] /= scale; /* Use scaled a's for transformation. */
					h += a[i*n+k]*a[i*n+k]; /* Form sigma in h. */
				}
				f=a[i*n+l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g; /* Now h is equation (11.2.4). */
				a[i*n+l]=f-g; /* Store u in the ith row of a. */
				f=0.0;
				for (j=0;j<=l;j++) {
					/* Next statement can be omitted if eigenvectors not wanted */
					//a[j*n+i]=a[i*n+j]/h; /* Store u/H in ith column of a. */
					g=0.0; /* Form an element of A.u in g. */
					for (k=0;k<=j;k++)
						g += a[j*n+k]*a[i*n+k];
					for (k=j+1;k<=l;k++)
						g += a[k*n+j]*a[i*n+k];
					e[j]=g/h; /* Form element of p in temporarily unused element of e. */
					f += e[j]*a[i*n+j];
				}
				hh=f/(h+h); /* Form K, equation (11.2.11). */
				for (j=0;j<=l;j++) { /* Form q and store in e overwriting p. */
					f=a[i*n+j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<=j;k++) /* Reduce a, equation (11.2.13). */
						a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
				}
			}
		} else
			e[i]=a[i*n+l];
		d[i]=h;
		}
		//printf("reach1");
		/* Next statement can be omitted if eigenvectors not wanted */
		d[0]=0.0;
		e[0]=0.0;
		/* Contents of this loop can be omitted if eigenvectors not
		   wanted except for statement d[i]=a[i][i]; */
		for (i=0;i<n;i++) { /* Begin accumulation of transformation matrices. */
		//	l=i-1;
		//if (d[i]) { /* This block skipped when i=1. */
		//	for (j=0;j<=l;j++) {
		//		g=0.0;
		//		for (k=0;k<=l;k++) /* Use u and u/H stored in a to form P.Q. */
		//			g += a[i*n+k]*a[k*n+j];
		//		for (k=0;k<=l;k++)
		//			a[k*n+j] -= g*a[k*n+i];
		//	}
		//}
		d[i]=a[i*n+i]; /* This statement remains. */
		//a[i*n+i]=1.0; /* Reset row and column of a to identity matrix for next iteration. */
		//for (j=0;j<=l;j++) a[j*n+i]=a[i*n+j]=0.0;
		//printf("reach2");
	}
}

void tqli_lite(double d[], double e[], int n, double *z)
{
	double pythag(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i]; /* Convenient to renumber the elements of e. */
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) { /* Look for a single small subdiagonal element to split the matrix. */
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) printf("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]); /* Form shift. */
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); /* This is dm - ks. */
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { /* A plane rotation as in the original QL, followed by Givens */
					f=s*e[i];          /* rotations to restore tridiagonal form.                     */
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) { /* Recover from underflow. */
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted */
					//for (k=1;k<=n;k++) { /* Form eigenvectors. */
					//	f=z[(k-1)*n+i+1-1];
					//	z[(k-1)*n+i+1-1]=s*z[(k-1)*n+i-1]+c*f;
					//	z[(k-1)*n+i-1]=c*z[(k-1)*n+i-1]-s*f;
					//}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
	printf("one done");
}


double pythag(double a, double b)
/*******************************************************************************
Computes (a2 + b2)1/2 without destructive underflow or overflow.
*******************************************************************************/
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

