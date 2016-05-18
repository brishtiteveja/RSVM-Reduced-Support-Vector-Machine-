#include "svm.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* LEVEL 1 BLAS */
extern double ddot_(int *, double *, int *, double *, int *);
/* LEVEL 2 BLAS */
extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
/* MINPACK 2 */
extern void dtron(int, double *, double *, double *, double, double, double, double, int, double);

int nfev;
int inc = 1;
double one = 1, zero = 0;
double ALPHA = 5;
int l=100;
double *Q=NULL;
double ga=1;

double getNorm1(int n, double *g){
	double d=0;
	int i;
	double j;
	for(i=0;i<n;i++)
		if((j =fabs(g[i]))>d)
			d = j;
	return d;
}

void newvBp(double **vB, double **p) {
	printf("ALPHA=%f\n", ALPHA);
	*vB = (double *)malloc(sizeof(double)*l);
	*p = (double *)malloc(sizeof(double)*l);
}

void freevBp(double *vB, double *p) {
	free(vB);
	free(p);
}

int vavb(int m, double *x, double *vA, double *vB)
{
	int i;

/* 1st step: evaluate vA = g0-Q*x (dimQ=l*n;dimx=dimb=n(=mbar+1)) */
/* 2nd step: evaluate vB = exp(-ALPHA*vA) */

	for(i=0;i<l;++i) {
		vA[i] = 1-ddot_(&m, &(Q[i]), &l, x, &inc);
		vB[i] = exp(-vA[i]*ALPHA);
	}

	return 0;
}

int upv(int m, double *x, double *vB, double *p)
{
	/* first call vavb */
	int i;
	double *vA = (double *) malloc(sizeof(double)*l);
	vavb(m,x,vA,vB);

	/* evaluate the function vector p = log(1+vB)/ALPHA + vA  */
	for(i=0;i<l;++i)
		p[i]=log(1+vB[i])/ALPHA + vA[i];

	free(vA);
	return 0;
}

int uhes(int n, double *p, double *vB, double *h)
{
	/* evaluate hession Hij=gamma*kSigma{Qik*Qjk*(1+ALPHA*vB[k]*p[k])/(1+vB[k])^2}+delta(ij) */
	double *t = (double *) malloc(sizeof(double)*l);
	int i,j,k;
	for(k=0;k<l;++k)
		t[k] = (1+ALPHA*vB[k]*p[k])/(1+vB[k])/(1+vB[k]);	

	for(i=0;i<n;++i)
		for(j=0;j<n;++j) {
			h[i*n+j]= (i==j)?1/ga:0;
			for(k=0;k<l;++k)
				h[i*n+j] += Q[i*l+k]*Q[j*l+k]*t[k]; 
		}

	free(t); 
	return 0;
}

int ugrad(int m, double *x, double *p, double *vB, double *g)
{
	int i;
	double *t = (double *) malloc(sizeof(double)*l);

	/* evaluate gradian g=-ga*Q'*[p/(1+vB)]+x  (dimQ'=n*l;dimvB=dimp=l) */
	memcpy(g, x, sizeof(double)*m);

	for(i=0;i<l;++i)
		t[i]=p[i]/(1+vB[i]);
	
	for(i=0;i<m;++i) {
		g[i]/=ga;
		g[i]-= ddot_(&l, &(Q[i*l]), &inc, t, &inc); 
	}

	free(t);	
	return 0;
}

int ufv(int m, double *x, double *p, double *f)
{
	*f = (ddot_(&l, p, &inc, p, &inc) + ddot_(&m, x, &inc, x,&inc)/ga)/2;
	return ++nfev;
}
/* return the objective value */

void solvebqp(struct BQP *qp, double alfa)
{
	/* driver for positive semidefinite quadratic programing version
	of tron */
	int i, n, maxfev;
	double *x, *xl, *xu;
	double frtol, fatol, fmin, gtol, cgtol;

	ALPHA = alfa;
	n = qp->n;
	l = qp->l;
	ga = qp->C[0];
	maxfev = 1000; /* ? */
	nfev = 0;

	x = qp->x;
	xu = (double *) malloc(sizeof(double)*n);
	Q = qp->Q;
	xl = (double *) malloc(sizeof(double)*n);
	for (i=0;i<n;i++) {
		xu[i] = 1e+20;
		xl[i] = -xu[i];
	}

	fatol = 0;
	frtol = 1e-12;
	fmin = -1e+32;
	cgtol = 0.1;	
	gtol = qp->eps;

	dtron(n, x, xl, xu, gtol, frtol, fatol, fmin, maxfev, cgtol);

	free(xl);
	
}
