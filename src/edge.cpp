/* Robert Hijmans, November 2011 

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"


SEXP _edge(SEXP d, SEXP dim, SEXP classes, SEXP type, SEXP directions) {

	R_len_t i, j;
	SEXP val;
	int nrow, ncol, n, cell, k;
	int *xd, *xval;

	int class = INTEGER(classes)[0];
	int edgetype = INTEGER(type)[0];
	int dirs = INTEGER(directions)[0];
	int falseval = 0;
	
	nrow = INTEGER(dim)[0];
	ncol = INTEGER(dim)[1];
	n = nrow * ncol;
	
	
	PROTECT(d = coerceVector(d, INTSXP));
	xd = INTEGER(d);

	PROTECT( val = allocVector( INTSXP, n) );
	xval = INTEGER(val);

	int r[8] = { -1,0,0,1 , -1,-1,1,1};
	int c[8] = { 0,-1,1,0 , -1,1,-1,1};	
	
	if (class == 0) {
		if (edgetype == 0) { // inner
			for (i = 1; i < (nrow-1); i++) {
				for (j = 1; j < (ncol-1); j++) {
					cell = i*ncol+j;
					xval[cell] = R_NaInt;
					if ( xd[cell] != R_NaInt ) {
						xval[cell] = falseval;
						for (k=0; k< dirs; k++) {
							if ( xd[cell + r[k] * ncol + c[k]] == R_NaInt ) {
								xval[cell] = 1;
								break;
							}
						}
					}
				}
			}
		
		} else if (edgetype == 1) { //outer
			for (i = 1; i < (nrow-1); i++) {
				for (j = 1; j < (ncol-1); j++) {
					cell = i*ncol+j;
					xval[cell] = falseval;
					if ( xd[cell] == R_NaInt ) {
						xval[cell] = R_NaInt;
						for (k=0; k < dirs; k++) {			
							if ( xd[cell+ r[k] * ncol + c[k] ] != R_NaInt ) {
								xval[cell] = 1;
								break;
							}
						}
					}
				}
			}
		} 
	} else { // by class
		int test;
		for (i = 1; i < (nrow-1); i++) {
			for (j = 1; j < (ncol-1); j++) {
				cell = i*ncol+j;
				test = xd[ cell+r[0]*ncol+c[0] ];
				if (test == R_NaInt) {
					xval[cell] = R_NaInt;			
				} else {
					xval[cell] = falseval;
				}
				for (k=1; k < dirs; k++) {
					if (test != xd[ cell+r[k]*ncol +c[k] ]) {
						xval[cell] = 1;
						break;
					}
				}
			}
		}

	}
	
	UNPROTECT(2);
	return(val);
}


*/