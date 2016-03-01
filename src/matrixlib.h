/**
 * Copyright (C) 2004 Harald SCHMIDBAUER - Vehbi Sinan TUNALIOGLU
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * A copy of the GNU General Public License is available via WWW at
 * http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 * writing to the Free Software Foundation, Inc., 59 Temple Place,
 * Suite 330, Boston, MA  02111-1307  USA.
 */

#include <math.h>

#ifndef MATRIXLIB
#define MATRIXLIB

/**
 * FUNCTION DEFINITIONS
 */

/* a macro for getting the maximum value of two given values */
#define maxval(p,q)				(p > q ? p : q)

/* memory allocation routine for a given type. quantity to be allocated is 1 */
#define NEW_MM(type)				((type *) calloc(1, sizeof(type)))

/* memory allocation routine for a given type. quantity to be allocated is user specified*/
#define NEW_A(quantity,type)	((type *) calloc((unsigned int)(quantity), sizeof(type)))

/**
 * set value of a matrix at the given row and column.
 * also check the boundries.
 */
#define m_set_val(A, i, j, val)							\
						((A)->me[i][j] = (				\
						(i) < 0 ||						\
						(j) < 0 ||						\
						(i) >= (A)->rows ||				\
						(j) >= (A)->cols ?				\
						err_exit(ERR_BOUNDS, "m_set_val"), 0.0 : (val)))

/**
 * set the m_foutput's file pointer to the stdout
 */
#define m_output(a)		(m_foutput(stdout, a))


/**
 * ERROR CODE DEFINITIONS
 */
#define ERR_DIMNEG				"Dimensions should be greater then 0"
#define ERR_MEM					"Cannot allocate memory for matrix initialization"
#define ERR_BOUNDS				"Setting value out of boundries"
#define ERR_SRCNULL				"Source matrix cannot be NULL"
#define ERR_DESTNULL			"Destination matrix should be pre-initialized"
#define ERR_DIFFDIMS			"Matrices differ in dimensions"
#define ERR_OPSNULL				"Operands cannot be NULL"
#define ERR_RESNULL				"Result matrix should be pre-initialized"
#define ERR_TRANSPDIFFDIMS		"Cannot transpose matrix. Dimensions don't match."
#define ERR_OPDIFFDIMS			"Cannot execute the operation. Dimensions don't match."
#define ERR_OPRESSAME			"Operands and result cannot be same"


/**
 * STRUCTURE DEFINITIONS
 */
typedef	struct
{
	unsigned int rows;	/* number of rows */
	unsigned int cols;	/* number of columns */
	double **me;		/* an array of dimension 'rows' by 'cols' */
} MAT;

int err_exit(char * err_message, char * func_name);
void mem_copy(int * source, int * destination, unsigned int length);
void m_transp(MAT * in, MAT * out);
void m_foutput(FILE * fp, MAT * a);
MAT * m_get(int m, int n);
void m_free(MAT * will_be_freed);
void m_copy(MAT * source, MAT * destination);
void m_add(MAT * mat1, MAT * mat2, MAT * result);
void m_mlt(MAT * A, MAT * B, MAT * OUT);
double m_det(MAT *mat);
double Determinant(double **a,int n);
void m_inverse(MAT *feed, MAT *result);
void CoFactor(double **a,int n,double **b);
void Transpose(double **a,int n);
void m_scalar_divide(MAT * arg, MAT * result, double divider);
void sm_mlt(double scalar, MAT * B, MAT * OUT);
#endif
