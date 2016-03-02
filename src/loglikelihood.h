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

#include <stdio.h>
#include <math.h>
#include "matrixlib.h"

#define VERYBIGNUMBER  1.0e+10
#define ERROR_TERM     1.0e+10

#define maxval(p,q)  (p > q ? p : q)
int count_triangular(int dimension);
void initialize_C(double * par_array, MAT *C);
void initialize_AG(double *par_array, int series_count, MAT *A, int count);
void initialize_H(int dimension, MAT *H);
void loglikelihood(
			double *params,		/* the parameter list */
			int    *fixed_indexs,	/* index list for fixed parameters */
			double *fixed_values,	/* value list for fixed parameters */
			int    *fixed_length,	/* length fixed parameters */
			double *seriestl,	/* series array */
			int    *series_count,	/* number of series */
			int    *series_length,	/* length of the series */
			int    *switcher,	/* the switcher array for BEKK's p and q */
			double *returnval	/* value to be returned at the end */
	);
