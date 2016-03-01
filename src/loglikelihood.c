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

#include <R.h>
#include <math.h>
#include "matrixlib.h"
#include "loglikelihood.h"

/**
 * count_triangular counts the nonzero elements of the matrix.
 * example:
 * for a 2x2 upper triangular matrix:
 * 1 1
 * 0 1
 * ==> 3 'nonzero' element
 * for a 3x3 upper triangular matrix:
 * 1 1 1
 * 0 1 1
 * 0 0 1
 * ==> 6 'nonzero' elements
 */
int count_triangular(int dimension)
{
	if(dimension <= 0)
	{
		return 0;
	}
	else
	{
		return dimension + count_triangular(dimension - 1);
	}
}
/**
 * initializes the C matrix
 */
void initialize_C(double * par_array, MAT *C)
{
	int i, j, count;
	count = 0;
	for(i = 0; i < C->rows; i++)
	{
		for(j = 0; j < C->cols; j++)
		{
			if(i >= j)
			{
				m_set_val(C, j, i, par_array[count]);
				count++;
			}
		}
	}


}

/**
 * initializes the A and G matrices.
 */
void initialize_AG(double *par_array, int series_count, MAT *A, int count)
{
	int i, j, num;
	num = 0;
	for(i = 0; i < A->rows; i++)
	{
		for(j = 0; j < A->cols; j++)
		{
			m_set_val(A, j, i, par_array[count * series_count * series_count + num + count_triangular(series_count)]);
			num++;
		}
	}
}

/**
 * initialized the initial H matrix.
 */
void initialize_H(int dimension, MAT *H)
{
	/**
	 * TODO:
	 * ask what this H matrix should look like
	 */
	int i, j;
	for(i = 0; i < dimension; i++)
	{
		for(j = 0; j < dimension; j++)
		{
			if(i == j)
			{
				m_set_val(H, j, i, 1);
			}
			else
			{
				m_set_val(H, j, i, 0);
			}
		}

	}
}

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
	)
{
	/**
	 * TODO
	 * all malloced things should be freed at the end.
	 */

	/**
	 * VARIABLE DECLARATIONS
	 */

	/* Paramter Matrices */
	MAT *C = NULL;			/* declare the C parameter matrix*/
	MAT *C_t = NULL;		/* declare the transposed C parameter matrix*/
	MAT *C_term = NULL;		/* C'C term */
	MAT **A = NULL, **G = NULL;	/* declare the arch and garch parameter matrix arrays*/
	MAT **A_t = NULL, **G_t = NULL;	/* declare the transposed arch and garch parameter matrix arrays*/

	/* H Matrices */
	MAT *H = NULL;			/* H term */
	MAT **HOLD = NULL;		/* an array holding former H terms */

	/* temporary matrices used for calculations */

	MAT *E0 = NULL;
	MAT *E1 = NULL;
	MAT *E1_tmp = NULL;
	MAT *E2 = NULL;
	MAT *ETERM = NULL;
	MAT *TEMP  = NULL;
	MAT *ATERM  = NULL;
	MAT *GTERM  = NULL;

	int count, counttemp, i, j;	/* will be used for loop counter */
	double temp, buffer;		/* will be used for temporary calculations */
	double detcomp;			/* the inner part of sqrt() function for determinant calculation */

	double * buffer_params;
	int shifted;
	int temp_length;
	int check_point;

	/**
	 * get the seriestl series and convert it into a multidimensional series array
	 */
	double **series;
	series = (double **)malloc(*series_count * sizeof(double *));
	count = 0;

	for(i = 0; i < *series_count; i++)
	{
		series[i] = (double *)malloc(*series_length * sizeof(double));
	}
	for(j = 0; j < *series_length; j++)
	{
		for(i = 0; i < *series_count; i++)
		{
			series[i][j] = seriestl[count];
			count++;
		}
	}
/*
	for(j = 0; j < *series_length; j++)
	{
		for(i = 0; i < *series_count; i++)
		{
			Rprintf("%f - ", series[i][j]);
		}
		Rprintf("\n");
	}
*/


	/**
	 * TODO
	 * check all the function arguments
	 */

	/**
	 * TODO
	 * check the model specification,
	 * whether p or q can be `0' or not
	 *
	 * check the switcher values of BEKK model.
	 * these are p and q values of BEKK(p,q,K) model.
	 * These can not be less then 0
	 */

	if(switcher[0] < 0)	/* check p switch */
	{
		Rprintf("ERROR! invalid p switch. ( p = %d )\n", switcher[0]);
		*returnval = ERROR_TERM;
		return;
	}

	if(switcher[1] < 0)	/* check q switch */
	{
		Rprintf("ERROR! invalid q switch. ( q = %d )\n", switcher[1]);
		*returnval = ERROR_TERM;
		return;
	}
	/**
	 * initialize the variables
	 */

	/* memory allocation for the parameter matrix arrays */
	A   = (MAT **) malloc(switcher[1] * sizeof(MAT *));
	A_t = (MAT **) malloc(switcher[1] * sizeof(MAT *));
	G   = (MAT **) malloc(switcher[0] * sizeof(MAT *));
	G_t = (MAT **) malloc(switcher[0] * sizeof(MAT *));

	/* memory allocation for the array that holds the previous H terms */
	HOLD   = (MAT **) malloc(switcher[0] * sizeof(MAT *));

	/**
	 * ATTENTION!
	 * HERE, WE WILL INSERT THE FIXED PARAMETERS
	 */

	if(*fixed_length)
	{
		shifted = 0;
		/**
		 * TODO
		 * ask the C matrix??
		 */
		temp_length = (count_triangular(*series_count) + *series_count * *series_count * switcher[0] + *series_count * *series_count * switcher[1]);
		buffer_params = (double *) malloc(temp_length * sizeof(double));

		for(count = 0; count < temp_length; count++)
		{
			check_point = 0;
			for(counttemp = 0; counttemp < *fixed_length; counttemp++)
			{
				if(count == (fixed_indexs[counttemp] - 1))
				{
					check_point = 1;
					shifted++;
					buffer_params[count] = fixed_values[counttemp];
					break;
				}
			}
			if(check_point == 0)
			{
				buffer_params[count] = params[count - shifted];
			}
		}

		/**
		 * TODO
		 * learn! what about the params array?
		 * is the stack automatically freed although we refer to another place
		 * with params symbol?
		 * Anyway, no big deal...
		 */
		params = buffer_params; /* copy the modified parameter array into params */
	}

	/* initialize the C parameter matrix*/
	C = m_get(*series_count, *series_count);
	initialize_C(params, C);
	/* m_output(C); */

	/* initialize the transposed C parameter matrix*/
	C_t = m_get(*series_count, *series_count);
	m_transp(C, C_t);

	/* initialize the C term by multiplying the C and C_t matrices */
	C_term = m_get(*series_count, *series_count);
	m_mlt(C_t, C, C_term);

	/* We have the C term. Now, do the critical A and G parameter matrices allocation*/
	for(count = 0; count < switcher[1]; count++)
	{
		A[count] = m_get(*series_count, *series_count);
		initialize_AG(params, *series_count, A[count], count);
		/* m_output(A[count]); */

		A_t[count] = m_get(*series_count, *series_count);
		m_transp(A[count], A_t[count]);
	}

	counttemp = count;	/* keep the counter where we left in the params matrix */

	for(count = 0; count < switcher[0]; count++)
	{
		G[count] = m_get(*series_count, *series_count);
		initialize_AG(params, *series_count, G[count], counttemp + count);
		/* m_output(G[count]); */

		G_t[count] = m_get(*series_count, *series_count);
		m_transp(G[count], G_t[count]);
	}

	/* initialize the H matrix. Then, do the critical buffer HOLD array for previous H terms */
	H = m_get(*series_count, *series_count);
	initialize_H(*series_count, H);
	/* m_output(H); */

	for(count = 0; count < switcher[0]; count++)
	{
		HOLD[count] = m_get(*series_count, *series_count);
	}

	/**
	 * TEMPORARY MATRICES FOR CALCULATIONS
	 */
	E0	= m_get(1, 1); /* an 1 by 1 matrix */
	E1	= m_get(*series_count, 1); /* an n by 1 matrix */
	E1_tmp	= m_get(*series_count, 1); /* an n by 1 matrix */
	E2	= m_get(1, *series_count); /* an 1 by n matrix (transpose of the E1 matrix) */
	ETERM	= m_get(*series_count, *series_count); /* an n by n matrix (E1 * E2) */
	TEMP	= m_get(*series_count, *series_count); /* an n by n matrix for temp usage */
	ATERM	= m_get(*series_count, *series_count); /* an n by n matrix */
	GTERM	= m_get(*series_count, *series_count); /* an n by n matrix */

	/**
	 * now doing the calculations
	 */

	*returnval = 0.0;					/* reset the return value */
	/**
	 * ATTENTION:
	 * we will now set the counter to its initial value.
	 * it is the max of the two switchers
	 * count = max(p,q)
	 */


	count = maxval(switcher[0], switcher[1]);
	while(count < *series_length) /* begin the loop */
	{
		/**
		 * CRITICAL!
		 * shift the H TERMS to previous HOLD items.
		 * HOLD[1] becomes H, HOLD[2] becomes HOLD[1] and so on
		 */

		for(counttemp = switcher[0] - 1; counttemp > 0; counttemp--)
		{
			m_copy(HOLD[counttemp - 1], HOLD[counttemp]);
		}

		/**
		 * apply if GARCH parameter is greater then zero
		 */
		if(switcher[0] > 0)
		{
			m_copy(H, HOLD[0]);
		}

		/**
		 * a bit complicated but following explanation will be useful hopefully
		 * H = (C')x(C) + (A')(E_t-1)(E_t-1')(A) + (B')(E_t-2)(E_t-2')(B) + ... +  (G')(H_t-1)(G) + (L')(H_t-2)(L) + ...
		 *                    |_____________|          |_____________|             |____________|   |____________| |_____|
		 *                        E1 TERM                  E2 TERM                     G1 TERM         G2 TERM     G3.G4..
		 *                |____________________|   |____________________| |_____|
		 *                        A1 TERM                  A2 TERM        A3.A4..
		 *     |______|  |_____________________________________________________|  |______________________________________|
		 *      C TERM                         A TERM                                              G TERM
		 */

		m_copy(C_term, H);

		/**
		 * first, do the A TERM part...
		 */
		for(counttemp = 0; counttemp < switcher[1]; counttemp++)
		{
			/**
			 * calculate the (E_t-x %*% E_t-x') term
			 */
			/**
			 * first initialize the E matrix
			 */
			for(i = 0; i < *series_count; i++)
			{
				m_set_val(E1, i, 0, series[i][count - 1 - counttemp]);
			}
			m_transp(E1, E2);
			m_mlt(E1, E2, ETERM); /* we got the term in ETERM*/

			/**
			 * calculate the (A_x' %*% E_t-x %*% E_t-x' %*% A_x) term
			 */

			m_mlt(A_t[counttemp], ETERM, TEMP);
			m_mlt(TEMP, A[counttemp], ATERM); /* we got the term in ATERM */

			/**
			 * accumulate the ATERM term into H
			 */

			m_add(H, ATERM, TEMP); /* addition completed, but assign the TEMP back to H */
			m_copy(TEMP, H);
		}

		/**
		 * A term calculated and accumulated to H term.
		 * Continue with the G term
		 */

		for(counttemp = 0; counttemp < switcher[0]; counttemp++)
		{
			/**
			 * calculate the (G_x' %*% H_x %*% G_x) term
			 */

			m_mlt(G_t[counttemp], HOLD[counttemp], TEMP);
			m_mlt(TEMP, G[counttemp], GTERM); /* we got the term in GTERM */

			/**
			 * accumulate the GTERM term into H
			 */

			m_add(H, GTERM, TEMP); /* addition completed, but assign the TEMP back to H */
			m_copy(TEMP, H);
		}

		/**
		 * TODO THIS PART OF THE CODE IS LEFT FROM THE PREVIOUS VERSION
		 * we got the H TERM.
		 *
		 * CRITICAL PART
		 * multivariate density normal distribution (???)
		 *
		 * MULTIVARIATE NORMAL DENSITY FUNCTION
		 *
		 * f(x) = (1 / (2*pi)^(n/2) * |E|^(1/2)) * (e^( -(x - m)' * E^(-1) * (x - m) / 2))
		 *
		 * MULTIVARIATE LOGNORMAL DENSITY FUNCTION
		 *
		 * log[f(x)] = - log[(2 * pi)^(n/2) *  sqrt(|E|)] - (x - m)' * E^(-1) * (x - m) / 2
		 *                                                       |____________ ___________|
		 *                                                                    V
		 *                                                        this will be stored in temp
		 *
		 */

		/**
		 * CALCULATING TEMP VALUE
		 */
		temp = 0;	/* reset the temp variable */

		/**
		 * CRITICAL!
		 * check whether the det(H) is 0 or not.
		 * if it is 0, that means that the H term is not invertable.
		 * if it is not invertable, we will add following matrix to the
		 * H term: (CHANGED: WE ONLY MULTIPLY THE H[0,0] WITH 1.01)
		 *
		 * NEW_H_TERM = H * | 0.01  0.00 |
		 *                  | 0.00  0.00 |
		 */
		detcomp = m_det(H);

		/**
		 * DEBUG
		 * printf("det(H)² = %f\n", detcomp);
		 */

		if(detcomp == 0)
		{
			/**
			 * H is not invertable
			 */
			Rprintf("H IS SINGULAR!...\n");

			H->me[0][0] *= 1.01;
			detcomp = m_det(H);
		}
		else if(detcomp < 0)
		{
			/**
			 * it should not be negative... else sqrt(detcomp) = undefined...
			 */
			*returnval = VERYBIGNUMBER;
			return;
		}


		/**
		 * calculate  E^(-1)
		 * TODO
		 * that doesn't work sometimes. Thus we do it by hand:
		 * m_pow(H, -1, T3);
		 */

		m_inverse(H, TEMP);

		/**
		 * a loop here
		 */
		for(i = 0; i < *series_count; i++)
		{
			/* initializing (x - m) term (since m = 0; x)*/
			m_set_val(E1, i, 0, series[i][count]);
		}

		m_transp(E1, E2);

		m_mlt(TEMP, E1, E1_tmp);
		m_mlt(E2, E1_tmp, E0); /* ??? */

		temp = E0->me[0][0];

		buffer = (-1) * log(2 * M_PI * sqrt(detcomp)) - (temp / 2);
		*returnval += buffer;

		count++;					/* next period */
	}

	*returnval *= -1;
	if(isnan(*returnval))
	{
		*returnval = VERYBIGNUMBER;
	}
	/**
	 * DEBUG
	 * Rprintf("\nreturnval : %f\n", *returnval);
	 */
	/**
	 * free malloced variables
	 */
	for(i = 0; i < *series_count; i++)
	{
		free(series[i]);
	}
	free(series);

	for(i = 0; i < switcher[1]; i++)
	{
		m_free(A[i]);
		m_free(A_t[i]);
	}
	free(A);
	free(A_t);

	for(i = 0; i < switcher[0]; i++)
	{
		m_free(G[i]);
		m_free(G_t[i]);
		m_free(HOLD[i]);
	}
	free(G);
	free(G_t);
	free(HOLD);

	m_free(C);
	m_free(C_t);
	m_free(C_term);
	m_free(H);
	m_free(E0);
	m_free(E1);
	m_free(E1_tmp);
	m_free(E2);
	m_free(ETERM);
	m_free(TEMP);
	m_free(ATERM);
	m_free(GTERM);

}
