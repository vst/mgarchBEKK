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
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"

#define VERYBIGNUMBER  1.0e+10
#define ERROR_TERM     1.0e+10

#define maxval(p,q)  (p > q ? p : q)

void loglikelihood_GJR(
					double *params,			/* the parameter list*/
					int    *fixed_indexs,	/* index list for fixed parameters	*/
					double *fixed_values,	/* value list for fixed parameters */
					int    *fixed_length,	/* length fixed parameters */
					double *series_1,		/* first series	*/
					double *series_2,		/* second series */
					int    *series_length,	/* second series */
					int    *switcher,		/* the switcher array for BEKK's p and q */
					double *returnval		/* value to be returned at the end */
				)
{
	/**
	 * VARIABLE DECLARATIONS
	 */

	/* Paramter Matrices */
	MAT *C = NULL;					/* declare the C parameter matrix*/
	MAT *C_t = NULL;				/* declare the transposed C parameter matrix*/
	MAT *C_term = NULL;				/* C'C term */
	MAT **A = NULL, **G = NULL;		/* declare the arch and garch parameter matrix arrays*/
	MAT **A_t = NULL, **G_t = NULL;	/* declare the transposed arch and garch parameter matrix arrays*/

	MAT **R = NULL, **R_t = NULL;	/* declare the GJR parameter matrix arrays*/

	/* H Matrices */
	MAT *H = NULL;					/* H term */
	MAT **HOLD = NULL;				/* an array holding former H terms */

	/* temporary matrices used for calculations */
	MAT *T0 = NULL;							/* 1x1 matrices */
	MAT *T1 = NULL, *T7 =NULL;				/* 2x1 matrices */
	MAT *T3 = NULL, *T4 = NULL, *T5 = NULL;	/* 2x2 matrices */
	MAT *T6 = NULL, *T2 = NULL, *T8 = NULL;	/* 1x2 matrices */

	int count, counttemp, counttemp2;	/* will be used for loop counter */
	double temp, buffer;				/* will be used for temporary calculations */
	double detcomp;						/* the inner part of sqrt() function for determinant calculation */
	double w;							/* last parameter for the GJR component */
	double S; 								/* GJR pareto, 0 or 1 */

	double * buffer_params;
	int shifted;
	int temp_length;
	int check_point;

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

	if(switcher[2] < 0)	/* check g switch */
	{
		Rprintf("ERROR! invalid g switch. ( g = %d )\n", switcher[2]);
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
	R   = (MAT **) malloc(switcher[2] * sizeof(MAT *));
	R_t = (MAT **) malloc(switcher[2] * sizeof(MAT *));

	/* memory allocation for the array that holds the previous H terms */
	HOLD   = (MAT **) malloc(switcher[0] * sizeof(MAT *));

	/**
	 * ATTENTION!
	 * HERE, WE WILL INSERT THE FIXED PARAMETERS
	 */

	if(*fixed_length)
	{
		shifted = 0;
		temp_length = (3 + 4 * switcher[0] + 4 * switcher[1] + 4 * switcher[2] + 1);
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

		params = buffer_params; /* copy the modified parameter array into params*/
	}

	/* initialize the C parameter matrix*/
	C = m_get(2, 2);				/* make a 2 by 2 matrix */
	m_set_val(C, 0, 0, params[0]);
	m_set_val(C, 1, 0, 0);
	m_set_val(C, 0, 1, params[1]);
	m_set_val(C, 1, 1, params[2]);

	/* initialize the transposed C parameter matrix*/
	C_t = m_get(2, 2);				/* make a 2 by 2 matrix */
	m_transp(C, C_t);				/* transpose the matrix */

	/* initialize the C term by multiplying the C and C_t matrices */
	C_term = m_get(2, 2);			/* make a 2 by 2 matrix */
	m_mlt(C_t, C, C_term);			/* initialize the C'C term */

	/* We have the C term. Now, do the critical A and G parameter matrices allocation*/
	for(count = 0; count < switcher[1]; count++)
	{
		A[count] = m_get(2, 2);				/* make a 2 by 2 matrix */
		m_set_val(A[count], 0, 0, params[3 + (4 * count)]);
		m_set_val(A[count], 1, 0, params[4 + (4 * count)]);
		m_set_val(A[count], 0, 1, params[5 + (4 * count)]);
		m_set_val(A[count], 1, 1, params[6 + (4 * count)]);

		A_t[count] = m_get(2, 2);			/* make a 2 by 2 matrix */
		m_transp(A[count], A_t[count]);		/* transpose the matrix */
	}

	counttemp = count;	/* keep the counter where we left in the params matrix */

	for(count = 0; count < switcher[0]; count++)
	{
		G[count] = m_get(2, 2);				/* make a 2 by 2 matrix */
		m_set_val(G[count], 0, 0, params[3 + (4 * count) + (4 * counttemp)]);
		m_set_val(G[count], 1, 0, params[4 + (4 * count) + (4 * counttemp)]);
		m_set_val(G[count], 0, 1, params[5 + (4 * count) + (4 * counttemp)]);
		m_set_val(G[count], 1, 1, params[6 + (4 * count) + (4 * counttemp)]);

		G_t[count] = m_get(2, 2);			/* make a 2 by 2 matrix */
		m_transp(G[count], G_t[count]);		/* transpose the matrix */
	}

	counttemp2 = counttemp;	/* keep the counter where we left in the params matrix */


	for(count = 0; count < switcher[2]; count++)
	{
		R[count] = m_get(2, 2);				/* make a 2 by 2 matrix */
		m_set_val(R[count], 0, 0, params[3 + (4 * count) + (4 * counttemp) + (4 * counttemp2)]);
		m_set_val(R[count], 1, 0, params[4 + (4 * count) + (4 * counttemp) + (4 * counttemp2)]);
		m_set_val(R[count], 0, 1, params[5 + (4 * count) + (4 * counttemp) + (4 * counttemp2)]);
		m_set_val(R[count], 1, 1, params[6 + (4 * count) + (4 * counttemp) + (4 * counttemp2)]);

		R_t[count] = m_get(2, 2);			/* make a 2 by 2 matrix */
		m_transp(R[count], R_t[count]);		/* transpose the matrix */
	}

	/* initialize the H matrix. Then, do the critical buffer HOLD array for previous H terms */
	H = m_get(2, 2);				/* make a 2 by 2 matrix */
	m_set_val(H, 0, 0, 1);
	m_set_val(H, 0, 1, 0);
	m_set_val(H, 1, 0, 0);
	m_set_val(H, 1, 1, 1);

	for(count = 0; count < switcher[0]; count++)
	{
		HOLD[count] = m_get(2, 2);				/* make a 2 by 2 matrix */
		m_set_val(HOLD[count], 0, 0, 1);
		m_set_val(HOLD[count], 1, 0, 0);
		m_set_val(HOLD[count], 0, 1, 0);
		m_set_val(HOLD[count], 1, 1, 1);
	}

	/* initialize the temporary matrices */
	T0 = m_get(1, 1);				/* make a 1 by 1 matrix */

	T1 = m_get(2, 1);				/* make a 2 by 1 matrix */
	T2 = m_get(1, 2);				/* make a 1 by 2 matrix */

	T3 = m_get(2, 2);				/* make a 2 by 2 matrix */
	T4 = m_get(2, 2);				/* make a 2 by 2 matrix */
	T5 = m_get(2, 2);				/* make a 2 by 2 matrix */

	T6 = m_get(1, 2);				/* make a 1 by 2 matrix */
	T7 = m_get(2, 1);				/* make a 2 by 1 matrix */
	T8 = m_get(1, 2);				/* make a 1 by 2 matrix */

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
		 * here is the BEKK model. We add a further operand: R TERM.
		 * it is similar to the A TERM, but each R TERM is multiplied with an S value.
		 * S value is either 0 or 1.
		 *
		 * If (w)(e[1][t-1]) + (1-w)([e]2[t-1]) > 0, S = 1 otherwise S = 0
		 *
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

			m_set_val(T1, 0, 0, series_1[count - 1 - counttemp]);
			m_set_val(T1, 1, 0, series_2[count - 1 - counttemp]);
			m_transp(T1, T2);
			m_mlt(T1, T2, T3); /* we got the term in T3*/

			/**
			 * calculate the (A_x' %*% E_t-x %*% E_t-x' %*% A_x) term
			 */

			m_mlt(A_t[counttemp], T3, T4);
			m_mlt(T4, A[counttemp], T3); /* we got the term in T3 */

			/**
			 * accumulate the T3 term in H
			 */

			m_add(H, T3, T4); /* addition completed, but assign the T4 back to H */
			m_copy(T4, H);
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

			m_mlt(G_t[counttemp], HOLD[counttemp], T4);
			m_mlt(T4, G[counttemp], T3); /* we got the term in T3 */

			/**
			 * accumulate the T3 term in H
			 */

			m_add(H, T3, T4); /* addition completed, but assign the T4 back to H */
			m_copy(T4, H);
		}

		/**
		 * Now the R TERM part...
		 */
		for(counttemp = 0; counttemp < switcher[2]; counttemp++)
		{
			/**
			 * calculate the (E_t-x %*% E_t-x') term
			 */

			m_set_val(T1, 0, 0, series_1[count - 1 - counttemp]);
			m_set_val(T1, 1, 0, series_2[count - 1 - counttemp]);
			m_transp(T1, T2);
			m_mlt(T1, T2, T3); /* we got the term in T3*/

			/**
			 * calculate the (R_x' %*% E_t-x %*% E_t-x' %*% R_x) term
			 */

			m_mlt(R_t[counttemp], T3, T4);
			m_mlt(T4, R[counttemp], T3); /* we got the term in T3 */


			/**
			 * Critical condition
			 */

			w = params[3 + 4*switcher[0] + 4*switcher[1] + 4*switcher[2]];

			/** if(w*series_1[count - 1 - counttemp] + (1 - w)*series_2[count - 1- counttemp] > 0)
			{
				S = 1.0;
			}
			else
			{
				S = 0.0;
			} */
			/*
                            NEW DEFINITION OF S --- HARALD, 2004-06-18
                        */
                        S = cos(0.7853982 + w) * series_1[count - 1 - counttemp] + sin(0.7853982 + w) * series_2[count - 1 - counttemp];
                        S = S/sqrt((series_1[count - 1 - counttemp])*(series_1[count - 1 - counttemp]) + (series_2[count - 1 - counttemp])*(series_2[count - 1 - counttemp]))+1;
                        S = 1 - 0.5*S;

			sm_mlt(S, T3, T5);

			/**
			 * accumulate the T3 term in H
			 */

			m_add(H, T5, T4); /* addition completed, but assign the T4 back to H */
			m_copy(T4, H);
		}

		/**
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
		 * log[f(x)] = - log[(2 * pi)^(n/2) * (1 / sqrt(|E|))] - (x - m)' * E^(-1) * (x - m) / 2
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
		detcomp = H->me[0][0] * H->me[1][1] - H->me[0][1] * H->me[1][0];

		/**
		 * DEBUG
		 * printf("det(H)² = %f\n", detcomp);
		 */

		if( detcomp == 0 )
		{
			/**
			 * H is not invertable
			 */
			Rprintf("H IS SINGULAR!...\n");

			H->me[0][0] *= 1.01;
			detcomp = H->me[0][0] * H->me[1][1] - H->me[0][1] * H->me[1][0];
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

		m_set_val(T3, 0, 0, H->me[1][1] / detcomp);
		m_set_val(T3, 0, 1, -1 * H->me[0][1] / detcomp);
		m_set_val(T3, 1, 0, -1 * H->me[1][0] / detcomp);
		m_set_val(T3, 1, 1, H->me[0][0] / detcomp);


		m_set_val(T1, 0, 0, series_1[count]); /* initializing (x - m) term (since m = 0; x)*/
		m_set_val(T1, 1, 0, series_2[count]); /* initializing (x - m) term (since m = 0; x)*/

		m_transp(T1, T8);

		m_mlt(T3, T1, T7);
		m_mlt(T8, T7, T0);

		temp = T0->me[0][0];

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
	 printf("\nreturnval : %f\n", *returnval);
	 */

}
