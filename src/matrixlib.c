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

/**
 * functions to be written:
 * !m_get(int rows, int cols);
 * !m_set_val(MAT * matrix_to_be_filled, int row_number, int col_number, double value);
 * m_transp(MAT * matrix_to_be_transposed);
 * m_mlt(MAT * first, MAT * second, MAT * result);
 * !m_copy(MAT * source, MAT * destination);
 * !m_add(MAT * first, MAT * second, MAT * result);
 */

int err_exit(char * err_message, char * func_name)
{
	Rprintf("\nERROR in matrixlib:\n\t %s", err_message);
	Rprintf("\n(in function %s)", func_name);
	Rprintf("\nExiting program\n");
        error("Exiting...");
}

void mem_copy(int * source, int * destination, unsigned int length)
{
	int i;
	if(source < destination)
	{
		for(i = 0; i < length; i++)
		{
			*(destination++) = *(source++);
		}
	}
	else
	{
		source += length;
		destination += length;

		for(i = 0; i < length; i++)
		{
			*(--destination) = *(--source);
		}
	}
	return;
}

void m_free(MAT * mat)
{
	int i;
	if(mat==(MAT *)NULL || (int)(mat->rows) < 0 || (int)(mat->cols) < 0)
	{
		return;
	}
	for(i = 0; i < mat->rows; i++)
	{
		free((char *)(mat->me[i]));
	}
	free((char *)(mat->me));
	free((char *)mat);
}

void m_transp(MAT * in, MAT * out)
{
	int	i, j;
	int	same;
	double tmp;

	if(in == NULL)
	{
		err_exit(ERR_SRCNULL, "m_transp");
	}

	if(out == NULL)
	{
		err_exit(ERR_DESTNULL, "m_transp");
	}

	same = ( in == out );

	if(same && in->rows != in->cols)
	{
		err_exit(ERR_TRANSPDIFFDIMS, "m_transp");
	}

	if(out->rows != in->cols || out->cols != in->rows)
	{
		err_exit(ERR_TRANSPDIFFDIMS, "m_transp");
	}

	if(!same)
	{
		for(i = 0; i < in->rows; i++)
		{
			for(j = 0; j < in->cols; j++)
			{
				out->me[j][i] = in->me[i][j];
			}
		}
	}
	else
	{
		for(i = 1; i < in->rows; i++)
		{
			for (j = 0; j < i; j++)
			{	tmp = in->me[i][j];
				in->me[i][j] = in->me[j][i];
				in->me[j][i] = tmp;
			}
		}
	}

	return;
}

void m_foutput(FILE * fp, MAT * a)
{
	unsigned int i, j, tmp;

	if(a == NULL) /* check whether the matrix exist or not*/
	{
		fprintf(fp, "MATRIX: NULL\n");
		return;
	}

	fprintf(fp, "Matrix: %d by %d\n", a->rows, a->cols);

	if(a->me == NULL) /* check the content array */
	{
		fprintf(fp, "NULL\n");
	}

	for(i = 0; i < a->rows; i++) /* for each row */
	{
		fprintf(fp, "row %u: ", i);
		for(j = 0, tmp = 2; j < a->cols; j++, tmp++) /* for each column */
		{
			fprintf(fp, "%f ", a->me[i][j]);
			if(!(tmp % 5))
			{
				putc('\n', fp);
			}
		}
		if((tmp % 5) != 1)
		{
			putc('\n', fp);
		}
	}
}

MAT * m_get(int m, int n)
{
	MAT * matrix;	/* matrix to be returned */
	int	i;			/* temporary variable */

	/* check the dimensions, both should be equal or greater then zero */
	if(m < 0 || n < 0)
	{
		err_exit(ERR_DIMNEG, "m_get");
	}

	/* atempt to allocate memory for the matrix */
	matrix = NEW_MM(MAT);
	if(matrix == (MAT *) NULL)
	{
		err_exit(ERR_MEM, "m_get");
	}

	matrix->rows = m;
	matrix->cols = n;

	matrix->me = NEW_A(m, double*);
	if(matrix->me == NULL)
	{
		free(matrix);
		err_exit(ERR_MEM, "m_get");
	}

	for(i = 0; i < m; i++)
	{
		matrix->me[i] = NEW_A(n, double);
		if(matrix->me[i] == NULL)
		{
			free(matrix->me);
			free(matrix);
			err_exit(ERR_MEM, "m_get");
		}
	}
	return matrix;
}

void m_copy(MAT * source, MAT * destination)
{
	unsigned int i, j;

	if(source == NULL) /* source cannot be NULL */
	{
		err_exit(ERR_SRCNULL, "m_copy");
	}

	if(destination == NULL) /* destination should be initialized */
	{
		err_exit(ERR_DESTNULL, "m_copy");
	}

	if(
		(source->cols != destination->cols) ||
		(source->rows != destination->cols)) /* not same-dimensioned matrices */
	{
		err_exit(ERR_DIFFDIMS, "m_copy");
	}

	if(source == destination) /* nothing to be done */
	{
		return;
	}

	for( i = 0; i < source->rows; i++)
	{
		/**
		 * TODO
		 * THERE IS A SIGNIFICANT BUG IN MEM_COPY FUNCTION, OR JUST HERE IN THE BELOW LINE
		 * mem_copy((int *)&(source->me[i][0]), (int *)&(destination->me[i][0]), source->cols * sizeof(double));
		 */
		for(j = 0; j < source->cols; j++)
		{
			destination->me[i][j] = source->me[i][j];
		}
	}

	return;
}

void m_add(MAT * mat1, MAT * mat2, MAT * result)
{
	unsigned int i, j;

	if(mat1 == NULL || mat2 == NULL) /* operands cannot be NULL */
	{
		err_exit(ERR_OPSNULL, "m_add");
	}

	if(result == NULL) /* result should be initialized */
	{
		err_exit(ERR_RESNULL, "m_add");
	}

	if(
		(mat1->cols != mat2->cols) ||
		(mat1->cols != result->cols) ||
		(mat1->rows != mat2->cols) ||
		(mat1->rows != result->rows)) /* not same-dimensioned matrices */
	{
		err_exit(ERR_DIFFDIMS, "m_add");
	}

	/**
	 * TODO
	 * Think about whether the values of the operands would be overwritten
	 * or not if we allow that the result is the same matrix of one of the
	 * operands. It seems ok, since the mathematical operations like
	 * "a = a + b" are ok.
	 */

	for( i = 0; i < mat1->rows; i++)
	{
		for(j = 0; j < mat1->cols; j++)
		{
			result->me[i][j] = mat1->me[i][j] + mat2->me[i][j];
		}
	}

	return;
}

void m_mlt(MAT * A, MAT * B, MAT * OUT)
{
	unsigned int i, j, k, m, n, p;
	double **A_v, **B_v, sum;

	if(A == NULL || B == NULL)
	{
		err_exit(ERR_OPSNULL, "m_mlt");
	}

	if(OUT == NULL)
	{
		err_exit(ERR_RESNULL, "m_mlt");
	}

 	if(A->cols != B->rows)
	{
		/**
		 * DEBUG
		 * printf("A->cols = %d\n", A->cols);
		 * printf("B->rows = %d\n", B->rows);
		 */
		err_exit(ERR_OPDIFFDIMS, "m_mlt");
	}

	if(A == OUT || B == OUT)
	{
		err_exit(ERR_OPRESSAME, "m_mlt");
	}

	m = A->rows;
	n = A->cols;
	p = B->cols;

	A_v = A->me;
	B_v = B->me;

	if(OUT->rows != A->rows || OUT->cols != B->cols)
	{
		/**
		 * DEBUG
		 * printf("A->rows = %d, B->cols = %d,\n OUT->rows = %d,\n OUT->cols = %d\n", A->rows, B->cols, OUT->rows, OUT->cols);
		 */
		err_exit(ERR_OPDIFFDIMS, "m_mlt");
	}

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < p; j++)
		{
			sum = 0.0;
			for(k = 0; k < n; k++)
			{
				sum += A_v[i][k] * B_v[k][j];
			}
			OUT->me[i][j] = sum;
		}
	}

	return;
}

void m_scalar_divide(MAT * arg, MAT * result, double divider)
{
	int i, j;
	for(i = 0; i < arg->rows; i++)
	{
		for(j = 0; j < arg->cols; j++)
		{
			result->me[i][j] = arg->me[i][j] / divider;
		}
	}
}

/**
 * DETERMINANT
 */
/**
 * TODO
 * WRITE THIS PART ON YOUR OWN
 */
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return det;
}

double m_det(MAT *mat)
{
	return Determinant(mat->me, mat->rows);
}

/**
 * INVERSE
 */
/**
 * TODO
 * WRITE THIS PART ON YOUR OWN
 */
/*
   Transpose of a square matrix, do it in place
*/
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}
/*
   Find the cofactor matrix of a square matrix
*/
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}

void m_inverse(MAT *feed, MAT *result)
{
	MAT * tmp;
	tmp = m_get(feed->rows, feed->cols);
	CoFactor(feed->me, feed->cols, tmp->me);
	Transpose(tmp->me, tmp->cols);
	m_scalar_divide(tmp, result, m_det(feed));
	m_free(tmp);
}

void sm_mlt(double scalar, MAT * B, MAT * OUT)
{
	unsigned int j, k, m, p;
	double **B_v/*, mul */;

	if(B == NULL)
	{
		err_exit(ERR_OPSNULL, "sm_mlt");
	}

	if(OUT == NULL)
	{
		err_exit(ERR_RESNULL, "sm_mlt");
	}

	m = B->rows;
	p = B->cols;

	B_v = B->me;

	if(OUT->rows != B->rows || OUT->cols != B->cols)
	{
		/**
		 * DEBUG
		 * printf("A->rows = %d, B->cols = %d,\n OUT->rows = %d,\n OUT->cols = %d\n", A->rows, B->cols, OUT->rows, OUT->cols);
		 */
		err_exit(ERR_OPDIFFDIMS, "sm_mlt");
	}

	for(j = 0; j < m; j++)
	{
		/* mul = 0.0; */
		for(k = 0; k < p; k++)
		{
			OUT->me[j][k] = scalar * B_v[j][k];
		}
	}

	return;
}
