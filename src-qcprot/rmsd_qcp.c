#include "rmsd_qcp.h"

void qcp_(double *weight, double *coords1, double *coords2, int *len, int *mode, double *rmsd, double *rotmat, double *travec, int *ier)
{  
 double wsum_a, wsum_b, sum_a[3], sum_b[3];
 
 int i, j;
 
 double rmsdx;
 
 double **frag_a, **frag_b;
 
 frag_a = MatInit(3, *len);
 frag_b = MatInit(3, *len);
 
 sum_a[0] = 0.0;
 sum_a[1] = 0.0;
 sum_a[2] = 0.0;
 
 sum_b[0] = 0.0;
 sum_b[1] = 0.0;
 sum_b[2] = 0.0;
 
 wsum_a = 0.0;
 wsum_b = 0.0;
 
 for ( i = 0; i < *len; i++ )
 {
  for ( j = 0; j < 3; j++ )
  {
   sum_a[j] += weight[i] * (*(coords1+i*3+j));
   sum_b[j] += weight[i] * (*(coords2+i*3+j));
  }
  
  wsum_a += weight[i];
  wsum_b += weight[i];
 }
 
 sum_a[0] /= wsum_a;
 sum_a[1] /= wsum_a;
 sum_a[2] /= wsum_a;
 
 sum_b[0] /= wsum_b;
 sum_b[1] /= wsum_b;
 sum_b[2] /= wsum_b;
 
 for ( i = 0; i < *len; i++ )
  for ( j = 0; j < 3; j++ )
  {
   frag_a[j][i] = (*(coords1+i*3+j)) - sum_a[j];
   
   frag_b[j][i] = (*(coords2+i*3+j)) - sum_b[j];
  }
 
 rmsdx = CalcRMSDRotationalMatrix((double **) frag_a, (double **) frag_b, *len, rotmat, NULL);
 
 *rmsd = rmsdx*rmsdx*(*len);
 
 *(travec+0) = ((sum_b[0] - (rotmat[0] * sum_a[0])) - (rotmat[3] * sum_a[1])) - (rotmat[6] * sum_a[2]);
 *(travec+1) = ((sum_b[1] - (rotmat[1] * sum_a[0])) - (rotmat[4] * sum_a[1])) - (rotmat[7] * sum_a[2]);
 *(travec+2) = ((sum_b[2] - (rotmat[2] * sum_a[0])) - (rotmat[5] * sum_a[1])) - (rotmat[8] * sum_a[2]);
 
 MatDestroy(&frag_a);
 MatDestroy(&frag_b);
}

double **MatInit(const int rows, const int cols)
{
    int             i;
    double        **matrix = NULL;
    double         *matspace = NULL;

    matspace = (double *) calloc((rows * cols), sizeof(double));
    if (matspace == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate matrix space in MatInit(): (%d x %d)\n", rows, cols);
        exit(EXIT_FAILURE);
    }

    /* allocate room for the pointers to the rows */
    matrix = (double **) malloc(rows * sizeof(double *));
    if (matrix == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate room for row pointers in MatInit(): (%d)\n", rows);
        exit(EXIT_FAILURE);
    }

    /*  now 'point' the pointers */
    for (i = 0; i < rows; i++)
        matrix[i] = matspace + (i * cols);

    return(matrix);
}

void MatDestroy(double ***matrix_ptr)
{
    double **matrix = *matrix_ptr;

    if (matrix != NULL)
    {
        if (matrix[0] != NULL)
        {
            free(matrix[0]);
            matrix[0] = NULL;
        }

        free(matrix);
        *matrix_ptr = NULL;
    }
}
