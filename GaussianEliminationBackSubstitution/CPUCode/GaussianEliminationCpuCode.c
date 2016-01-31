#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <MaxSLiCInterface.h>
#include "Maxfiles.h"
#define ABS(x) (x < 0 ? -(x) : (x))
#define EPS 0.00001
#define TRUE  1
#define FALSE 0

void WriteSolution(float *a,int n,float *x)
{
   int j,k;

   for (j=0;j<n;j++) {
      for (k=0;k<n+1;k++) {
         printf("%10.3f ",a[k*n+j]);
      }
      printf(" | %10.3f\n",x[j]);
   }
}

int GSolve(float *a,int n,float *x, int streamSize, int vectorSize)
{
   int i,j,k,maxrow;
   float tmp;
   int redBr1;
   int redBr2;
   size_t sizeBytes = streamSize * vectorSize * sizeof(float);
   float *dataIn = malloc(sizeBytes);
   float *dataOut = malloc(sizeBytes);


   for (i=0;i<n;i++) {

	  for(int i = 0; i < (streamSize * vectorSize); i++) {
	   		dataOut[i] = 0;
	  }

      /* Find the row with the largest first value */
      maxrow = i;
      for (j=i+1;j<n;j++) {
         if (ABS(a[i*n+j]) > ABS(a[i*n+maxrow]))
            maxrow = j;
      }

      redBr1=0;
      redBr2=1;
      for (k=i;k<n+1;k++) {
    	  dataIn[redBr1]=a[k*n+i];
    	  dataIn[redBr2]=a[k*n+maxrow];
    	  redBr1=redBr1+2;
    	  redBr2=redBr2+2;
      }

      /* Swap the maxrow and ith row */
      GaussianElimination(streamSize, dataIn, dataOut);

      redBr1=0;
      redBr2=1;
      for (k=i;k<n+1;k++) {
    	  a[k*n+i]=dataOut[redBr1];
    	  a[k*n+maxrow]=dataOut[redBr2];
    	  redBr1=redBr1+2;
    	  redBr2=redBr2+2;
      }

      /* Singular matrix? */
      if (ABS(a[i*n+i]) < EPS)
         return(FALSE);

      /* Eliminate the ith element of the jth row */
      for (j=i+1;j<n;j++) {
         for (k=n;k>=i;k--) {
            a[k*n+j] -= a[k*n+i] * a[i*n+j] / a[i*n+i];
         }
      }
   }

   /* Do the back substitution */
   for (j=n-1;j>=0;j--) {
      tmp = 0;
      for (k=j+1;k<n;k++)
         tmp += a[k*n+j] * x[k];
      x[j] = (a[n*n+j] - tmp) / a[j*n+j];
   }

   return(TRUE);
}

int main()
{
	const int streamSize = 4;
	const int vectorSize = GaussianElimination_vectorSize;
	int n = 3;
	float x[3] = {0.0,0.0,0.0};
	float *a;

	a = (float *)malloc(n*(n+1)*sizeof(float));

	a[0] = 1;   a[3] = 1;   a[6] = 1;  a[9] = 0;
	a[1] = 2;   a[4] = 1;   a[7] = 1;  a[10] = 1;
	a[2] = 1;   a[5] = 2;   a[8] = 1;  a[11] = 15;

	printf("Start maxtrix:\n");
	WriteSolution(a,n,x);
	GSolve(a,n,x,streamSize,vectorSize);
	printf("End matrix:\n");
	WriteSolution(a,n,x);

	return 0;
}

