#include <stdio.h> /* for printf */
#include <stdlib.h> /* for malloc,drand48 */
#include <sys/times.h> /* for function times and struct tms */
#include <unistd.h> /* for function sysconf */
int main(int argc,char **argv) {
/* argc is number of words in the command line, and
 * argv is an array of character strings of the words on the command line
 */
  int i,j,k,m,n,incx,incy,incAcolumn;
  double elapsed,xj,x[512],y[512]; /* 64-bit floating point */
  double A[512][512]; /* A[column][row] */
  double *A_dynamic=0; /* pointer to double for dynamic memory allocation */
  struct tms usage; /* in sys/time.h */
  clock_t start; /* clock_t in sys/time.h */
/* other data types:
 *   float complex c;
 *   double complex z;
 */
  double alpha=1.,beta=0.;
  char trans='N';

  m=n=512;
/* dynamic memory allocation: */
  A_dynamic = (double*) malloc(m*n*sizeof(double));
  incx=incy=incAcolumn=1;
/*fill array A with random numbers */
  for (j=0;j<n;j++) {
    for (i=0;i<m;i++) {
      A[j][i]=drand48(); // in stdlib.h
      A_dynamic[i+j*m]=drand48();
    }
    x[j]=drand48();
  }
  for (i=0;i<m;i++) y[i]=0.;

/*times(&usage) will return
 * usage.tms_utime = User CPU time
 * usage.tms_stime = System CPU time
 * usage.tms_cutime = User CPU time of dead children
 * usage.tms_cstime = System CPU time of dead children
 * all are measured in CLK_TCKths of a second
 */
  times(&usage);
  start=usage.tms_utime;
  for (k=0;k<4096;k++) {
    for (j=0;j<n;j++) {
      xj=x[j];
      double *Aj=A[j]; /* save array access time */
      for (i=0;i<m;i++) {
        y[i]+=xj*Aj[i];
        /* alternatively,
         * y[i] += x[j] * A[j][i];
         * y[i] += A_dynamic[i+j*m] * x[j];
         */
      }
    }
  }
  times(&usage);
  elapsed=(double)(usage.tms_utime-start)
    /((double)(sysconf(_SC_CLK_TCK))*4096.);
  printf("time for double loop j over i = %.3g\n",elapsed);

  times(&usage);
  start=usage.tms_utime;
  for (k=0;k<4096;k++) {
    for (j=0;j<n;j++) {
      daxpy_(&m,&x[j],&A[j][0],&incAcolumn,y,&incy);
/*    note that we did not need to copy the j'th column of A to a vector */
/*    alternatively, we could use A_dynamic:
 *    daxpy_(&m,&x[j],&A_dynamic[j*m],&incAcolumn,y,&incy);
 */
    }
  }
  times(&usage);
  elapsed=(double)(usage.tms_utime-start)
    /((double)(sysconf(_SC_CLK_TCK))*4096.);
  printf("time for loop j over daxpy = %.3g\n",elapsed);

  times(&usage);
  start=usage.tms_utime;
  for (k=0;k<4096;k++) {
    for (i=0;i<m;i++) {
      y[i]=ddot_(&n,&A[0][i],&m,x,&incx);
/*    note that we did not need to copy the i'th row of A to a vector */
/*    alternatively, we could use A_dynamic:
 *    y[i]=ddot_(&n,&A_dynamic[i],&m,x,&incx);
 */
    }
  }
  times(&usage);
  elapsed=(double)(usage.tms_utime-start)
    /((double)(sysconf(_SC_CLK_TCK))*4096.);
  printf("time for loop i over ddot = %.3g\n",elapsed);

  times(&usage);
  start=usage.tms_utime;
  for (k=0;k<4096;k++) {
    dgemv_(&trans,&m,&n,&alpha,A,&m,x,&incx,&beta,y,&incy);
  }
/* alternatively, we could use A_dynamic:
 * dgemv_(&trans,&m,&n,&alpha,A_dynamic,&m,x,&incx,&beta,y,&incy);
 */
  times(&usage);
  elapsed=(double)(usage.tms_utime-start)
    /((double)(sysconf(_SC_CLK_TCK))*4096.);
  printf("time for dgemv = %.3g\n",elapsed);

  free(A_dynamic); A_dynamic=0; /* deallocation */
}
