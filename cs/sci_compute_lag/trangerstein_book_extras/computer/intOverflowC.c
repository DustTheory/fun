#include <stdio.h> /* for FILE, scanf, printf, fprintf */
#include <stdlib.h> /* for EXIT_SUCCESS */
#include <values.h> /* for INT_MAX */
int main(int argc,char **argv) {
  int i;
  short s;
  long l;
  long long ll;
  unsigned int uold,unew;
//
  printf("sizeof(int) = %d\n",sizeof(int));
  i=-1;
  printf("(int) -1 = %d %#x\n",i,i);
  i=INT_MAX;
  printf("INT_MAX = %d %#x\n",i,i);
  i=INT_MIN;
  printf("INT_MIN = %d %#x\n\n",i,i);

  printf("sizeof(short) = %d\n",sizeof(short));
  s=-1;
  printf("(short) -1 = %d %#x\n",s,s);
  s=SHRT_MAX;
  printf("SHRT_MAX = %d %#x\n",s,s);
  s=SHRT_MIN;
  printf("SHRT_MIN = %d %#x\n\n",s,s);

  printf("sizeof(long) = %d\n",sizeof(long));
  l=-1;
  printf("(long) -1 = %d %#x\n",l,l);
  l=LONG_MAX;
  printf("LONG_MAX = %d %#x\n",l,l);
  l=LONG_MIN;
  printf("LONG_MIN = %d %#x\n\n",l,l);

  printf("sizeof(long long) = %d\n\n",sizeof(long long));
  ll=-1;
/*
  printf("(long long) -1 = %d %#x\n\n",ll,ll);
  ll=LONG_MAX;
  printf("LONG_MAX = %d %#x\n",ll,ll);
  ll=LONG_MIN;
  printf("LONG_MIN = %d %#x\n",ll,ll);
*/
//
  i=1;
  printf("compute j_n = 2*j_{n-1} + 1, j_0 = 1 while j_n > 0\n");
  while (i>0) {
    printf("%d %#x\n",i,i);
    i=2*i+1;
  }
  printf("%d %#x\n\n",i,i);

  printf("compute j_n = 2*j_{n-1} - 1, j_0 = -1 while j_n < 0\n");
  i=-1;
  while (i<0) {
    printf("%d %#x\n",i,i);
    i=2*i-1;
  }
  printf("%d %#x\n\n",i,i);

  uold=0,unew=1;
  printf("compute u_n = 2*u_{n-1} + 1, u_0 = 0 while u_n > u_{n-1}\n");
  while (unew>uold) {
    uold=unew;
    printf("%d %#x\n",unew,unew);
    unew=2*uold+1;
  }
  printf("%d %#x\n\n",unew,unew);
  return EXIT_SUCCESS;
}
