#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc,char** argv) {
  char c;
  char *str="abcdefghijklmnopqrstuvwxyz";
  bool b=true;
  short s=1;
  int i=2;
  long l=3L;
  float f=4.;
  double d=5.;
  FILE *fout;
  FILE *fin;
  FILE *uout;
  FILE *uin;
  char cc;
  char *strstr=0;
  bool bb;
  short ss;
  int ii;
  long ll;
  float ff;
  double dd;
  int n,n2;
  size_t r;
  char *str2=(char*) malloc(128*sizeof(char));
  char *str3=(char*) malloc(128*sizeof(char));
  char c2;

/*formatted to stdout*/
  printf("in cio\n");
  for (c='A';c<='Z';c++) putchar(c);
  putchar('\n');
  puts(str);
  printf("bool = %s\n",b ? "true" : "false");
  printf("short = %hd\n",s);
  printf("int = %d\n",i);
  printf("long = %ld\n",l);
  printf("float = %2.0f\n",f);
  printf("double = %2.0f\n",d);

/*formatted to FILE*/
  printf("\nformatted write to file\n");
  fout=fopen("cio_formatted_output","w");
  for (c='A';c<='Z';c++) fputc(c,fout);
  fputc('\n',fout);
  n=strlen(str);
  fprintf(fout,"strlen = %d\n",n);
  fputs(str,fout);
  fputc('\n',fout);
  fprintf(fout,"bool = %s\n",b ? "true" : "false");
/*printf("after write bool, ftell = %ld\n",ftell(fout)); */
  fprintf(fout,"short = %hd\n",s);
/*printf("after write short, ftell = %ld\n",ftell(fout)); */
  fprintf(fout,"int = %d\n",i);
  fprintf(fout,"long = %ld\n",l);
  fprintf(fout,"float = %2.0f\n",f);
  fprintf(fout,"double = %2.0f\n",d);
  fclose(fout);

/*formatted from FILE*/
  printf("\nformatted read from file\n");
  fin=fopen("cio_formatted_output","r");
  while (true) {
    cc=fgetc(fin);
    putchar(cc);
    if (cc=='\n') break;
  }
  r=fscanf(fin,"%s %c %d",str2,&c2,&n2);
  cc=fgetc(fin);
  printf("strlen = %d\n",n2);
  strstr=(char*) malloc((n2+1)*sizeof(char));
  fgets(strstr,n2,fin);
  cc=fgetc(fin);
  puts(strstr);
/*putchar('\n'); */
  r=fscanf(fin,"%s %c %s",str2,&c2,str3);
  cc=fgetc(fin);
/*printf("after read bool = %s, ftell = %ld\n",str3,ftell(fin)); */
  bb=(strcmp(str3,"true")==0 ? true : false);
  printf("bool = %s\n",bb ? "true" : "false");
  fscanf(fin,"%s %c %hd\n",str2,&c2,&ss);
  cc=fgetc(fin);
  printf("short = %hd\n",ss);
  fscanf(fin,"%s %c %d\n",str2,&c2,&ii);
  cc=fgetc(fin);
  printf("int = %d\n",ii);
  fscanf(fin,"%s %c %ld\n",str2,&c2,&ll);
  cc=fgetc(fin);
  printf("long = %ld\n",ll);
  fscanf(fin,"%s %c %f\n",str2,&c2,&ff);
  cc=fgetc(fin);
  printf("float = %2.0f\n",ff);
  fscanf(fin,"%s %c %lf\n",str2,&c2,&dd);
  cc=fgetc(fin);
  printf("double = %2.0f\n",dd);
  fclose(fin);
  free(strstr); strstr=0;

/*unformatted to FILE*/
  printf("\nunformatted write to file\n");
  uout=fopen("cio_unformatted_output","w");
  for (c='A';c<='Z';c++) fwrite(&c,sizeof(char),1,uout);
  c='\n';
  fwrite(&c,sizeof(char),1,uout);
  n=strlen(str);
  fwrite(&n,sizeof(int),1,uout);
  fwrite(str,sizeof(char),strlen(str)+1,uout);
  fwrite(&b,sizeof(bool),1,uout);
  fwrite(&s,sizeof(short),1,uout);
  fwrite(&i,sizeof(int),1,uout);
  fwrite(&l,sizeof(long),1,uout);
  fwrite(&f,sizeof(float),1,uout);
  fwrite(&d,sizeof(double),1,uout);
  fclose(uout);

/*unformatted from FILE*/
  printf("\nunformatted read from file\n");
  uin=fopen("cio_unformatted_output","r");
  while (true) {
    r=fread(&cc,sizeof(char),1,uin);
    putchar(cc);
    if (cc=='\n') break;
  }
  r=fread(&n,sizeof(int),1,uin);
  strstr=(char*) malloc((n+1)*sizeof(char));
  fread(strstr,sizeof(char),n+1,uin);
  puts(strstr);
/*putchar('\n'); */
  fread(&bb,sizeof(bool),1,uin);
  printf("bool = %s\n",bb ? "true" : "false");
  fread(&ss,sizeof(short),1,uin);
  printf("short = %hd\n",ss);
  fread(&ii,sizeof(int),1,uin);
  printf("int = %d\n",ii);
  fread(&ll,sizeof(long),1,uin);
  printf("long = %ld\n",ll);
  fread(&ff,sizeof(float),1,uin);
  printf("float = %2.0f\n",ff);
  fread(&dd,sizeof(double),1,uin);
  printf("double = %2.0f\n",dd);
  fclose(uin);
  free(strstr); strstr=0;

  free(str2); str2=0;
  free(str3); str3=0;
  return EXIT_SUCCESS;
}
