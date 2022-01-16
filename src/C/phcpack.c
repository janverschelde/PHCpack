/* file phcpack.c contains the definitions of the functions in phcpack.h */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "phcpack.h"

#define v 0 /* verbose flag */

/* most BASIC operations in PHCpack : */

int version_string ( int *n, char *s )
{
   int fail,i;
   int b[40];
   double *c;

   fail = _ada_use_c2phc(999,n,b,c,0);

   for(i=0; i<(*n); i++) s[i] = (char) b[i];
   s[*n] = '\0';

   return fail;
}

int set_seed ( int seed )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(998,&seed,b,c,0);

   return fail;
}

int get_seed ( int *seed )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(997,seed,b,c,0);

   return fail;
}

int solve_standard_system
 ( int *root_count, int silent, int *nrcs, char *rocos, int nbtasks,
   int mvfocus, int vrb )
{
   int fail,i;
   int a[2];
   int b[1024];
   double *c;

   b[0] = silent;
   b[1] = nbtasks;
   b[2] = mvfocus;

   fail = _ada_use_c2phc(77,a,b,c,vrb);

   *root_count = a[0];

   if(silent == 0)
   {
      *nrcs = a[1];
      for(i=0; i<(*nrcs); i++) rocos[i] = (char) b[i];
      rocos[*nrcs] = '\0';
   }
   return fail;
}

int solve_dobldobl_system
 ( int *root_count, int silent, int *nrcs, char *rocos, int nbtasks, int vrb )
{
   int fail,i;
   int a[2];
   int b[1024];
   double *c;

   b[0] = silent;
   b[1] = nbtasks;

   fail = _ada_use_c2phc(700,a,b,c,vrb);

   *root_count = a[0];

   if(silent == 0)
   {
      *nrcs = a[1];
      for(i=0; i<(*nrcs); i++) rocos[i] = (char) b[i];
      rocos[*nrcs] = '\0';
   }
   return fail;
}

int solve_quaddobl_system
 ( int *root_count, int silent, int *nrcs, char *rocos, int nbtasks, int vrb )
{
   int fail,i;
   int a[2];
   int b[1024];
   double *c;

   b[0] = silent;
   b[1] = nbtasks;

   fail = _ada_use_c2phc(702,a,b,c,vrb);

   *root_count = a[0];

   if(silent == 0)
   {
      *nrcs = a[1];
      for(i=0; i<(*nrcs); i++) rocos[i] = (char) b[i];
      rocos[*nrcs] = '\0';
   }
   return fail;
}

int solve_standard_Laurent_system
 ( int *root_count, int silent, int *nrcs, char *rocos, int nbtasks,
   int mvfocus, int vrb )
{
   int fail,i;
   int a[2];
   int b[1024];
   double *c;

   b[0] = silent;
   b[1] = nbtasks;
   b[2] = mvfocus;

   fail = _ada_use_c2phc(75,a,b,c,vrb);

   *root_count = a[0];

   if(silent == 0)
   {
      *nrcs = a[1];
      for(i=0; i<(*nrcs); i++) rocos[i] = (char) b[i];
      rocos[*nrcs] = '\0';
   }
   return fail;
}

int solve_dobldobl_Laurent_system
 ( int *root_count, int silent, int *nrcs, char *rocos, int nbtasks, int vrb )
{
   int fail,i;
   int a[2];
   int b[1024];
   double *c;

   b[0] = silent;
   b[1] = nbtasks;

   fail = _ada_use_c2phc(701,a,b,c,vrb);

   *root_count = a[0];

   if(silent == 0)
   {
      *nrcs = a[1];
      for(i=0; i<(*nrcs); i++) rocos[i] = (char) b[i];
      rocos[*nrcs] = '\0';
   }
   return fail;
}

int solve_quaddobl_Laurent_system
 ( int *root_count, int silent, int *nrcs, char *rocos, int nbtasks, int vrb )
{
   int fail,i;
   int a[2];
   int b[1024];
   double *c;

   b[0] = silent;
   b[1] = nbtasks;

   fail = _ada_use_c2phc(703,a,b,c,vrb);

   *root_count = a[0];

   if(silent == 0)
   {
      *nrcs = a[1];
      for(i=0; i<(*nrcs); i++) rocos[i] = (char) b[i];
      rocos[*nrcs] = '\0';
   }
   return fail;
}

int set_gamma_constant
 ( double regamma, double imgamma, int precision, int vrb )
{
   int fail = 1;
   int *b;
   double c[2];
  
   c[0] = regamma;
   c[1] = imgamma;

   fail = _ada_use_c2phc(996,&precision,b,c,vrb);

   return fail;
}

int get_gamma_constant
 ( double *regamma, double *imgamma, int precision, int vrb )
{
   int fail = 1;
   int *b;
   double c[2];
  
   fail = _ada_use_c2phc(995,&precision,b,c,vrb);

   *regamma = c[0];
   *imgamma = c[1];

   return fail;
}

int mixed_volume ( int *mv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(78,mv,b,c,0);

   return fail;
}

int stable_mixed_volume ( int *mv, int *smv )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(79,mv,smv,c,0);

   return fail;
}

int mixed_volume_by_demics ( int *mv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(843,mv,b,c,0);

   return fail;
}

int stable_mixed_volume_by_demics ( int *mv, int *smv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(844,mv,smv,c,0);

   return fail;
}

int standard_deflate
 ( int maxitr, int maxdef, double tolerr, double tolres, double tolrnk )
{
   int fail;
   double tols[3];

   tols[0] = tolerr;
   tols[1] = tolres;
   tols[2] = tolrnk;

   fail = _ada_use_c2phc(196,&maxitr,&maxdef,tols,0);

   return fail;
}

int dobldobl_deflate
 ( int maxitr, int maxdef, double tolerr, double tolres, double tolrnk )
{
   int fail;
   double tols[3];

   tols[0] = tolerr;
   tols[1] = tolres;
   tols[2] = tolrnk;

   fail = _ada_use_c2phc(249,&maxitr,&maxdef,tols,0);

   return fail;
}

int quaddobl_deflate
 ( int maxitr, int maxdef, double tolerr, double tolres, double tolrnk )
{
   int fail;
   double tols[3];

   tols[0] = tolerr;
   tols[1] = tolres;
   tols[2] = tolrnk;

   fail = _ada_use_c2phc(250,&maxitr,&maxdef,tols,0);

   return fail;
}

int standard_Newton_step ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(199,a,b,c,0);

   return fail;
}

int dobldobl_Newton_step ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(198,a,b,c,0);

   return fail;
}

int quaddobl_Newton_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(197,a,b,c,0);
   return fail;
}

int multprec_Newton_step ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(195,&deci,b,c,0);
   return fail;
}

int standard_Newton_Laurent_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(326,a,b,c,0);
   return fail;
}

int dobldobl_Newton_Laurent_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(327,a,b,c,0);
   return fail;
}

int quaddobl_Newton_Laurent_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(328,a,b,c,0);
   return fail;
}

int multprec_Newton_Laurent_step ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(329,&deci,b,c,0);
   return fail;
}

char *read_equations_from_file
 ( FILE *fp, int nq, int k, int *len, char *accu )
{
   int i,lenterm;
   char term[256];
   char *result;

   fscanf(fp,"%s\n",term);
   lenterm = strlen(term);
   for(i=0; i<lenterm; i++) if(term[i] == ';') k = k+1;

   result = (char*)calloc(*len+lenterm+1,sizeof(char));
   for(i=0; i<*len; i++) result[i] = accu[i];
   result[*len] = ' '; /* add space to separate terms */
   for(i=0; i<lenterm; i++) result[*len+1+i] = term[i];
   *len = *len + lenterm + 1;
   if(k < nq)
      return read_equations_from_file(fp,nq,k,len,result);
   else
   {
      free(accu);
      return result;
   }
}

int scan_number_of_variables ( int nc, char *eqs, int *dim )
{
   int fail = 0;
   int i,b[nc];
   double *c;

   for(i=0; i<nc; i++) b[i] = (int) eqs[i];

   *dim = nc;

   fail = _ada_use_c2phc(439,dim,b,c,0);

   return fail;
}

char *read_polynomials_from_file
 ( int nc, char *name, int *len, int *nq, int *nv, int *fail )
{
   char *result;
   FILE *fp;
   
   fp = fopen(name,"r");
   if(fp == NULL)
   {
      printf("File with name %s could not be opened for reading!\n",name);
      *fail = 1;
   }
   else
   {
      char c;
      char *acc;
      fscanf(fp,"%d",nq);
      /* printf("Number of equations : %d\n",*nq); */
      c = getc(fp);
      if(c == '\n')
         *nv = *nq;
      else
         fscanf(fp,"%d",nv);
      /* printf("Number of variables : %d\n",*nv); */
      acc = (char*)calloc(1,sizeof(char));
      acc[0] = '\0';
      *len = 0;
      result = read_equations_from_file(fp,*nq,0,len,acc);
      *fail = 0;
   }
   fclose(fp);

   return result;
}

int skip_lines ( FILE *fp, int k )
{
   char c;
   int go_on = 1;
   int cnt = k;   /* number of new lines to read */

   do
   {
      go_on = fscanf(fp,"%c",&c);
      if(go_on != 1) break;         /* at end of file */
      if(c == '\n') cnt = cnt - 1;  /* one less line to read */
      go_on = (cnt > 0);
   }
   while(go_on == 1);

   return cnt;
}

char *buffered_line_reader 
 ( FILE *fp, int k, int n, int *len, char *accu )
{
   const int SIZE = 80;
   char *result;
   char buffer[SIZE];
   char c;
   int i;
   int go_on = 1;
   int idx = 0;
   int cnt = k;

   for(i=0; i<SIZE; i++) buffer[i] = ' ';

   do
   {
      int rc = fgetc(fp);
      if(rc == EOF) break;
      c = (char) rc;
      if(c == '\n') cnt = cnt + 1;
      buffer[idx++] = (char) rc;
   }
   while((idx < SIZE) && (cnt < n));

   result = (char*)calloc(*len+idx+1,sizeof(char));
   /* The initialization of result seems critical !!! */
   for(i=0; i<*len+idx+1; i++) result[i] = ' ';
   for(i=0; i<*len; i++) result[i] = accu[i];
   for(i=0; i<idx; i++) result[*len+i] = buffer[i];
   *len = *len + idx;
   result[*len] = '\0';

   if(cnt < n)
      return buffered_line_reader(fp,cnt,n,len,result);
   else
   {
      free(accu);
      return result;
   }
}

char *store_lines ( FILE *fp, int k )
{
   char *result;
   char *acc;
   int len = 0;
 
   acc = (char*)calloc(1,sizeof(char));
   acc[0] = '\0';
   
   result = buffered_line_reader(fp,0,k,&len,acc);

   return result;
}

int read_solution_banner ( FILE *fp, int *len, int *dim )
{
   const int size = 13;
   const char banner[] = "THE SOLUTIONS";
   char buffer[size];
   char c;
   int go_on = 1;
   int cnt = 0;
   int found = 0;

   do
   {
      go_on = fscanf(fp,"%c",&c);
      if(go_on != 1) break;         /* at end of file */
      if(cnt < size)
         buffer[cnt++] = c;         /* just add character to buffer */
      else
      {                             /* shift the characters in buffer */
         int i;
         for(i=0; i<size-1; i++) buffer[i] = buffer[i+1];
         buffer[size-1] = c;
         found = 1;                 /* assume we have a match */
         for(i=0; i<size; i++)
            if(buffer[i] != banner[i])
            { 
               found = 0; break;
            }
         go_on = (found == 0);
      }
   }
   while(go_on == 1);

   if(found == 0)
      printf("Did not find the solution banner!\n");
   else
   {
      go_on = skip_lines(fp,1);     /* skip rest of the banner line */
      if(go_on == 0)
      {
         fscanf(fp,"%d",len);
         fscanf(fp,"%d",dim);
         printf("Length : %d and dimension : %d\n",*len,*dim);
         go_on = skip_lines(fp,2);  /* skip rest of line and next line */
      }
   }

   return go_on;
}

char *read_solution_string ( FILE *fp, int k, int len, int dim )
{
   char *result;
   const int lines = dim + 5;     /* number of lines a solution occupies */
   int fail;

   if(k > 1) fail = skip_lines(fp,(k-1)*lines);

   if(fail == 0)
   {
      fail = skip_lines(fp,1);        /* skip line with solution counter */
      result = store_lines(fp,lines-1);
   }

   return result;
}

char *read_solution_banner_and_string ( FILE *fp, int k, int *len, int *dim )
{
   int fail;
   char *result;

   fail = read_solution_banner(fp,len,dim);
   if(fail == 0) result = read_solution_string(fp,k,*len,*dim);

   return result;
}

int varbprec_Newton_Laurent_step
 ( int dim, int wanted, int maxitr, int maxprc, int ns, char *s )
{
   int i,fail;
   int a[5];
   int b[ns];
   double *c;

   a[0] = dim;
   a[1] = ns;
   a[2] = wanted;
   a[3] = maxitr;
   a[4] = maxprc;
   for(i=0; i<ns; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(179,a,b,c,0);

   return fail;
}

/* wrapping the operations from C_to_PHCpack */

int read_standard_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(11,a,b,c,0);
   return fail;
}

int read_standard_target_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(540,&n,b,c,0);

   return fail;
}

int read_dobldobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(231,a,b,c,0);
   return fail;
}

int read_dobldobl_target_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(541,&n,b,c,0);

   return fail;
}

int read_quaddobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(241,a,b,c,0);
   return fail;
}

int read_quaddobl_target_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(542,&n,b,c,0);

   return fail;
}

int read_multprec_target_system ( int decimals )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(491,&decimals,b,c,0);
   return fail;
}

int read_multprec_target_system_from_file
  ( int decimals, int n, char* filename )
{
   int b[n],i,fail;
   int a[2];
   double *c;

   a[0] = n;
   a[1] = decimals;
   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(543,a,b,c,0);

   return fail;
}

int write_standard_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(12,a,b,c,0);
   return fail;
}

int write_dobldobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(232,a,b,c,0);
   return fail;
}

int write_quaddobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(242,a,b,c,0);
   return fail;
}

int write_multprec_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(492,a,b,c,0);
   return fail;
}

int read_standard_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(13,a,b,c,0);
   return fail;
}

int read_standard_start_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(544,&n,b,c,0);

   return fail;
}

int read_dobldobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(233,a,b,c,0);
   return fail;
}

int read_dobldobl_start_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(545,&n,b,c,0);

   return fail;
}

int read_quaddobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(243,a,b,c,0);
   return fail;
}

int read_quaddobl_start_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(546,&n,b,c,0);

   return fail;
}

int read_multprec_start_system ( int decimals )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(493,&decimals,b,c,0);
   return fail;
}

int read_multprec_start_system_from_file
 ( int decimals, int n, char* filename )
{
   int b[n],i,fail;
   int a[2];
   double *c;

   a[0] = n;
   a[1] = decimals;
   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc(547,a,b,c,0);

   return fail;
}

int write_standard_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(14,a,b,c,0);
   return fail;
}

int write_dobldobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(234,a,b,c,0);
   return fail;
}

int write_quaddobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(244,a,b,c,0);
   return fail;
}

int write_multprec_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(494,a,b,c,0);
   return fail;
}

int read_standard_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(759,a,b,c,0);
   return fail;
}

int write_standard_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(760,a,b,c,0);
   return fail;
}

int read_standard_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(761,a,b,c,0);
   return fail;
}

int write_standard_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(762,a,b,c,0);
   return fail;
}

int read_dobldobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(763,a,b,c,0);
   return fail;
}

int write_dobldobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(764,a,b,c,0);
   return fail;
}

int read_dobldobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(765,a,b,c,0);
   return fail;
}

int write_dobldobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(766,a,b,c,0);
   return fail;
}

int read_quaddobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(767,a,b,c,0);
   return fail;
}

int write_quaddobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(768,a,b,c,0);
   return fail;
}

int read_quaddobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(769,a,b,c,0);
   return fail;
}

int write_quaddobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(770,a,b,c,0);
   return fail;
}

int write_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(15,a,b,c,0);
   return fail;
}

int write_dobldobl_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(235,a,b,c,0);
   return fail;
}

int write_quaddobl_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(245,a,b,c,0);
   return fail;
}

int write_multprec_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(495,a,b,c,0);
   return fail;
}

int tune_continuation_parameters ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(70,a,b,c,0);
   return fail;
}

int determine_output_during_continuation ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(71,a,b,c,0);
   return fail;
}

int retrieve_continuation_parameters ( double *c )
{
   int *a,*b,fail;
   fail = _ada_use_c2phc(72,a,b,c,0);
   return fail;
}

int set_continuation_parameters ( double *c )
{
   int *a,*b,fail;
   fail = _ada_use_c2phc(73,a,b,c,0);
   return fail;
}

int autotune_continuation_parameters
 ( int difficulty_level, int digits_of_precision )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(193,&difficulty_level,&digits_of_precision,c,0);
   return fail;
}

int show_continuation_parameters ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(194,a,b,c,0);
   return fail;
}

int get_value_of_continuation_parameter ( int k, double *val )
{
   int *b,fail;
   fail = _ada_use_c2phc(189,&k,b,val,0);
   return fail;
}

int set_value_of_continuation_parameter ( int k, double *val )
{
   int *b,fail;
   fail = _ada_use_c2phc(190,&k,b,val,0);
   return fail;
}

int create_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(152,a,b,c,0);
   return fail;
}

int create_dobldobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(172,a,b,c,0);
   return fail;
}

int create_quaddobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(182,a,b,c,0);
   return fail;
}

int create_multprec_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(522,a,b,c,0);
   return fail;
}

int create_standard_Laurent_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(792,a,b,c,0);
   return fail;
}

int create_dobldobl_Laurent_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(793,a,b,c,0);
   return fail;
}

int create_quaddobl_Laurent_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(794,a,b,c,0);
   return fail;
}

int create_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im, int pwt )
{
   int *b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc(153,&pwt,b,c,0);
   return fail;
}

int create_dobldobl_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im, int pwt )
{
   int *b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc(173,&pwt,b,c,0);
   return fail;
}

int create_quaddobl_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im, int pwt )
{
   int *b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc(183,&pwt,b,c,0);
   return fail;
}

int create_multprec_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im, int pwt )
{
   int *b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc(523,&pwt,b,c,0);
   return fail;
}

int clear_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(154,a,b,c,0);
   return fail;
}

int clear_dobldobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(174,a,b,c,0);
   return fail;
}

int clear_quaddobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(184,a,b,c,0);
   return fail;
}

int clear_multprec_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(524,a,b,c,0);
   return fail;
}

int refine_root ( int n, int *m, double *c )
{
   int *a,fail;
   int b[2];

   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc(149,a,b,c,0);
   *m = b[1];

   return fail;
}

int solve_by_standard_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(16,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_dobldobl_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(236,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_quaddobl_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(246,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_multprec_homotopy_continuation ( int decimals )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(496,&decimals,b,c,0);
   return fail;
}

int solve_by_standard_Laurent_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(774,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_dobldobl_Laurent_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(775,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_quaddobl_Laurent_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(776,&number_of_tasks,b,c,0);
   return fail;
}

int write_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(17,a,b,c,0);
   return fail;
}

int write_dobldobl_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(237,a,b,c,0);
   return fail;
}

int write_quaddobl_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(247,a,b,c,0);
   return fail;
}

int write_multprec_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(497,a,b,c,0);
   return fail;
}

int clear_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(18,a,b,c,0);
   return fail;
}

int clear_dobldobl_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(238,a,b,c,0);
   return fail;
}

int clear_quaddobl_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(248,a,b,c,0);
   return fail;
}

int clear_multprec_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(498,a,b,c,0);
   return fail;
}

int clear_standard_Laurent_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(771,a,b,c,0);
   return fail;
}

int clear_dobldobl_Laurent_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(772,a,b,c,0);
   return fail;
}

int clear_quaddobl_Laurent_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(773,a,b,c,0);
   return fail;
}

int define_output_file ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(19,a,b,c,0);
   return fail;
}

int define_output_file_with_string ( int n, char *s )
{
   int i,b[n],fail;
   double *c;
   for(i=0; i<n; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc(191,&n,b,c,0);
   return fail;
}

int close_output_file ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(192,a,b,c,0);
   return fail;
}

int write_string_to_defined_output_file ( int n, char *s )
{
   int i,b[n],fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(158,&n,b,c,0);

   return fail;
}

int write_integers_to_defined_output_file ( int n, int *a )
{
   int i,fail;
   double *c;

   fail = _ada_use_c2phc(159,&n,a,c,0);

   return fail;
}

int write_doubles_to_defined_output_file ( int n, double *a )
{
   int i,*b,fail;

   fail = _ada_use_c2phc(160,&n,b,a,0);

   return fail;
}

/* TRANSFER of data between PHCpack and the containers : */

int copy_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(1,a,b,c,0);
   return fail;
}

int copy_dobldobl_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(251,a,b,c,0);
   return fail;
}

int copy_quaddobl_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(261,a,b,c,0);
   return fail;
}

int copy_multprec_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(281,a,b,c,0);
   return fail;
}

int copy_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(2,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(252,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(262,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(282,a,b,c,0);
   return fail;
}

int copy_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(3,a,b,c,0);
   return fail;
}

int copy_dobldobl_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(253,a,b,c,0);
   return fail;
}

int copy_quaddobl_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(263,a,b,c,0);
   return fail;
}

int copy_multprec_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(283,a,b,c,0);
   return fail;
}

int copy_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(4,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(254,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(264,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(284,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(777,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(778,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(779,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(780,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(781,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(782,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(783,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(784,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(785,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(786,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(787,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(788,a,b,c,0);
   return fail;
}

int copy_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(5,a,b,c,0);
   return fail;
}

int copy_dobldobl_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(255,a,b,c,0);
   return fail;
}

int copy_quaddobl_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(265,a,b,c,0);
   return fail;
}

int copy_multprec_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(285,a,b,c,0);
   return fail;
}

int copy_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(6,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(256,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(266,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(286,a,b,c,0);
   return fail;
}

int copy_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(7,a,b,c,0);
   return fail;
}

int copy_dobldobl_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(257,a,b,c,0);
   return fail;
}

int copy_quaddobl_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(267,a,b,c,0);
   return fail;
}

int copy_multprec_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(287,a,b,c,0);
   return fail;
}

int copy_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(8,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(258,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(268,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(288,a,b,c,0);
   return fail;
}

int standard_condition_report
      ( int maxit, double tolres, double tolerr, double tolsing,
        int nbc, char *name, int *cntfail, int *cntreal, int *cntcmplx,
        int *cntregu, int *cntsing, int *cntclus,
        int *t_err, int *t_rco, int *t_res, int verbose )
{
   int fail,idx;
   int a[6];
   double c[3];

   a[0] = maxit; a[1] = verbose; a[2] = nbc;
   c[0] = tolres; c[1] = tolerr; c[2] = tolsing;

   if(nbc == 0)
   {
      int b[48];

      fail = _ada_use_c2phc(920,a,b,c,0);

      idx = 0;
      for(int i=0; i<16; i++) t_err[i] = b[idx++];
      for(int i=0; i<16; i++) t_rco[i] = b[idx++];
      for(int i=0; i<16; i++) t_res[i] = b[idx++];
   }
   else if(nbc > 48)
   {
      int b[nbc];

      for(int i=0; i<nbc; i++) b[i] = (int) name[i];

      fail = _ada_use_c2phc(920,a,b,c,0);

      idx = 0;
      for(int i=0; i<16; i++) t_err[i] = b[idx++];
      for(int i=0; i<16; i++) t_rco[i] = b[idx++];
      for(int i=0; i<16; i++) t_res[i] = b[idx++];
   }
   else
   {
      int b[48];

      for(int i=0; i<nbc; i++) b[i] = (int) name[i];

      fail = _ada_use_c2phc(920,a,b,c,0);

      idx = 0;
      for(int i=0; i<16; i++) t_err[i] = b[idx++];
      for(int i=0; i<16; i++) t_rco[i] = b[idx++];
      for(int i=0; i<16; i++) t_res[i] = b[idx++];
   }
   *cntfail = a[0];
   *cntreal = a[1];
   *cntcmplx = a[2];
   *cntregu = a[3];
   *cntsing = a[4];
   *cntclus = a[5];

   return fail;
}

/* OPERATIONS on data in the containers : */

int validate_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(9,a,b,c,0);
   return fail;
}

int print_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(21,a,b,c,0);
   return fail;
}
