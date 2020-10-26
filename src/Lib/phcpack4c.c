/* file phcpack4c.c contains the definitions of the functions 
 * in phcpack4c.h */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "phcpack4c.h"

#define v 0 /* verbose flag */

/* most BASIC operations in PHCpack : */

int version_string ( int *n, char *s )
{
   int fail,i;
   int b[40];
   double *c;

   fail = _ada_use_c2phc4c(999,n,b,c,0);

   for(i=0; i<(*n); i++) s[i] = (char) b[i];
   s[*n] = '\0';

   return fail;
}

int set_seed ( int seed )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(998,&seed,b,c,0);

   return fail;
}

int get_seed ( int *seed )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(997,seed,b,c,0);

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

   fail = _ada_use_c2phc4c(77,a,b,c,vrb);

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

   fail = _ada_use_c2phc4c(700,a,b,c,vrb);

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

   fail = _ada_use_c2phc4c(702,a,b,c,vrb);

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

   fail = _ada_use_c2phc4c(75,a,b,c,vrb);

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

   fail = _ada_use_c2phc4c(701,a,b,c,vrb);

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

   fail = _ada_use_c2phc4c(703,a,b,c,vrb);

   *root_count = a[0];

   if(silent == 0)
   {
      *nrcs = a[1];
      for(i=0; i<(*nrcs); i++) rocos[i] = (char) b[i];
      rocos[*nrcs] = '\0';
   }
   return fail;
}

int mixed_volume ( int *mv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(78,mv,b,c,0);

   return fail;
}

int stable_mixed_volume ( int *mv, int *smv )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(79,mv,smv,c,0);

   return fail;
}

int mixed_volume_by_demics ( int *mv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(843,mv,b,c,0);

   return fail;
}

int stable_mixed_volume_by_demics ( int *mv, int *smv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(844,mv,smv,c,0);

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

   fail = _ada_use_c2phc4c(196,&maxitr,&maxdef,tols,0);

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

   fail = _ada_use_c2phc4c(249,&maxitr,&maxdef,tols,0);

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

   fail = _ada_use_c2phc4c(250,&maxitr,&maxdef,tols,0);

   return fail;
}

int standard_Newton_step ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(199,a,b,c,0);

   return fail;
}

int dobldobl_Newton_step ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(198,a,b,c,0);

   return fail;
}

int quaddobl_Newton_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(197,a,b,c,0);
   return fail;
}

int multprec_Newton_step ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(195,&deci,b,c,0);
   return fail;
}

int standard_Newton_Laurent_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(326,a,b,c,0);
   return fail;
}

int dobldobl_Newton_Laurent_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(327,a,b,c,0);
   return fail;
}

int quaddobl_Newton_Laurent_step ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(328,a,b,c,0);
   return fail;
}

int multprec_Newton_Laurent_step ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(329,&deci,b,c,0);
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

   fail = _ada_use_c2phc4c(439,dim,b,c,0);

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

   fail = _ada_use_c2phc4c(179,a,b,c,0);

   return fail;
}

/* wrapping the operations from C_to_PHCpack */

int read_standard_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(11,a,b,c,0);
   return fail;
}

int read_standard_target_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc4c(540,&n,b,c,0);

   return fail;
}

int read_dobldobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(231,a,b,c,0);
   return fail;
}

int read_dobldobl_target_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc4c(541,&n,b,c,0);

   return fail;
}

int read_quaddobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(241,a,b,c,0);
   return fail;
}

int read_quaddobl_target_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc4c(542,&n,b,c,0);

   return fail;
}

int read_multprec_target_system ( int decimals )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(491,&decimals,b,c,0);
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

   fail = _ada_use_c2phc4c(543,a,b,c,0);

   return fail;
}

int write_standard_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(12,a,b,c,0);
   return fail;
}

int write_dobldobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(232,a,b,c,0);
   return fail;
}

int write_quaddobl_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(242,a,b,c,0);
   return fail;
}

int write_multprec_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(492,a,b,c,0);
   return fail;
}

int read_standard_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(13,a,b,c,0);
   return fail;
}

int read_standard_start_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc4c(544,&n,b,c,0);

   return fail;
}

int read_dobldobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(233,a,b,c,0);
   return fail;
}

int read_dobldobl_start_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc4c(545,&n,b,c,0);

   return fail;
}

int read_quaddobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(243,a,b,c,0);
   return fail;
}

int read_quaddobl_start_system_from_file ( int n, char* filename )
{
   int b[n],i,fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) filename[i];

   fail = _ada_use_c2phc4c(546,&n,b,c,0);

   return fail;
}

int read_multprec_start_system ( int decimals )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(493,&decimals,b,c,0);
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

   fail = _ada_use_c2phc4c(547,a,b,c,0);

   return fail;
}

int write_standard_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(14,a,b,c,0);
   return fail;
}

int write_dobldobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(234,a,b,c,0);
   return fail;
}

int write_quaddobl_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(244,a,b,c,0);
   return fail;
}

int write_multprec_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(494,a,b,c,0);
   return fail;
}

int read_standard_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(759,a,b,c,0);
   return fail;
}

int write_standard_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(760,a,b,c,0);
   return fail;
}

int read_standard_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(761,a,b,c,0);
   return fail;
}

int write_standard_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(762,a,b,c,0);
   return fail;
}

int read_dobldobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(763,a,b,c,0);
   return fail;
}

int write_dobldobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(764,a,b,c,0);
   return fail;
}

int read_dobldobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(765,a,b,c,0);
   return fail;
}

int write_dobldobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(766,a,b,c,0);
   return fail;
}

int read_quaddobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(767,a,b,c,0);
   return fail;
}

int write_quaddobl_start_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(768,a,b,c,0);
   return fail;
}

int read_quaddobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(769,a,b,c,0);
   return fail;
}

int write_quaddobl_target_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(770,a,b,c,0);
   return fail;
}

int write_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(15,a,b,c,0);
   return fail;
}

int write_dobldobl_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(235,a,b,c,0);
   return fail;
}

int write_quaddobl_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(245,a,b,c,0);
   return fail;
}

int write_multprec_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(495,a,b,c,0);
   return fail;
}

int tune_continuation_parameters ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(70,a,b,c,0);
   return fail;
}

int determine_output_during_continuation ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(71,a,b,c,0);
   return fail;
}

int retrieve_continuation_parameters ( double *c )
{
   int *a,*b,fail;
   fail = _ada_use_c2phc4c(72,a,b,c,0);
   return fail;
}

int set_continuation_parameters ( double *c )
{
   int *a,*b,fail;
   fail = _ada_use_c2phc4c(73,a,b,c,0);
   return fail;
}

int autotune_continuation_parameters
 ( int difficulty_level, int digits_of_precision )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(193,&difficulty_level,&digits_of_precision,c,0);
   return fail;
}

int show_continuation_parameters ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(194,a,b,c,0);
   return fail;
}

int get_value_of_continuation_parameter ( int k, double *val )
{
   int *b,fail;
   fail = _ada_use_c2phc4c(189,&k,b,val,0);
   return fail;
}

int set_value_of_continuation_parameter ( int k, double *val )
{
   int *b,fail;
   fail = _ada_use_c2phc4c(190,&k,b,val,0);
   return fail;
}

int create_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(152,a,b,c,0);
   return fail;
}

int create_dobldobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(172,a,b,c,0);
   return fail;
}

int create_quaddobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(182,a,b,c,0);
   return fail;
}

int create_multprec_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(522,a,b,c,0);
   return fail;
}

int create_standard_Laurent_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(792,a,b,c,0);
   return fail;
}

int create_dobldobl_Laurent_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(793,a,b,c,0);
   return fail;
}

int create_quaddobl_Laurent_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(794,a,b,c,0);
   return fail;
}

int create_homotopy_with_given_gamma ( double gamma_re, double gamma_im )
{
   int *a,*b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc4c(153,a,b,c,0);
   return fail;
}

int create_dobldobl_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im )
{
   int *a,*b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc4c(173,a,b,c,0);
   return fail;
}

int create_quaddobl_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im )
{
   int *a,*b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc4c(183,a,b,c,0);
   return fail;
}

int create_multprec_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im )
{
   int *a,*b,fail;
   double c[2];
   c[0] = gamma_re; c[1] = gamma_im;
   fail = _ada_use_c2phc4c(523,a,b,c,0);
   return fail;
}

int clear_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(154,a,b,c,0);
   return fail;
}

int clear_dobldobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(174,a,b,c,0);
   return fail;
}

int clear_quaddobl_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(184,a,b,c,0);
   return fail;
}

int clear_multprec_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(524,a,b,c,0);
   return fail;
}

int refine_root ( int n, int *m, double *c )
{
   int *a,fail;
   int b[2];

   b[0] = n; b[1] = *m;
   fail = _ada_use_c2phc4c(149,a,b,c,0);
   *m = b[1];

   return fail;
}

int solve_by_standard_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(16,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_dobldobl_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(236,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_quaddobl_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(246,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_multprec_homotopy_continuation ( int decimals )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(496,&decimals,b,c,0);
   return fail;
}

int solve_by_standard_Laurent_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(774,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_dobldobl_Laurent_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(775,&number_of_tasks,b,c,0);
   return fail;
}

int solve_by_quaddobl_Laurent_homotopy_continuation ( int number_of_tasks )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(776,&number_of_tasks,b,c,0);
   return fail;
}

int write_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(17,a,b,c,0);
   return fail;
}

int write_dobldobl_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(237,a,b,c,0);
   return fail;
}

int write_quaddobl_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(247,a,b,c,0);
   return fail;
}

int write_multprec_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(497,a,b,c,0);
   return fail;
}

int clear_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(18,a,b,c,0);
   return fail;
}

int clear_dobldobl_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(238,a,b,c,0);
   return fail;
}

int clear_quaddobl_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(248,a,b,c,0);
   return fail;
}

int clear_multprec_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(498,a,b,c,0);
   return fail;
}

int clear_standard_Laurent_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(771,a,b,c,0);
   return fail;
}

int clear_dobldobl_Laurent_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(772,a,b,c,0);
   return fail;
}

int clear_quaddobl_Laurent_data ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(773,a,b,c,0);
   return fail;
}

int define_output_file ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(19,a,b,c,0);
   return fail;
}

int define_output_file_with_string ( int n, char *s )
{
   int i,b[n],fail;
   double *c;
   for(i=0; i<n; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(191,&n,b,c,0);
   return fail;
}

int close_output_file ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(192,a,b,c,0);
   return fail;
}

int write_string_to_defined_output_file ( int n, char *s )
{
   int i,b[n],fail;
   double *c;

   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(158,&n,b,c,0);

   return fail;
}

int write_integers_to_defined_output_file ( int n, int *a )
{
   int i,fail;
   double *c;

   fail = _ada_use_c2phc4c(159,&n,a,c,0);

   return fail;
}

int write_doubles_to_defined_output_file ( int n, double *a )
{
   int i,*b,fail;

   fail = _ada_use_c2phc4c(160,&n,b,a,0);

   return fail;
}

/* TRANSFER of data between PHCpack and the containers : */

int copy_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(1,a,b,c,0);
   return fail;
}

int copy_dobldobl_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(251,a,b,c,0);
   return fail;
}

int copy_quaddobl_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(261,a,b,c,0);
   return fail;
}

int copy_multprec_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(281,a,b,c,0);
   return fail;
}

int copy_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(2,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(252,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(262,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(282,a,b,c,0);
   return fail;
}

int copy_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(3,a,b,c,0);
   return fail;
}

int copy_dobldobl_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(253,a,b,c,0);
   return fail;
}

int copy_quaddobl_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(263,a,b,c,0);
   return fail;
}

int copy_multprec_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(283,a,b,c,0);
   return fail;
}

int copy_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(4,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(254,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(264,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(284,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(777,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(778,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_container_to_start_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(779,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(780,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(781,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_container_to_target_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(782,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(783,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(784,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_start_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(785,a,b,c,0);
   return fail;
}

int copy_standard_Laurent_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(786,a,b,c,0);
   return fail;
}

int copy_dobldobl_Laurent_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(787,a,b,c,0);
   return fail;
}

int copy_quaddobl_Laurent_target_system_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(788,a,b,c,0);
   return fail;
}

int copy_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(5,a,b,c,0);
   return fail;
}

int copy_dobldobl_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(255,a,b,c,0);
   return fail;
}

int copy_quaddobl_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(265,a,b,c,0);
   return fail;
}

int copy_multprec_target_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(285,a,b,c,0);
   return fail;
}

int copy_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(6,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(256,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(266,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_target_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(286,a,b,c,0);
   return fail;
}

int copy_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(7,a,b,c,0);
   return fail;
}

int copy_dobldobl_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(257,a,b,c,0);
   return fail;
}

int copy_quaddobl_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(267,a,b,c,0);
   return fail;
}

int copy_multprec_start_solutions_to_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(287,a,b,c,0);
   return fail;
}

int copy_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(8,a,b,c,0);
   return fail;
}

int copy_dobldobl_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(258,a,b,c,0);
   return fail;
}

int copy_quaddobl_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(268,a,b,c,0);
   return fail;
}

int copy_multprec_container_to_start_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(288,a,b,c,0);
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

      fail = _ada_use_c2phc4c(920,a,b,c,0);

      idx = 0;
      for(int i=0; i<16; i++) t_err[i] = b[idx++];
      for(int i=0; i<16; i++) t_rco[i] = b[idx++];
      for(int i=0; i<16; i++) t_res[i] = b[idx++];
   }
   else if(nbc > 48)
   {
      int b[nbc];

      for(int i=0; i<nbc; i++) b[i] = (int) name[i];

      fail = _ada_use_c2phc4c(920,a,b,c,0);

      idx = 0;
      for(int i=0; i<16; i++) t_err[i] = b[idx++];
      for(int i=0; i<16; i++) t_rco[i] = b[idx++];
      for(int i=0; i<16; i++) t_res[i] = b[idx++];
   }
   else
   {
      int b[48];

      for(int i=0; i<nbc; i++) b[i] = (int) name[i];

      fail = _ada_use_c2phc4c(920,a,b,c,0);

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
   fail = _ada_use_c2phc4c(9,a,b,c,0);
   return fail;
}

int print_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(21,a,b,c,0);
   return fail;
}
