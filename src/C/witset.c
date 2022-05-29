/* file witset.c contains the definitions of the functions in witset.h */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phcpack.h"
#include "syscon.h"
#include "solcon.h"
#include "witset.h"

#define verbose 0 /* verbose flag */

/* some basic OPERATIONS on witness sets */

int embed_system ( int d, int precision, int vrblvl )
{
   int fail;

   if(precision == 0) fail = embed_standard_system(d,vrblvl);
   if(precision == 1) fail = embed_dobldobl_system(d,vrblvl);
   if(precision == 2) fail = embed_quaddobl_system(d,vrblvl);

   return fail;
}

int embed_standard_system ( int d, int vrblvl )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(66,&d,b,c,vrblvl);
   return fail;
}

int embed_dobldobl_system ( int d, int vrblvl )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(129,&d,b,c,vrblvl);
   return fail;
}

int embed_quaddobl_system ( int d, int vrblvl )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(260,&d,b,c,vrblvl);
   return fail;
}

int embed_standard_Laurent_system ( int d, int vrblvl )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(625,&d,b,c,vrblvl);
   return fail;
}

int embed_dobldobl_Laurent_system ( int d, int vrblvl )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(626,&d,b,c,vrblvl);
   return fail;
}

int embed_quaddobl_Laurent_system ( int d, int vrblvl )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(627,&d,b,c,vrblvl);
   return fail;
}

int read_witness_set ( int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

   fail = _ada_use_c2phc4c(41,n,d,c,0);    /* read the witness set */
   *dim = d[0];
   *deg = d[1];
   if(verbose>0) printf("The ambient dimension is %d.\n", *n); 
   if(verbose>0) printf("The dimension of the solution set is %d.\n", d[0]); 
   if(verbose>0) printf("The degree of the solution set is %d.\n", d[1]); 
   if(verbose>1)
   {
      printf("The embedded system :\n");
      fail = syscon_write_standard_system();
   }
   if(verbose>1)
   {
      printf("The following solutions are in container :\n");
      fail = solcon_write_standard_solutions();
   }
   return fail;
}

int read_dobldobl_witness_set ( int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

   fail = _ada_use_c2phc4c(631,n,d,c,0);    /* read the witness set */
   *dim = d[0];
   *deg = d[1];
   if(verbose>0) printf("The ambient dimension is %d.\n", *n); 
   if(verbose>0) printf("The dimension of the solution set is %d.\n", d[0]); 
   if(verbose>0) printf("The degree of the solution set is %d.\n", d[1]); 
   if(verbose>1)
   {
      printf("The embedded system :\n");
      fail = syscon_write_dobldobl_system();
   }
   if(verbose>1)
   {
      printf("The following solutions are in container :\n");
      fail = solcon_write_dobldobl_solutions();
   }
   return fail;
}

int read_quaddobl_witness_set ( int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

   fail = _ada_use_c2phc4c(661,n,d,c,0);    /* read the witness set */
   *dim = d[0];
   *deg = d[1];
   if(verbose>0) printf("The ambient dimension is %d.\n", *n); 
   if(verbose>0) printf("The dimension of the solution set is %d.\n", d[0]); 
   if(verbose>0) printf("The degree of the solution set is %d.\n", d[1]); 
   if(verbose>1)
   {
      printf("The embedded system :\n");
      fail = syscon_write_quaddobl_system();
   }
   if(verbose>1)
   {
      printf("The following solutions are in container :\n");
      fail = solcon_write_quaddobl_solutions();
   }
   return fail;
}

int read_witness_set_from_file ( int m, char *s, int *n, int *dim, int *deg )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   *n = m;
   fail = _ada_use_c2phc4c(64,n,b,c,0);
   *dim = b[0];
   *deg = b[1];

   return fail;
}

int read_dobldobl_witness_set_from_file
 ( int m, char *s, int *n, int *dim, int *deg )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   *n = m;
   fail = _ada_use_c2phc4c(654,n,b,c,0);
   *dim = b[0];
   *deg = b[1];

   return fail;
}

int read_quaddobl_witness_set_from_file
 ( int m, char *s, int *n, int *dim, int *deg )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   *n = m;
   fail = _ada_use_c2phc4c(684,n,b,c,0);
   *dim = b[0];
   *deg = b[1];

   return fail;
}

int read_standard_Laurent_witness_set ( int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

   fail = _ada_use_c2phc4c(798,n,d,c,0);    /* read the witness set */
   *dim = d[0];
   *deg = d[1];
   if(verbose>0) printf("The ambient dimension is %d.\n", *n); 
   if(verbose>0) printf("The dimension of the solution set is %d.\n", d[0]); 
   if(verbose>0) printf("The degree of the solution set is %d.\n", d[1]); 
   if(verbose>1)
   {
      printf("The embedded system :\n");
      fail = syscon_write_standard_system();
   }
   if(verbose>1)
   {
      printf("The following solutions are in container :\n");
      fail = solcon_write_standard_solutions();
   }
   return fail;
}

int read_dobldobl_Laurent_witness_set ( int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

   fail = _ada_use_c2phc4c(799,n,d,c,0);    /* read the witness set */
   *dim = d[0];
   *deg = d[1];
   if(verbose>0) printf("The ambient dimension is %d.\n", *n); 
   if(verbose>0) printf("The dimension of the solution set is %d.\n", d[0]); 
   if(verbose>0) printf("The degree of the solution set is %d.\n", d[1]); 
   if(verbose>1)
   {
      printf("The embedded system :\n");
      fail = syscon_write_dobldobl_system();
   }
   if(verbose>1)
   {
      printf("The following solutions are in container :\n");
      fail = solcon_write_dobldobl_solutions();
   }
   return fail;
}

int read_quaddobl_Laurent_witness_set ( int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

   fail = _ada_use_c2phc4c(800,n,d,c,0);    /* read the witness set */
   *dim = d[0];
   *deg = d[1];
   if(verbose>0) printf("The ambient dimension is %d.\n", *n); 
   if(verbose>0) printf("The dimension of the solution set is %d.\n", d[0]); 
   if(verbose>0) printf("The degree of the solution set is %d.\n", d[1]); 
   if(verbose>1)
   {
      printf("The embedded system :\n");
      fail = syscon_write_quaddobl_system();
   }
   if(verbose>1)
   {
      printf("The following solutions are in container :\n");
      fail = solcon_write_quaddobl_solutions();
   }
   return fail;
}

int read_standard_Laurent_witness_set_from_file
 ( int m, char *s, int *n, int *dim, int *deg )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   *n = m;
   fail = _ada_use_c2phc4c(801,n,b,c,0);
   *dim = b[0];
   *deg = b[1];

   return fail;
}

int read_dobldobl_Laurent_witness_set_from_file
 ( int m, char *s, int *n, int *dim, int *deg )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   *n = m;
   fail = _ada_use_c2phc4c(802,n,b,c,0);
   *dim = b[0];
   *deg = b[1];

   return fail;
}

int read_quaddobl_Laurent_witness_set_from_file
 ( int m, char *s, int *n, int *dim, int *deg )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   *n = m;
   fail = _ada_use_c2phc4c(803,n,b,c,0);
   *dim = b[0];
   *deg = b[1];

   return fail;
}

int write_witness_set_to_file ( int m, char *s )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(65,&m,b,c,0);

   return fail;
}

int write_dobldobl_witness_set_to_file ( int m, char *s )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(655,&m,b,c,0);

   return fail;
}

int write_quaddobl_witness_set_to_file ( int m, char *s )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(685,&m,b,c,0);

   return fail;
}

int read_a_witness_set ( int k, int *n, int *dim, int *deg )
{
   int fail;
   double *c;
   int d[2];

   *n = k;
   fail = _ada_use_c2phc4c(166,n,d,c,0);
   *dim = d[0];
   *deg = d[1];

   return fail;
}

int standard_witness_set_to_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(47,a,b,c,0);   /* copy system to container */
   return fail;
}

int dobldobl_witness_set_to_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(637,a,b,c,0);   /* copy system to container */
   return fail;
}

int quaddobl_witness_set_to_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(667,a,b,c,0);   /* copy system to container */
   return fail;
}

int standard_witness_set_to_Laurent_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(807,a,b,c,0);   /* copy system to container */
   return fail;
}

int dobldobl_witness_set_to_Laurent_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(808,a,b,c,0);   /* copy system to container */
   return fail;
}

int quaddobl_witness_set_to_Laurent_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(809,a,b,c,0);   /* copy system to container */
   return fail;
}

int swap_symbols_for_standard_witness_set ( int nvr, int dim )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(816,&nvr,&dim,c,0);
   return fail;
}

int swap_symbols_for_dobldobl_witness_set ( int nvr, int dim )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(817,&nvr,&dim,c,0);
   return fail;
}

int swap_symbols_for_quaddobl_witness_set ( int nvr, int dim )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(818,&nvr,&dim,c,0);
   return fail;
}

int swap_symbols_for_standard_Laurent_witness_set ( int nvr, int dim )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(819,&nvr,&dim,c,0);
   return fail;
}

int swap_symbols_for_dobldobl_Laurent_witness_set ( int nvr, int dim )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(820,&nvr,&dim,c,0);
   return fail;
}

int swap_symbols_for_quaddobl_Laurent_witness_set ( int nvr, int dim )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(821,&nvr,&dim,c,0);
   return fail;
}

int create_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(164,a,b,c,0);
   return fail;
}

int create_dobldobl_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(178,a,b,c,0);
   return fail;
}

int create_quaddobl_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(188,a,b,c,0);
   return fail;
}

int create_standard_Laurent_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(789,a,b,c,0);
   return fail;
}

int create_dobldobl_Laurent_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(790,a,b,c,0);
   return fail;
}

int create_quaddobl_Laurent_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
    fail = _ada_use_c2phc4c(791,a,b,c,0);
   return fail;
}

/* OPERATIONS to intersect witness sets */

int standard_diagonal_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(165,&a,&b,c,0);
   return fail;
}

int dobldobl_diagonal_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(289,&a,&b,c,0);
   return fail;
}

int quaddobl_diagonal_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(290,&a,&b,c,0);
   return fail;
}

int standard_diagonal_Laurent_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(810,&a,&b,c,0);
   return fail;
}

int dobldobl_diagonal_Laurent_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(811,&a,&b,c,0);
   return fail;
}

int quaddobl_diagonal_Laurent_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(812,&a,&b,c,0);
   return fail;
}

int standard_diagonal_cascade_solutions ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(271,&a,&b,c,0);
   return fail;
}

int dobldobl_diagonal_cascade_solutions ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(297,&a,&b,c,0);
   return fail;
}

int quaddobl_diagonal_cascade_solutions ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(298,&a,&b,c,0);
   return fail;
}

int extrinsic_top_diagonal_dimension
 ( int n1, int n2, int a, int b, int *d )
{
   int fail,alpha[2],beta[2];
   double *c;
   alpha[0] = n1; alpha[1] = n2;
   beta[0] = a; beta[1] = b;

   fail = _ada_use_c2phc4c(168,alpha,beta,c,0);

   *d = alpha[0];

   return fail;
}

int hypersurface_witness_set ( int k, int n, char *s )
{
   int a[2],b[n],i,fail;
   double *c;

   a[0] = k;
   a[1] = n;
   for(i=0; i<n; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc4c(169,a,b,c,0);
   
   return fail;
}

int standard_witset_of_hypersurface ( int nv, int nc, char *p )
{
   int i,fail,a[2],b[nc];
   double *c;

   /* printf("nv = %d, nc = %d\n",nv,nc); */
   /* printf("p = %s\n",p); */
   a[0] = nv;
   a[1] = nc;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(270,a,b,c,0);

   return fail;
}

int dobldobl_witset_of_hypersurface ( int nv, int nc, char *p )
{
   int i,fail,a[2],b[nc];
   double *c;

   /* printf("nv = %d, nc = %d\n",nv,nc); */
   /* printf("p = %s\n",p); */
   a[0] = nv;
   a[1] = nc;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(259,a,b,c,0);

   return fail;
}

int quaddobl_witset_of_hypersurface ( int nv, int nc, char *p )
{
   int i,fail,a[2],b[nc];
   double *c;

   /* printf("nv = %d, nc = %d\n",nv,nc); */
   /* printf("p = %s\n",p); */
   a[0] = nv;
   a[1] = nc;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(269,a,b,c,0);

   return fail;
}

int standard_witset_of_Laurent_hypersurface ( int nv, int nc, char *p )
{
   int i,fail,a[2],b[nc];
   double *c;

   /* printf("nv = %d, nc = %d\n",nv,nc); */
   /* printf("p = %s\n",p); */
   a[0] = nv;
   a[1] = nc;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(813,a,b,c,0);

   return fail;
}

int dobldobl_witset_of_Laurent_hypersurface ( int nv, int nc, char *p )
{
   int i,fail,a[2],b[nc];
   double *c;

   /* printf("nv = %d, nc = %d\n",nv,nc); */
   /* printf("p = %s\n",p); */
   a[0] = nv;
   a[1] = nc;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(814,a,b,c,0);

   return fail;
}

int quaddobl_witset_of_Laurent_hypersurface ( int nv, int nc, char *p )
{
   int i,fail,a[2],b[nc];
   double *c;

   /* printf("nv = %d, nc = %d\n",nv,nc); */
   /* printf("p = %s\n",p); */
   a[0] = nv;
   a[1] = nc;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(815,a,b,c,0);

   return fail;
}

int diagonal_symbols_doubler ( int n, int d, int nc, char *s )
{
   int k,fail = 0;
   int a[3];
   int b[nc];
   double *c;

   a[0] = n;
   a[1] = d;
   a[2] = nc;

   for(k=0; k<nc; k++) b[k] = (int) s[k];
 
   fail = _ada_use_c2phc4c(230,a,b,c,0);

   return fail;
}

int standard_collapse_diagonal ( int k, int d )
{
   int a[2],*b,fail;
   double *c;

   a[0] = k; a[1] = d;
   fail = _ada_use_c2phc4c(170,a,b,c,0);

   return fail;
}

int dobldobl_collapse_diagonal ( int k, int d )
{
   int a[2],*b,fail;
   double *c;

   a[0] = k; a[1] = d;
   fail = _ada_use_c2phc4c(299,a,b,c,0);

   return fail;
}

int quaddobl_collapse_diagonal ( int k, int d )
{
   int a[2],*b,fail;
   double *c;

   a[0] = k; a[1] = d;
   fail = _ada_use_c2phc4c(312,a,b,c,0);

   return fail;
}

int remove_last_slack ( int k )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(171,&k,b,c,0);

   return fail;
}

/* OPERATIONS needed in the monodromy factorization */

int list2str ( int n, int *d, char *s )
{
   int bufsize = 16;
   int cnt = 0;
   int i,j;

   s[cnt++] = '[';
   for(i=0; i<n; i++)
   {
      char buf[bufsize];
      for(j=0; j<bufsize; j++) buf[j] = ' ';
      sprintf(buf,"%d",d[i]);
      for(j=0; j<bufsize; j++)
      {
         if((buf[j] == '\0') || (buf[j] == ' ')) break;
         s[cnt++] = buf[j];
      }
      if(i < n-1) s[cnt++] = ',';
   }
   s[cnt++] = ']';
   s[cnt++] = '\0';

   return cnt;
}

int str2list ( int n, char *s, int *d )
{
   int nb;
   int i=0;
   int ind=0;

   while(i < n) if(s[i++] == '[') break;
   while(i < n)
   {
      sscanf(&s[i],"%d",&nb);
      d[ind++] = nb;
      while(i<n)       /* skip the number */
      {
         if(s[i] == ',') break;
         if(s[i] == ']') break;
         i = i + 1;
      }
      i = i + 1;       /* skip ',' or ']' */
   }
   return ind;
}

int set_standard_state_to_silent ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(39,a,b,c,0);

   return fail;
}

int set_dobldobl_state_to_silent ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(658,a,b,c,0);

   return fail;
}

int set_quaddobl_state_to_silent ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(688,a,b,c,0);

   return fail;
}

int set_standard_state_to_verbose ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(630,a,b,c,0);

   return fail;
}

int set_dobldobl_state_to_verbose ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(660,a,b,c,0);

   return fail;
}

int set_quaddobl_state_to_verbose ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc4c(690,a,b,c,0);

   return fail;
}

int assign_labels ( int n, int nbsols, int precision )
{
   if(precision == 0) return standard_assign_labels(n,nbsols);
   if(precision == 1) return dobldobl_assign_labels(n,nbsols);
   if(precision == 2) return quaddobl_assign_labels(n,nbsols);

   return -1;
}

int standard_assign_labels ( int n, int nbsols )
{
   int i,j,m,fail;
   double x[2*n+5];

   for(i=1; i<=nbsols; i++)
   {
      fail = solcon_retrieve_standard_solution(n,i,&m,x);
      m = i;
      fail = solcon_replace_standard_solution(n,i,m,x);
   }

   return fail;
}

int dobldobl_assign_labels ( int n, int nbsols )
{
   int i,j,m,fail;
   double x[4*n+10];

   for(i=1; i<=nbsols; i++)
   {
      fail = solcon_retrieve_dobldobl_solution(n,i,&m,x);
      m = i;
      fail = solcon_replace_dobldobl_solution(n,i,m,x);
   }

   return fail;
}

int quaddobl_assign_labels ( int n, int nbsols )
{
   int i,j,m,fail;
   double x[8*n+20];

   for(i=1; i<=nbsols; i++)
   {
      fail = solcon_retrieve_quaddobl_solution(n,i,&m,x);
      m = i;
      fail = solcon_replace_quaddobl_solution(n,i,m,x);
   }

   return fail;
}

int initialize_standard_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(42,&dim,b,c,0);        /* initialize sampler */

   return fail;
}

int initialize_dobldobl_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(632,&dim,b,c,0);        /* initialize sampler */

   return fail;
}

int initialize_quaddobl_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(662,&dim,b,c,0);        /* initialize sampler */

   return fail;
}

int initialize_standard_Laurent_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(804,&dim,b,c,0);        /* initialize sampler */

   return fail;
}

int initialize_dobldobl_Laurent_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(805,&dim,b,c,0);        /* initialize sampler */

   return fail;
}

int initialize_quaddobl_Laurent_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(806,&dim,b,c,0);        /* initialize sampler */

   return fail;
}

int initialize_standard_monodromy ( int n, int d, int k )
{
   int fail;
   double *c;
   int dk[2];

   dk[0] = d;
   dk[1] = k;
   fail = _ada_use_c2phc4c(50,&n,dk,c,0); 
   /* initialize Monodromy_Permutations */

   return fail;
}

int initialize_dobldobl_monodromy ( int n, int d, int k )
{
   int fail;
   double *c;
   int dk[2];

   dk[0] = d;
   dk[1] = k;
   fail = _ada_use_c2phc4c(640,&n,dk,c,0);

   return fail;
}

int initialize_quaddobl_monodromy ( int n, int d, int k )
{
   int fail;
   double *c;
   int dk[2];

   dk[0] = d;
   dk[1] = k;
   fail = _ada_use_c2phc4c(670,&n,dk,c,0);

   return fail;
}

int standard_trace_grid_diagnostics ( double *err, double *dis )
{
   int *a,*b,fail;
   double c[2];

   fail = _ada_use_c2phc4c(56,a,b,c,0);
   *err = c[0];
   *dis = c[1];

   return fail;
}

int dobldobl_trace_grid_diagnostics ( double *err, double *dis )
{
   int *a,*b,fail;
   double c[2];

   fail = _ada_use_c2phc4c(646,a,b,c,0);
   *err = c[0];
   *dis = c[1];

   return fail;
}

int quaddobl_trace_grid_diagnostics ( double *err, double *dis )
{
   int *a,*b,fail;
   double c[2];

   fail = _ada_use_c2phc4c(676,a,b,c,0);
   *err = c[0];
   *dis = c[1];

   return fail;
}

void random_complex ( double *re, double *im )
{
   double angle = ((double) rand())/RAND_MAX;
   *re = cos(angle);
   *im = sin(angle);
}

int random_standard_complex ( double *re, double *im )
{
   int *a,*b,fail;
   double c[2];
   
   fail = _ada_use_c2phc4c(280,a,b,c,0);
   *re = c[0];
   *im = c[1];

   return fail;
}

int random_dobldobl_complex ( double *re, double *im )
{
   int *a,*b,fail;
   double c[4];
   
   fail = _ada_use_c2phc4c(659,a,b,c,0);
   re[0] = c[0]; re[1] = c[1];
   im[0] = c[2]; im[1] = c[3];

   return fail;
}

int random_quaddobl_complex ( double *re, double *im )
{
   int *a,*b,fail;
   double c[8];
   
   fail = _ada_use_c2phc4c(689,a,b,c,0);
   re[0] = c[0]; re[1] = c[1]; re[2] = c[2]; re[3] = c[3];
   im[0] = c[4]; im[1] = c[5]; im[2] = c[6]; im[3] = c[7];

   return fail;
}

int store_standard_gamma ( int n, double *re_gamma, double *im_gamma )
{
   int fail,i;
   int *b;
   double r[2];

   for(i=0; i<n; i++)
   {
      r[0] = re_gamma[i];
      r[1] = im_gamma[i];
      fail = _ada_use_c2phc4c(44,&i,b,r,0);  /* store gamma */
   }

   return fail;
}

int store_dobldobl_gamma ( int n, double *re_gamma, double *im_gamma )
{
   int fail,i;
   int *b;
   double r[4];

   for(i=0; i<n; i++)
   {
      r[0] = re_gamma[2*i]; r[1] = re_gamma[2*i+1];
      r[2] = im_gamma[2*i]; r[3] = im_gamma[2*i+1];
      fail = _ada_use_c2phc4c(634,&i,b,r,0);  /* store gamma */
   }

   return fail;
}

int store_quaddobl_gamma ( int n, double *re_gamma, double *im_gamma )
{
   int fail,i;
   int *b;
   double r[8];

   for(i=0; i<n; i++)
   {
      r[0] = re_gamma[4*i];   r[1] = re_gamma[4*i+1];
      r[2] = re_gamma[4*i+2]; r[3] = re_gamma[4*i+3];
      r[4] = im_gamma[4*i];   r[5] = im_gamma[4*i+1];
      r[6] = im_gamma[4*i+2]; r[7] = im_gamma[4*i+3];
      fail = _ada_use_c2phc4c(664,&i,b,r,0);  /* store gamma */
   }

   return fail;
}

int assign_standard_coefficient_of_slice ( int i, int j, double *r )
{
   int fail;
   fail = _ada_use_c2phc4c(43,&i,&j,r,0);
   return fail;
}

int assign_dobldobl_coefficient_of_slice ( int i, int j, double *r )
{
   int fail;
   fail = _ada_use_c2phc4c(633,&i,&j,r,0);
   return fail;
}

int assign_quaddobl_coefficient_of_slice ( int i, int j, double *r )
{
   int fail;
   fail = _ada_use_c2phc4c(663,&i,&j,r,0);
   return fail;
}

int initialize_hyperplane_sections ( int m )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(59,&m,b,c,0);
   return fail;
}

int store_new_hyperplane_sections ( int m, int k, int n, double *c )
{
   int a[3],*b,fail;
   a[0] = m; a[1] = k; a[2] = n;
   fail = _ada_use_c2phc4c(60,a,b,c,0);
   return fail;
}

int retrieve_hyperplane_sections ( int m, int k, int n, int i, double *c )
{
   int a[3],fail;
   a[0] = m; a[1] = k; a[2] = n;
   fail = _ada_use_c2phc4c(61,a,&i,c,0);
   return fail;
}

int set_target_hyperplane_sections ( int i )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(62,&i,b,c,0);
   return fail;
}

int new_standard_slices ( int k, int n )
{
   int i,j,fail;
   double r[2];

   for(i=0; i<k; i++)
      for(j=0; j<n+1; j++)
      {
         random_complex(&r[0],&r[1]);
         fail = _ada_use_c2phc4c(43,&i,&j,r,0);  /* assign coefficient */
      }

   return fail;
}

int new_dobldobl_slices ( int k, int n )
{
   int i,j,fail;
   double r[4];

   for(i=0; i<k; i++)
      for(j=0; j<n+1; j++)
      {
         random_dobldobl_complex(&r[0],&r[2]);
         fail = _ada_use_c2phc4c(633,&i,&j,r,0);  /* assign coefficient */
      }

   return fail;
}

int new_quaddobl_slices ( int k, int n )
{
   int i,j,fail;
   double r[8];

   for(i=0; i<k; i++)
      for(j=0; j<n+1; j++)
      {
          random_quaddobl_complex(&r[0],&r[3]);
          fail = _ada_use_c2phc4c(663,&i,&j,r,0);  /* assign coefficient */
      }

   return fail;
}

int swap_standard_slices ( void )
{
   double *c;
   int fail,*a,*b;
   fail = _ada_use_c2phc4c(46,a,b,c,0); /* swap start with new slices */
   return fail;
}

int swap_dobldobl_slices ( void )
{
   double *c;
   int fail,*a,*b;
   fail = _ada_use_c2phc4c(636,a,b,c,0); /* swap start with new slices */
   return fail;
}

int swap_quaddobl_slices ( void )
{
   double *c;
   int fail,*a,*b;
   fail = _ada_use_c2phc4c(666,a,b,c,0); /* swap start with new slices */
   return fail;
}

int store_standard_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(51,a,b,c,0);   /* standard sols to permutations */
   return fail;
}

int store_dobldobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(641,a,b,c,0);  /* dobldobl sols to permutations */
   return fail;
}

int store_quaddobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(671,a,b,c,0);  /* quaddobl sols to permutations */
   return fail;
}

int restore_standard_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(48,a,b,c,0);   /* first solutions to container */ 
   return fail;
}

int restore_dobldobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(638,a,b,c,0); /* first dobldobl sols to container */ 
   return fail;
}

int restore_quaddobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(668,a,b,c,0); /* first quaddobl sols to container */ 
   return fail;
}

int retrieve_solutions_on_grid ( int i )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(49,&i,b,c,0);
   return fail;
}

int in_slice ( int label, int slice, int *position )
{
   int a[2],fail;
   double *c;

   a[0] = label;
   a[1] = slice;
   fail = _ada_use_c2phc4c(58,a,position,c,0);

   return fail;
}

int standard_sample_to_new_slices ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(45,a,b,c,0);   /* track paths */
   return fail;
}

int dobldobl_sample_to_new_slices ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(635,a,b,c,0);   /* track paths */
   return fail;
}

int quaddobl_sample_to_new_slices ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(665,a,b,c,0);   /* track paths */
   return fail;
}

int standard_track_paths ( int islaurent )
{
   int fail;

   fail = standard_sample_to_new_slices();            /* do path tracking */
   if(verbose>0) printf("Done tracking.\n");
   if(verbose>1) printf("Solutions computed :\n");
   if(verbose>1) fail = solcon_write_standard_solutions();
   fail = swap_standard_slices();           /* swap start with new slices */

   if(islaurent == 1)
   {
      fail = standard_witness_set_to_Laurent_system_container();
      fail = standard_Newton_Laurent_step();
   }
   else
   {
      fail = standard_witness_set_to_system_container();
      fail = standard_Newton_step();
   }
   return fail;
}

int dobldobl_track_paths ( int islaurent )
{
   int *a,*b,fail;
   double *c;

   fail = dobldobl_sample_to_new_slices();   /* do path tracking */
   if(verbose>0) printf("Done tracking.\n");
   if(verbose>1) printf("Solutions computed :\n");
   if(verbose>1) fail = solcon_write_dobldobl_solutions();
   fail = swap_dobldobl_slices();            /* swap start with new slices */

   if(islaurent == 1) 
   {
      fail = dobldobl_witness_set_to_Laurent_system_container();
      fail = dobldobl_Newton_Laurent_step();
   }
   else
   {
      fail = dobldobl_witness_set_to_system_container();
      fail = dobldobl_Newton_step();
   }
   return fail;
}

int quaddobl_track_paths ( int islaurent )
{
   int *a,*b,fail;
   double *c;

   fail = quaddobl_sample_to_new_slices();   /* do path tracking */
   if(verbose>0) printf("Done tracking.\n");
   if(verbose>1) printf("Solutions computed :\n");
   if(verbose>1) fail = solcon_write_quaddobl_solutions();
   fail = swap_quaddobl_slices();            /* swap start with new slices */

   if(islaurent == 1)
   {
      fail = quaddobl_witness_set_to_Laurent_system_container();
      fail = quaddobl_Newton_Laurent_step();
   }
   else
   {
      fail = quaddobl_witness_set_to_system_container();
      fail = quaddobl_Newton_step();
   }
   return fail;
}

int standard_sample_loop
 ( int start_slice, int target_slice,
   int start_label, int *target_label )
{
   int a[2],fail;
   double *c;

   a[0] = start_slice;
   a[1] = target_slice;
   *target_label = start_label;
   fail = _ada_use_c2phc4c(63,a,target_label,c,0);  /* calls sample loop */

   return fail;
}

int dobldobl_sample_loop
 ( int start_slice, int target_slice,
   int start_label, int *target_label )
{
   int a[2],fail;
   double *c;

   a[0] = start_slice;
   a[1] = target_slice;
   *target_label = start_label;
   fail = _ada_use_c2phc4c(653,a,target_label,c,0);  /* calls sample loop */

   return fail;
}

int quaddobl_sample_loop
 ( int start_slice, int target_slice,
   int start_label, int *target_label )
{
   int a[2],fail;
   double *c;

   a[0] = start_slice;
   a[1] = target_slice;
   *target_label = start_label;
   fail = _ada_use_c2phc4c(683,a,target_label,c,0);  /* calls sample loop */

   return fail;
}

int standard_trace_sum_difference ( int n, int *f, double *d )
{
   int fail;
   fail = _ada_use_c2phc4c(57,&n,f,d,0);
   return fail;
}

int dobldobl_trace_sum_difference ( int n, int *f, double *d )
{
   int fail;
   fail = _ada_use_c2phc4c(647,&n,f,d,0);
   return fail;
}

int quaddobl_trace_sum_difference ( int n, int *f, double *d )
{
   int fail;
   fail = _ada_use_c2phc4c(677,&n,f,d,0);
   return fail;
}

int number_of_standard_factors ( int *nf )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(68,nf,b,c,0);

   return fail;
}

int number_of_dobldobl_factors ( int *nf )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(656,nf,b,c,0);

   return fail;
}

int number_of_quaddobl_factors ( int *nf )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc4c(686,nf,b,c,0);

   return fail;
}

int witness_points_of_standard_factor ( int k, int *d, int *w )
{
   int fail;
   double *c;

   *d = k;
   fail = _ada_use_c2phc4c(69,d,w,c,0);

   return fail;
}

int witness_points_of_dobldobl_factor ( int k, int *d, int *w )
{
   int fail;
   double *c;

   *d = k;
   fail = _ada_use_c2phc4c(657,d,w,c,0);

   return fail;
}

int witness_points_of_quaddobl_factor ( int k, int *d, int *w )
{
   int fail;
   double *c;

   *d = k;
   fail = _ada_use_c2phc4c(687,d,w,c,0);

   return fail;
}

int permutation_after_standard_loop ( int d, int *permutation )
{
   int *a,fail;
   double *c;

   fail = _ada_use_c2phc4c(52,a,permutation,c,0);   /* compute permutation */

   return fail;
}

int permutation_after_dobldobl_loop ( int d, int *permutation )
{
   int *a,fail;
   double *c;

   fail = _ada_use_c2phc4c(642,a,permutation,c,0);   /* compute permutation */

   return fail;
}

int permutation_after_quaddobl_loop ( int d, int *permutation )
{
   int *a,fail;
   double *c;

   fail = _ada_use_c2phc4c(672,a,permutation,c,0);  /* compute permutation */

   return fail;
}

int update_standard_decomposition
 ( int d, int *permutation, int *nf, int *done )
{
   int fail;
   int *b;
   double *c;

   nf[0] = d; nf[1] = 0;
   fail = _ada_use_c2phc4c(53,nf,permutation,c,0);
   fail = _ada_use_c2phc4c(55,done,b,c,0);

   return fail;
}

int update_dobldobl_decomposition
 ( int d, int *permutation, int *nf, int *done )
{
   int fail;
   int *b;
   double *c;

   nf[0] = d; nf[1] = 0;
   fail = _ada_use_c2phc4c(643,nf,permutation,c,0);
   fail = _ada_use_c2phc4c(645,done,b,c,0);

   return fail;
}

int update_quaddobl_decomposition
 ( int d, int *permutation, int *nf, int *done )
{
   int fail;
   int *b;
   double *c;

   nf[0] = d; nf[1] = 0;
   fail = _ada_use_c2phc4c(673,nf,permutation,c,0);
   fail = _ada_use_c2phc4c(675,done,b,c,0);

   return fail;
}

int standard_monodromy_permutation ( int d, int *done )
{
   int *a,*b,fail,i;
   int permutation[d];
   int nf[2];
   double *c;

   fail = _ada_use_c2phc4c(52,a,permutation,c,0);   /* compute permutation */
   if(verbose>0)
   {
      printf("the permutation :");
      for (i=0; i<d; i++) printf(" %d",permutation[i]);
   }
   nf[0] = d;
   nf[1] = 0;
   fail = _ada_use_c2phc4c(53,nf,permutation,c,0); /* update decomposition */
   if(verbose>0) printf(" : %d -> %d\n",nf[0],nf[1]);
   fail = _ada_use_c2phc4c(55,done,b,c,0);        /* apply linear trace */
   /* fail = _ada_use_c2phc4c(54,a,b,c); */       /* write decomposition */

   return fail;
}

int dobldobl_monodromy_permutation ( int d, int *done )
{
   int *a,*b,fail,i;
   int permutation[d];
   int nf[2];
   double *c;

   fail = _ada_use_c2phc4c(642,a,permutation,c,0);  /* compute permutation */
   if(verbose>0)
   {
      printf("the permutation :");
      for (i=0; i<d; i++) printf(" %d",permutation[i]);
   }
   nf[0] = d;
   nf[1] = 0;
   fail = _ada_use_c2phc4c(643,nf,permutation,c,0); /* update decomposition */
   if(verbose>0) printf(" : %d -> %d\n",nf[0],nf[1]);
   fail = _ada_use_c2phc4c(645,done,b,c,0);        /* apply linear trace */
   /* fail = _ada_use_c2phc4c(644,a,b,c); */       /* write decomposition */

   return fail;
}

int quaddobl_monodromy_permutation ( int d, int *done )
{
   int *a,*b,fail,i;
   int permutation[d];
   int nf[2];
   double *c;

   fail = _ada_use_c2phc4c(672,a,permutation,c,0);  /* compute permutation */
   if(verbose>0)
   {
      printf("the permutation :");
      for (i=0; i<d; i++) printf(" %d",permutation[i]);
   }
   nf[0] = d;
   nf[1] = 0;
   fail = _ada_use_c2phc4c(673,nf,permutation,c,0); /* update decomposition */
   if(verbose>0) printf(" : %d -> %d\n",nf[0],nf[1]);
   fail = _ada_use_c2phc4c(675,done,b,c,0);        /* apply linear trace */
   /* fail = _ada_use_c2phc4c(674,a,b,c); */       /* write decomposition */

   return fail;
}

int standard_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[3];
   double cffs[2+2*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   dims[2] = nbtasks;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<2*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc4c(537,&vrb,dims,cffs,0);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int dobldobl_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[3];
   double cffs[2+4*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   dims[2] = nbtasks;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<4*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc4c(538,&vrb,dims,cffs,0);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int quaddobl_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[3];
   double cffs[2+8*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   dims[2] = nbtasks;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<8*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc4c(539,&vrb,dims,cffs,0);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int standard_Laurent_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[3];
   double cffs[2+2*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   dims[2] = nbtasks;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<2*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc4c(795,&vrb,dims,cffs,0);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int dobldobl_Laurent_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[3];
   double cffs[2+4*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   dims[2] = nbtasks;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<4*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc4c(796,&vrb,dims,cffs,0);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int quaddobl_Laurent_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[3];
   double cffs[2+8*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   dims[2] = nbtasks;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<8*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc4c(797,&vrb,dims,cffs,0);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int standard_homotopy_ismember
 ( int vrb, int nvr, int dim, int nbc, char *tpt,
   double restol, double homtol, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[5];
   double tols[2];
   int testsol[nbc];

   dims[0] = vrb;
   dims[1] = nvr;
   dims[2] = dim;
   dims[3] = nbc;
   dims[4] = nbtasks;
   tols[0] = restol;
   tols[1] = homtol;

   for(k=0; k<nbc; k++) testsol[k] = (int) tpt[k];

   fail = _ada_use_c2phc4c(822,dims,testsol,tols,0);

   *onsys = dims[0];
   *onset = testsol[0];

   return fail;
}

int dobldobl_homotopy_ismember
 ( int vrb, int nvr, int dim, int nbc, char *tpt,
   double restol, double homtol, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[5];
   double tols[2];
   int testsol[nbc];

   dims[0] = vrb;
   dims[1] = nvr;
   dims[2] = dim;
   dims[3] = nbc;
   dims[4] = nbtasks;
   tols[0] = restol;
   tols[1] = homtol;

   for(k=0; k<nbc; k++) testsol[k] = (int) tpt[k];

   fail = _ada_use_c2phc4c(823,dims,testsol,tols,0);

   *onsys = dims[0];
   *onset = testsol[0];

   return fail;
}

int quaddobl_homotopy_ismember
 ( int vrb, int nvr, int dim, int nbc, char *tpt,
   double restol, double homtol, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[5];
   double tols[2];
   int testsol[nbc];

   dims[0] = vrb;
   dims[1] = nvr;
   dims[2] = dim;
   dims[3] = nbc;
   dims[4] = nbtasks;
   tols[0] = restol;
   tols[1] = homtol;

   for(k=0; k<nbc; k++) testsol[k] = (int) tpt[k];

   fail = _ada_use_c2phc4c(824,dims,testsol,tols,0);

   *onsys = dims[0];
   *onset = testsol[0];

   return fail;
}

int standard_Laurent_homotopy_ismember
 ( int vrb, int nvr, int dim, int nbc, char *tpt,
   double restol, double homtol, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[5];
   double tols[2];
   int testsol[nbc];

   dims[0] = vrb;
   dims[1] = nvr;
   dims[2] = dim;
   dims[3] = nbc;
   dims[4] = nbtasks;
   tols[0] = restol;
   tols[1] = homtol;

   for(k=0; k<nbc; k++) testsol[k] = (int) tpt[k];

   fail = _ada_use_c2phc4c(825,dims,testsol,tols,0);

   *onsys = dims[0];
   *onset = testsol[0];

   return fail;
}

int dobldobl_Laurent_homotopy_ismember
 ( int vrb, int nvr, int dim, int nbc, char *tpt,
   double restol, double homtol, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[5];
   double tols[2];
   int testsol[nbc];

   dims[0] = vrb;
   dims[1] = nvr;
   dims[2] = dim;
   dims[3] = nbc;
   dims[4] = nbtasks;
   tols[0] = restol;
   tols[1] = homtol;

   for(k=0; k<nbc; k++) testsol[k] = (int) tpt[k];

   fail = _ada_use_c2phc4c(826,dims,testsol,tols,0);

   *onsys = dims[0];
   *onset = testsol[0];

   return fail;
}

int quaddobl_Laurent_homotopy_ismember
 ( int vrb, int nvr, int dim, int nbc, char *tpt,
   double restol, double homtol, int *onsys, int *onset, int nbtasks )
{
   int fail,k,idx;
   int dims[5];
   double tols[2];
   int testsol[nbc];

   dims[0] = vrb;
   dims[1] = nvr;
   dims[2] = dim;
   dims[3] = nbc;
   dims[4] = nbtasks;
   tols[0] = restol;
   tols[1] = homtol;

   for(k=0; k<nbc; k++) testsol[k] = (int) tpt[k];

   fail = _ada_use_c2phc4c(827,dims,testsol,tols,0);

   *onsys = dims[0];
   *onset = testsol[0];

   return fail;
}
