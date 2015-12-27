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

int embed_system ( int d, int precision )
{
   int fail;

   if(precision == 0) fail = embed_standard_system(d);
   if(precision == 1) fail = embed_dobldobl_system(d);
   if(precision == 2) fail = embed_quaddobl_system(d);

   return fail;
}

int embed_standard_system ( int d )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc(66,&d,b,c);
   return fail;
}

int embed_dobldobl_system ( int d )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc(129,&d,b,c);
   return fail;
}

int embed_quaddobl_system ( int d )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc(260,&d,b,c);
   return fail;
}

int read_witness_set ( int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

   fail = _ada_use_c2phc(41,n,d,c);    /* read the witness set */
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

   fail = _ada_use_c2phc(631,n,d,c);    /* read the witness set */
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

   fail = _ada_use_c2phc(661,n,d,c);    /* read the witness set */
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
   fail = _ada_use_c2phc(64,n,b,c);
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
   fail = _ada_use_c2phc(654,n,b,c);
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
   fail = _ada_use_c2phc(684,n,b,c);
   *dim = b[0];
   *deg = b[1];

   return fail;
}

int write_witness_set_to_file ( int m, char *s )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(65,&m,b,c);

   return fail;
}

int write_dobldobl_witness_set_to_file ( int m, char *s )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(655,&m,b,c);

   return fail;
}

int write_quaddobl_witness_set_to_file ( int m, char *s )
{
   int b[m],i,fail;
   double *c;

   for (i=0; i<m; i++) b[i] = (int) s[i];

   fail = _ada_use_c2phc(685,&m,b,c);

   return fail;
}

int read_a_witness_set ( int k, int *n, int *dim, int *deg )
{
   int fail;
   double *c;
   int d[2];

   *n = k;
   fail = _ada_use_c2phc(166,n,d,c);
   *dim = d[0];
   *deg = d[1];

   return fail;
}

int witness_set_to_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(47,a,b,c);   /* copy system to container */
   return fail;
}

int dobldobl_witness_set_to_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(637,a,b,c);   /* copy system to container */
   return fail;
}

int quaddobl_witness_set_to_system_container ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(667,a,b,c);   /* copy system to container */
   return fail;
}

int create_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(164,a,b,c);
   return fail;
}

int create_dobldobl_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(178,a,b,c);
   return fail;
}

int create_quaddobl_cascade_homotopy ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(188,a,b,c);
   return fail;
}

/* OPERATIONS to intersect witness sets */

int standard_diagonal_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(165,&a,&b,c);
   return fail;
}

int dobldobl_diagonal_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(289,&a,&b,c);
   return fail;
}

int quaddobl_diagonal_homotopy ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(290,&a,&b,c);
   return fail;
}

int standard_diagonal_cascade_solutions ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(271,&a,&b,c);
   return fail;
}

int dobldobl_diagonal_cascade_solutions ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(297,&a,&b,c);
   return fail;
}

int quaddobl_diagonal_cascade_solutions ( int a, int b )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(298,&a,&b,c);
   return fail;
}

int extrinsic_top_diagonal_dimension
 ( int n1, int n2, int a, int b, int *d )
{
   int fail,alpha[2],beta[2];
   double *c;
   alpha[0] = n1; alpha[1] = n2;
   beta[0] = a; beta[1] = b;

   fail = _ada_use_c2phc(168,alpha,beta,c);

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

   fail = _ada_use_c2phc(169,a,b,c);
   
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

   fail = _ada_use_c2phc(270,a,b,c);

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

   fail = _ada_use_c2phc(259,a,b,c);

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

   fail = _ada_use_c2phc(269,a,b,c);

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
 
   fail = _ada_use_c2phc(230,a,b,c);

   return fail;
}

int standard_collapse_diagonal ( int k, int d )
{
   int a[2],*b,fail;
   double *c;

   a[0] = k; a[1] = d;
   fail = _ada_use_c2phc(170,a,b,c);

   return fail;
}

int dobldobl_collapse_diagonal ( int k, int d )
{
   int a[2],*b,fail;
   double *c;

   a[0] = k; a[1] = d;
   fail = _ada_use_c2phc(299,a,b,c);

   return fail;
}

int quaddobl_collapse_diagonal ( int k, int d )
{
   int a[2],*b,fail;
   double *c;

   a[0] = k; a[1] = d;
   fail = _ada_use_c2phc(312,a,b,c);

   return fail;
}

int remove_last_slack ( int k )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(171,&k,b,c);

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

int set_state_to_silent ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(39,a,b,c);

   return fail;
}

int set_dobldobl_state_to_silent ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(658,a,b,c);

   return fail;
}

int set_quaddobl_state_to_silent ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(688,a,b,c);

   return fail;
}

int set_standard_state_to_verbose ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(630,a,b,c);

   return fail;
}

int set_dobldobl_state_to_verbose ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(660,a,b,c);

   return fail;
}

int set_quaddobl_state_to_verbose ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(690,a,b,c);

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

int initialize_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(42,&dim,b,c);        /* initialize sampler */

   return fail;
}

int initialize_dobldobl_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(632,&dim,b,c);        /* initialize sampler */

   return fail;
}

int initialize_quaddobl_sampler ( int dim )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(662,&dim,b,c);        /* initialize sampler */

   return fail;
}

int initialize_monodromy ( int n, int d, int k )
{
   int fail;
   double *c;
   int dk[2];

   dk[0] = d;
   dk[1] = k;
   fail = _ada_use_c2phc(50,&n,dk,c); /* initialize Monodromy_Permutations */

   return fail;
}

int initialize_dobldobl_monodromy ( int n, int d, int k )
{
   int fail;
   double *c;
   int dk[2];

   dk[0] = d;
   dk[1] = k;
   fail = _ada_use_c2phc(640,&n,dk,c);

   return fail;
}

int initialize_quaddobl_monodromy ( int n, int d, int k )
{
   int fail;
   double *c;
   int dk[2];

   dk[0] = d;
   dk[1] = k;
   fail = _ada_use_c2phc(670,&n,dk,c);

   return fail;
}

int trace_grid_diagnostics ( double *err, double *dis )
{
   int *a,*b,fail;
   double c[2];

   fail = _ada_use_c2phc(56,a,b,c);
   *err = c[0];
   *dis = c[1];

   return fail;
}

int dobldobl_trace_grid_diagnostics ( double *err, double *dis )
{
   int *a,*b,fail;
   double c[2];

   fail = _ada_use_c2phc(646,a,b,c);
   *err = c[0];
   *dis = c[1];

   return fail;
}

int quaddobl_trace_grid_diagnostics ( double *err, double *dis )
{
   int *a,*b,fail;
   double c[2];

   fail = _ada_use_c2phc(676,a,b,c);
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
   
   fail = _ada_use_c2phc(280,a,b,c);
   *re = c[0];
   *im = c[1];

   return fail;
}

int random_dobldobl_complex ( double *re, double *im )
{
   int *a,*b,fail;
   double c[4];
   
   fail = _ada_use_c2phc(659,a,b,c);
   re[0] = c[0]; re[1] = c[1];
   im[0] = c[2]; im[1] = c[3];

   return fail;
}

int random_quaddobl_complex ( double *re, double *im )
{
   int *a,*b,fail;
   double c[8];
   
   fail = _ada_use_c2phc(689,a,b,c);
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
      fail = _ada_use_c2phc(44,&i,b,r);  /* store gamma */
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
      fail = _ada_use_c2phc(634,&i,b,r);  /* store gamma */
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
      fail = _ada_use_c2phc(664,&i,b,r);  /* store gamma */
   }

   return fail;
}

int assign_coefficient_of_slice ( int i, int j, double *r )
{
   int fail;
   fail = _ada_use_c2phc(43,&i,&j,r);
   return fail;
}

int assign_dobldobl_coefficient_of_slice ( int i, int j, double *r )
{
   int fail;
   fail = _ada_use_c2phc(633,&i,&j,r);
   return fail;
}

int assign_quaddobl_coefficient_of_slice ( int i, int j, double *r )
{
   int fail;
   fail = _ada_use_c2phc(663,&i,&j,r);
   return fail;
}

int initialize_hyperplane_sections ( int m )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(59,&m,b,c);
   return fail;
}

int store_new_hyperplane_sections ( int m, int k, int n, double *c )
{
   int a[3],*b,fail;
   a[0] = m; a[1] = k; a[2] = n;
   fail = _ada_use_c2phc(60,a,b,c);
   return fail;
}

int retrieve_hyperplane_sections ( int m, int k, int n, int i, double *c )
{
   int a[3],fail;
   a[0] = m; a[1] = k; a[2] = n;
   fail = _ada_use_c2phc(61,a,&i,c);
   return fail;
}

int set_target_hyperplane_sections ( int i )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(62,&i,b,c);
   return fail;
}

int new_slices ( int k, int n )
{
   int i,j,fail;
   double r[2];

   for(i=0; i<k; i++)
      for(j=0; j<n+1; j++)
      {
          random_complex(&r[0],&r[1]);
          fail = _ada_use_c2phc(43,&i,&j,r);  /* assign coefficient */
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
          fail = _ada_use_c2phc(633,&i,&j,r);  /* assign coefficient */
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
          fail = _ada_use_c2phc(663,&i,&j,r);  /* assign coefficient */
      }

   return fail;
}

int swap_slices ( void )
{
   double *c;
   int fail,*a,*b;
   fail = _ada_use_c2phc(46,a,b,c); /* swap start with new slices */
   return fail;
}

int swap_dobldobl_slices ( void )
{
   double *c;
   int fail,*a,*b;
   fail = _ada_use_c2phc(636,a,b,c); /* swap start with new slices */
   return fail;
}

int swap_quaddobl_slices ( void )
{
   double *c;
   int fail,*a,*b;
   fail = _ada_use_c2phc(666,a,b,c); /* swap start with new slices */
   return fail;
}

int store_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(51,a,b,c);   /* standard sols to permutations */
   return fail;
}

int store_dobldobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(641,a,b,c);  /* dobldobl sols to permutations */
   return fail;
}

int store_quaddobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(671,a,b,c);  /* quaddobl sols to permutations */
   return fail;
}

int restore_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(48,a,b,c);   /* first solutions to container */ 
   return fail;
}

int restore_dobldobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(638,a,b,c);  /* first dobldobl sols to container */ 
   return fail;
}

int restore_quaddobl_solutions ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(668,a,b,c);  /* first quaddobl sols to container */ 
   return fail;
}

int retrieve_solutions_on_grid ( int i )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(49,&i,b,c);
   return fail;
}

int in_slice ( int label, int slice, int *position )
{
   int a[2],fail;
   double *c;

   a[0] = label;
   a[1] = slice;
   fail = _ada_use_c2phc(58,a,position,c);

   return fail;
}

int sample_to_new_slices ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(45,a,b,c);   /* track paths */
   return fail;
}

int dobldobl_sample_to_new_slices ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(635,a,b,c);   /* track paths */
   return fail;
}

int quaddobl_sample_to_new_slices ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(665,a,b,c);   /* track paths */
   return fail;
}

int track_paths ( void )
{
   int fail;

   fail = sample_to_new_slices();            /* do path tracking */
   if(verbose>0) printf("Done tracking.\n");
   if(verbose>1) printf("Solutions computed :\n");
   if(verbose>1) fail = solcon_write_standard_solutions();
   fail = swap_slices();                     /* swap start with new slices */
   fail = witness_set_to_system_container();
   fail = validate_solutions();

   return fail;
}

int dobldobl_track_paths ( void )
{
   int *a,*b,fail;
   double *c;

   fail = dobldobl_sample_to_new_slices();   /* do path tracking */
   if(verbose>0) printf("Done tracking.\n");
   if(verbose>1) printf("Solutions computed :\n");
   if(verbose>1) fail = solcon_write_dobldobl_solutions();
   fail = swap_dobldobl_slices();            /* swap start with new slices */
   fail = dobldobl_witness_set_to_system_container();
   fail = dobldobl_Newton_step();

   return fail;
}

int quaddobl_track_paths ( void )
{
   int *a,*b,fail;
   double *c;

   fail = quaddobl_sample_to_new_slices();   /* do path tracking */
   if(verbose>0) printf("Done tracking.\n");
   if(verbose>1) printf("Solutions computed :\n");
   if(verbose>1) fail = solcon_write_quaddobl_solutions();
   fail = swap_quaddobl_slices();            /* swap start with new slices */
   fail = quaddobl_witness_set_to_system_container();
   fail = quaddobl_Newton_step();

   return fail;
}

int sample_loop ( int start_slice, int target_slice,
                  int start_label, int *target_label )
{
   int a[2],fail;
   double *c;

   a[0] = start_slice;
   a[1] = target_slice;
   *target_label = start_label;
   fail = _ada_use_c2phc(63,a,target_label,c);  /* calls sample loop */

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
   fail = _ada_use_c2phc(653,a,target_label,c);  /* calls sample loop */

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
   fail = _ada_use_c2phc(683,a,target_label,c);  /* calls sample loop */

   return fail;
}

int trace_sum_difference ( int n, int *f, double *d )
{
   int fail;
   fail = _ada_use_c2phc(57,&n,f,d);
   return fail;
}

int dobldobl_trace_sum_difference ( int n, int *f, double *d )
{
   int fail;
   fail = _ada_use_c2phc(647,&n,f,d);
   return fail;
}

int quaddobl_trace_sum_difference ( int n, int *f, double *d )
{
   int fail;
   fail = _ada_use_c2phc(677,&n,f,d);
   return fail;
}

int number_of_irreducible_factors ( int *nf )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(68,nf,b,c);

   return fail;
}

int number_of_dobldobl_factors ( int *nf )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(656,nf,b,c);

   return fail;
}

int number_of_quaddobl_factors ( int *nf )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(686,nf,b,c);

   return fail;
}

int witness_points_of_irreducible_factor ( int k, int *d, int *w )
{
   int fail;
   double *c;

   *d = k;
   fail = _ada_use_c2phc(69,d,w,c);

   return fail;
}

int witness_points_of_dobldobl_factor ( int k, int *d, int *w )
{
   int fail;
   double *c;

   *d = k;
   fail = _ada_use_c2phc(657,d,w,c);

   return fail;
}

int witness_points_of_quaddobl_factor ( int k, int *d, int *w )
{
   int fail;
   double *c;

   *d = k;
   fail = _ada_use_c2phc(687,d,w,c);

   return fail;
}

int permutation_after_loop ( int d, int *permutation )
{
   int *a,fail;
   double *c;

   fail = _ada_use_c2phc(52,a,permutation,c);   /* compute permutation */

   return fail;
}

int permutation_after_dobldobl_loop ( int d, int *permutation )
{
   int *a,fail;
   double *c;

   fail = _ada_use_c2phc(642,a,permutation,c);   /* compute permutation */

   return fail;
}

int permutation_after_quaddobl_loop ( int d, int *permutation )
{
   int *a,fail;
   double *c;

   fail = _ada_use_c2phc(672,a,permutation,c);   /* compute permutation */

   return fail;
}

int update_decomposition ( int d, int *permutation, int *nf, int *done )
{
   int fail;
   int *b;
   double *c;

   nf[0] = d; nf[1] = 0;
   fail = _ada_use_c2phc(53,nf,permutation,c);
   fail = _ada_use_c2phc(55,done,b,c);

   return fail;
}

int update_dobldobl_decomposition
 ( int d, int *permutation, int *nf, int *done )
{
   int fail;
   int *b;
   double *c;

   nf[0] = d; nf[1] = 0;
   fail = _ada_use_c2phc(643,nf,permutation,c);
   fail = _ada_use_c2phc(645,done,b,c);

   return fail;
}

int update_quaddobl_decomposition
 ( int d, int *permutation, int *nf, int *done )
{
   int fail;
   int *b;
   double *c;

   nf[0] = d; nf[1] = 0;
   fail = _ada_use_c2phc(673,nf,permutation,c);
   fail = _ada_use_c2phc(675,done,b,c);

   return fail;
}

int monodromy_permutation ( int d, int *done )
{
   int *a,*b,fail,i;
   int permutation[d];
   int nf[2];
   double *c;

   fail = _ada_use_c2phc(52,a,permutation,c);   /* compute permutation */
   if(verbose>0)
   {
      printf("the permutation :");
      for (i=0; i<d; i++) printf(" %d",permutation[i]);
   }
   nf[0] = d;
   nf[1] = 0;
   fail = _ada_use_c2phc(53,nf,permutation,c);  /* update decomposition */
   if(verbose>0) printf(" : %d -> %d\n",nf[0],nf[1]);
   fail = _ada_use_c2phc(55,done,b,c);          /* apply linear trace */
   /* fail = _ada_use_c2phc(54,a,b,c); */       /* write decomposition */

   return fail;
}

int dobldobl_monodromy_permutation ( int d, int *done )
{
   int *a,*b,fail,i;
   int permutation[d];
   int nf[2];
   double *c;

   fail = _ada_use_c2phc(642,a,permutation,c);   /* compute permutation */
   if(verbose>0)
   {
      printf("the permutation :");
      for (i=0; i<d; i++) printf(" %d",permutation[i]);
   }
   nf[0] = d;
   nf[1] = 0;
   fail = _ada_use_c2phc(643,nf,permutation,c);  /* update decomposition */
   if(verbose>0) printf(" : %d -> %d\n",nf[0],nf[1]);
   fail = _ada_use_c2phc(645,done,b,c);          /* apply linear trace */
   /* fail = _ada_use_c2phc(644,a,b,c); */       /* write decomposition */

   return fail;
}

int quaddobl_monodromy_permutation ( int d, int *done )
{
   int *a,*b,fail,i;
   int permutation[d];
   int nf[2];
   double *c;

   fail = _ada_use_c2phc(672,a,permutation,c);   /* compute permutation */
   if(verbose>0)
   {
      printf("the permutation :");
      for (i=0; i<d; i++) printf(" %d",permutation[i]);
   }
   nf[0] = d;
   nf[1] = 0;
   fail = _ada_use_c2phc(673,nf,permutation,c);  /* update decomposition */
   if(verbose>0) printf(" : %d -> %d\n",nf[0],nf[1]);
   fail = _ada_use_c2phc(675,done,b,c);          /* apply linear trace */
   /* fail = _ada_use_c2phc(674,a,b,c); */       /* write decomposition */

   return fail;
}

int standard_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset )
{
   int fail,k,idx;
   int dims[2];
   double cffs[2+2*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<2*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc(537,&vrb,dims,cffs);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int dobldobl_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset )
{
   int fail,k,idx;
   int dims[2];
   double cffs[2+4*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<4*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc(538,&vrb,dims,cffs);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}

int quaddobl_homotopy_membership_test
 ( int vrb, int nvr, int dim, double restol, double homtol,
   double *tpt, int *onsys, int *onset )
{
   int fail,k,idx;
   int dims[2];
   double cffs[2+8*nvr];

   dims[0] = nvr;
   dims[1] = dim;
   cffs[0] = restol;
   cffs[1] = homtol;
   idx = 2;
   for(k=0; k<8*nvr; k++) cffs[idx++] = tpt[k];

   fail = _ada_use_c2phc(539,&vrb,dims,cffs);

   *onsys = vrb;
   *onset = dims[0];

   return fail;
}
