/* This file "syscon.c" contains the definitions of the operations
 * declared in the file "syscon.h". */

/* #include<stdio.h> only used for extra print statements */
#include <stdlib.h>
#include "syscon.h"

int syscon_read_standard_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(20,a,b,c,0);
   return fail;
}

int syscon_read_standard_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(120,a,b,c,0);
   return fail;
}

int syscon_read_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(330,a,b,c,0);
   return fail;
}

int syscon_read_dobldobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(550,a,b,c,0);
   return fail;
}

int syscon_read_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(380,a,b,c,0);
   return fail;
}

int syscon_read_quaddobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(560,a,b,c,0);
   return fail;
}

int syscon_read_multprec_system ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(440,&deci,b,c,0);
   return fail;
}

int syscon_read_multprec_Laurent_system ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(570,&deci,b,c,0);
   return fail;
}

int syscon_random_system ( int n, int m, int d, int c, int neq )
{
   int a[2],b[3],fail;
   double *tmp;

   a[0] = n;
   a[1] = neq;
   b[0] = m;
   b[1] = d;
   b[2] = c;

   fail = _ada_use_c2phc4c(109,a,b,tmp,0);

   return fail;
}

int syscon_dobldobl_random_system ( int n, int m, int d, int c, int neq )
{
   int a[2],b[3],fail;
   double *tmp;

   a[0] = n;
   a[1] = neq;
   b[0] = m;
   b[1] = d;
   b[2] = c;

   fail = _ada_use_c2phc4c(548,a,b,tmp,0);

   return fail;
}

int syscon_quaddobl_random_system ( int n, int m, int d, int c, int neq )
{
   int a[2],b[3],fail;
   double *tmp;

   a[0] = n;
   a[1] = neq;
   b[0] = m;
   b[1] = d;
   b[2] = c;

   fail = _ada_use_c2phc4c(549,a,b,tmp,0);

   return fail;
}

int syscon_write_standard_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(21,a,b,c,0);
   return fail;
}

int syscon_write_standard_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(121,a,b,c,0);
   return fail;
}

int syscon_write_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(331,a,b,c,0);
   return fail;
}

int syscon_write_dobldobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(551,a,b,c,0);
   return fail;
}

int syscon_write_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(381,a,b,c,0);
   return fail;
}

int syscon_write_quaddobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(561,a,b,c,0);
   return fail;
}

int syscon_write_multprec_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(441,a,b,c,0);
   return fail;
}

int syscon_write_multprec_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(571,a,b,c,0);
   return fail;
}

int syscon_number_of_standard_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(22,length,b,c,0);
   return fail;
}

int syscon_number_of_standard_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(122,length,b,c,0);
   return fail;
}

int syscon_number_of_dobldobl_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(332,length,b,c,0);
   return fail;
}

int syscon_number_of_dobldobl_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(552,length,b,c,0);
   return fail;
}

int syscon_number_of_quaddobl_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(382,length,b,c,0);
   return fail;
}

int syscon_number_of_quaddobl_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(562,length,b,c,0);
   return fail;
}

int syscon_number_of_multprec_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(442,length,b,c,0);
   return fail;
}

int syscon_number_of_multprec_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(572,length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_standard_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(23,&length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_standard_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(123,&length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_dobldobl_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(333,&length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_dobldobl_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(553,&length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_quaddobl_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(383,&length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_quaddobl_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(563,&length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_multprec_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(443,&length,b,c,0);
   return fail;
}

int syscon_initialize_number_of_multprec_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(573,&length,b,c,0);
   return fail;
}

int syscon_degree_of_standard_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(119,&k,d,c,0);
   return fail;
}

int syscon_degree_of_dobldobl_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(339,&k,d,c,0);
   return fail;
}

int syscon_degree_of_quaddobl_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(389,&k,d,c,0);
   return fail;
}

int syscon_degree_of_multprec_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc4c(449,&k,d,c,0);
   return fail;
}

int syscon_store_standard_polynomial ( int nc, int n, int k, char *p )
{
   int a[3],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(76,a,b,c,0);

   /* if(fail != 0) printf("Failed to store a polynomial.\n"); */

   return fail;
}

int syscon_store_dobldobl_polynomial ( int nc, int n, int k, char *p )
{
   int a[3],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(338,a,b,c,0);

   return fail;
}

int syscon_store_quaddobl_polynomial ( int nc, int n, int k, char *p )
{
   int a[3],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(388,a,b,c,0);

   return fail;
}

int syscon_store_multprec_polynomial
 ( int nc, int n, int k, int deci, char *p )
{
   int a[4],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   a[3] = deci;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(448,a,b,c,0);

   return fail;
}

int syscon_standard_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(600,&k,szl,c,0);

   return fail;
}

int syscon_dobldobl_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(601,&k,szl,c,0);

   return fail;
}

int syscon_quaddobl_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(602,&k,szl,c,0);

   return fail;
}

int syscon_multprec_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(603,&k,szl,c,0);

   return fail;
}

int syscon_standard_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(604,&k,szl,c,0);

   return fail;
}

int syscon_dobldobl_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(605,&k,szl,c,0);

   return fail;
}

int syscon_quaddobl_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(606,&k,szl,c,0);

   return fail;
}

int syscon_multprec_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc4c(607,&k,szl,c,0);

   return fail;
}

int syscon_load_standard_polynomial ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_standard_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(67,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_load_dobldobl_polynomial ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_dobldobl_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(106,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_load_quaddobl_polynomial ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_quaddobl_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(107,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_load_multprec_polynomial ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_multprec_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(108,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_store_standard_Laurential ( int nc, int n, int k, char *p )
{
   int a[3],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(74,a,b,c,0);

   return fail;
}

int syscon_store_dobldobl_Laurential ( int nc, int n, int k, char *p )
{
   int a[3],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(558,a,b,c,0);

   return fail;
}

int syscon_store_quaddobl_Laurential ( int nc, int n, int k, char *p )
{
   int a[3],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(568,a,b,c,0);

   return fail;
}

int syscon_store_multprec_Laurential
 ( int nc, int n, int k, int deci, char *p )
{
   int a[4],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   a[3] = deci;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc4c(578,a,b,c,0);

   return fail;
}

int syscon_load_standard_Laurential ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_standard_Laurent_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(128,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_load_dobldobl_Laurential ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_dobldobl_Laurent_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(559,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_load_quaddobl_Laurential ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_quaddobl_Laurent_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(569,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_load_multprec_Laurential ( int k, int *nc, char *p )
{
   int fail,i,szl;
   int *buffer;
   int size = k;
   double *c;

   fail = syscon_multprec_Laurent_size_limit(k,&szl);
   buffer = (int*)calloc(szl,sizeof(int));
   fail = _ada_use_c2phc4c(579,&size,buffer,c,0);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) buffer[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */
   free(buffer);

   return fail;
}

int syscon_create_evaluator ( void )
{
   int fail,*a,*b;
   double *c;
   fail = _ada_use_c2phc4c(147,a,b,c,0);
   return fail;
}

int syscon_create_Jacobian_evaluator ( void )
{
   int fail,*a,*b;
   double *c;
   fail = _ada_use_c2phc4c(148,a,b,c,0);
   return fail;
}

int syscon_number_of_standard_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc4c(24,a,b,c,0);
   *nt = a[0];
   return fail;
}

int syscon_number_of_standard_Laurent_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc4c(124,a,b,c,0);
   *nt = a[0];
   return fail;
}

int syscon_number_of_dobldobl_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc4c(334,a,b,c,0);
   *nt = a[0];
   return fail;
}

int syscon_number_of_dobldobl_Laurent_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc4c(554,a,b,c,0);
   *nt = a[0];
   return fail;
}

int syscon_number_of_quaddobl_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc4c(384,a,b,c,0);
   *nt = a[0];
   return fail;
}

int syscon_number_of_quaddobl_Laurent_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc4c(564,a,b,c,0);
   *nt = a[0];
   return fail;
}

int syscon_number_of_multprec_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc4c(444,a,b,c,0);
   *nt = a[0];
   return fail;
}

int syscon_retrieve_standard_term ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc4c(25,a,exp,c,0);  
   return fail;
}

int syscon_retrieve_dobldobl_term ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc4c(335,a,exp,c,0);  
   return fail;
}

int syscon_retrieve_dobldobl_Laurent_term
 ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc4c(555,a,exp,c,0);  
   return fail;
}

int syscon_retrieve_quaddobl_term
 ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc4c(385,a,exp,c,0);  
   return fail;
}

int syscon_retrieve_quaddobl_Laurent_term
 ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc4c(565,a,exp,c,0);  
   return fail;
}

int syscon_add_standard_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc4c(26,a,exp,c,0); 
   return fail;
}

int syscon_add_dobldobl_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc4c(336,a,exp,c,0); 
   return fail;
}

int syscon_add_dobldobl_Laurent_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc4c(556,a,exp,c,0); 
   return fail;
}

int syscon_add_quaddobl_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc4c(386,a,exp,c,0); 
   return fail;
}

int syscon_add_quaddobl_Laurent_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc4c(566,a,exp,c,0); 
   return fail;
}

int syscon_total_degree ( int *d )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(28,d,b,c,0);
   return fail;
}

int syscon_clear_standard_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(27,a,b,c,0);
   return fail;
}

int syscon_clear_standard_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(127,a,b,c,0);
   return fail;
}

int syscon_clear_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(337,a,b,c,0);
   return fail;
}

int syscon_clear_dobldobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(557,a,b,c,0);
   return fail;
}

int syscon_clear_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(387,a,b,c,0);
   return fail;
}

int syscon_clear_quaddobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(567,a,b,c,0);
   return fail;
}

int syscon_clear_multprec_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(447,a,b,c,0);
   return fail;
}

int syscon_clear_multprec_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(577,a,b,c,0);
   return fail;
}

int syscon_number_of_symbols ( int *n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(293,n,b,c,0);
   return fail;
}

int syscon_write_symbols ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(294,a,b,c,0);
   return fail;
}

int syscon_string_of_symbols ( int *n, char *s )
{
   int i,b[*n],fail;
   double *c;

   fail = _ada_use_c2phc4c(295,n,b,c,0);

   for(i=0; i<*n; i++)
      s[i] = (char)b[i];
   s[*n] = '\0';

   return fail;
}

int syscon_clear_symbol_table ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc4c(29,a,b,c,0);
   return fail;
}

int syscon_remove_symbol_from_table ( int i )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(291,&i,b,c,0);
   return fail;
}

int syscon_remove_symbol_name_from_table ( int nc, char *s )
{
   int i,fail;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(296,&nc,b,c,0);
   return fail;
}

int syscon_sort_embed_symbols ( int *nzz )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc4c(292,nzz,b,c,0);
   return fail;
}

int syscon_standard_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(306,&k,b,c,0);
   return fail;
}

int syscon_standard_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(309,&nc,b,c,0);
   return fail;
}

int syscon_dobldobl_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(307,&k,b,c,0);
   return fail;
}

int syscon_dobldobl_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(310,&nc,b,c,0);
   return fail;
}

int syscon_quaddobl_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(308,&k,b,c,0);
   return fail;
}

int syscon_quaddobl_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(311,&nc,b,c,0);
   return fail;
}

int syscon_standard_Laurent_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(828,&k,b,c,0);
   return fail;
}

int syscon_standard_Laurent_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(831,&nc,b,c,0);
   return fail;
}

int syscon_dobldobl_Laurent_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(829,&k,b,c,0);
   return fail;
}

int syscon_dobldobl_Laurent_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(832,&nc,b,c,0);
   return fail;
}

int syscon_quaddobl_Laurent_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc4c(830,&k,b,c,0);
   return fail;
}

int syscon_quaddobl_Laurent_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc4c(833,&nc,b,c,0);
   return fail;
}

int syscon_standard_one_homogenization ( int lintype )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(891,&lintype,b,c,0);

   return fail;
}

int syscon_dobldobl_one_homogenization ( int lintype )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(892,&lintype,b,c,0);

   return fail;
}

int syscon_quaddobl_one_homogenization ( int lintype )
{
   int fail;
   int *b;
   double *c;

   fail = _ada_use_c2phc4c(893,&lintype,b,c,0);

   return fail;
}

int syscon_standard_multi_homogenization
 ( int n, int m, int *idz, int lintype )
{
   int fail;
   int pars[3];
   double *c;

   pars[0] = n;
   pars[1] = m;
   pars[2] = lintype;

   fail = _ada_use_c2phc4c(904,pars,idz,c,0);

   return fail;
}

int syscon_dobldobl_multi_homogenization
 ( int n, int m, int *idz, int lintype )
{
   int fail;
   int pars[3];
   double *c;

   pars[0] = n;
   pars[1] = m;
   pars[2] = lintype;

   fail = _ada_use_c2phc4c(905,pars,idz,c,0);

   return fail;
}

int syscon_quaddobl_multi_homogenization
 ( int n, int m, int *idz, int lintype )
{
   int fail;
   int pars[3];
   double *c;

   pars[0] = n;
   pars[1] = m;
   pars[2] = lintype;

   fail = _ada_use_c2phc4c(906,pars,idz,c,0);

   return fail;
}

int syscon_add_symbol ( int nbc, char *name )
{
   int fail,idx;
   int chrs[nbc];
   double *c;

   for(idx=0; idx<nbc; idx++) chrs[idx] = (int) name[idx];

   fail = _ada_use_c2phc4c(897,&nbc,chrs,c,0);

   return fail;
}

int syscon_standard_one_affinization ( void )
{
   int fail,a,b;
   double c;

   fail = _ada_use_c2phc4c(901,&a,&b,&c,0);

   return fail;
}

int syscon_dobldobl_one_affinization ( void )
{
   int fail,a,b;
   double c;

   fail = _ada_use_c2phc4c(902,&a,&b,&c,0);

   return fail;
}

int syscon_quaddobl_one_affinization ( void )
{
   int fail,a,b;
   double c;

   fail = _ada_use_c2phc4c(903,&a,&b,&c,0);

   return fail;
}

int syscon_standard_multi_affinization ( int m )
{
   int fail,b;
   double c;

   fail = _ada_use_c2phc4c(907,&m,&b,&c,0);

   return fail;
}

int syscon_dobldobl_multi_affinization ( int m )
{
   int fail,b;
   double c;

   fail = _ada_use_c2phc4c(908,&m,&b,&c,0);

   return fail;
}

int syscon_quaddobl_multi_affinization ( int m )
{
   int fail,b;
   double c;

   fail = _ada_use_c2phc4c(909,&m,&b,&c,0);

   return fail;
}
