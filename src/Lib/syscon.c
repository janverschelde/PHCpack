/* This file "syscon.c" contains the definitions of the operations
 * declared in the file "syscon.h". */

/* #include<stdio.h> only used for extra print statements */
#include <stdlib.h>
#include "syscon.h"

int syscon_read_standard_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(20,a,b,c);
   return fail;
}

int syscon_read_standard_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(120,a,b,c);
   return fail;
}

int syscon_read_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(330,a,b,c);
   return fail;
}

int syscon_read_dobldobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(550,a,b,c);
   return fail;
}

int syscon_read_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(380,a,b,c);
   return fail;
}

int syscon_read_quaddobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(560,a,b,c);
   return fail;
}

int syscon_read_multprec_system ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(440,&deci,b,c);
   return fail;
}

int syscon_read_multprec_Laurent_system ( int deci )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(570,&deci,b,c);
   return fail;
}

int syscon_random_system ( int n, int m, int d, int c )
{
   int b[3],fail;
   double *tmp;
   b[0] = m;
   b[1] = d;
   b[2] = c;
   fail = _ada_use_c2phc(109,&n,b,tmp);
   return fail;
}

int syscon_write_standard_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(21,a,b,c);
   return fail;
}

int syscon_write_standard_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(121,a,b,c);
   return fail;
}

int syscon_write_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(331,a,b,c);
   return fail;
}

int syscon_write_dobldobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(551,a,b,c);
   return fail;
}

int syscon_write_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(381,a,b,c);
   return fail;
}

int syscon_write_quaddobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(561,a,b,c);
   return fail;
}

int syscon_write_multprec_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(441,a,b,c);
   return fail;
}

int syscon_write_multprec_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(571,a,b,c);
   return fail;
}

int syscon_number_of_standard_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(22,length,b,c);
   return fail;
}

int syscon_number_of_standard_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(122,length,b,c);
   return fail;
}

int syscon_number_of_dobldobl_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(332,length,b,c);
   return fail;
}

int syscon_number_of_dobldobl_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(552,length,b,c);
   return fail;
}

int syscon_number_of_quaddobl_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(382,length,b,c);
   return fail;
}

int syscon_number_of_quaddobl_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(562,length,b,c);
   return fail;
}

int syscon_number_of_multprec_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(442,length,b,c);
   return fail;
}

int syscon_number_of_multprec_Laurentials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(572,length,b,c);
   return fail;
}

int syscon_initialize_number_of_standard_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(23,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_standard_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(123,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_dobldobl_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(333,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_dobldobl_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(553,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_quaddobl_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(383,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_quaddobl_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(563,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_multprec_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(443,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_multprec_Laurentials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(573,&length,b,c);
   return fail;
}

int syscon_degree_of_standard_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(119,&k,d,c);
   return fail;
}

int syscon_degree_of_dobldobl_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(339,&k,d,c);
   return fail;
}

int syscon_degree_of_quaddobl_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(389,&k,d,c);
   return fail;
}

int syscon_degree_of_multprec_polynomial ( int k, int *d )
{
   int fail;
   double *c;
   fail = _ada_use_c2phc(449,&k,d,c);
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

   fail = _ada_use_c2phc(76,a,b,c);

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

   fail = _ada_use_c2phc(338,a,b,c);

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

   fail = _ada_use_c2phc(388,a,b,c);

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

   fail = _ada_use_c2phc(448,a,b,c);

   return fail;
}

int syscon_standard_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(600,&k,szl,c);

   return fail;
}

int syscon_dobldobl_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(601,&k,szl,c);

   return fail;
}

int syscon_quaddobl_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(602,&k,szl,c);

   return fail;
}

int syscon_multprec_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(603,&k,szl,c);

   return fail;
}

int syscon_standard_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(604,&k,szl,c);

   return fail;
}

int syscon_dobldobl_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(605,&k,szl,c);

   return fail;
}

int syscon_quaddobl_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(606,&k,szl,c);

   return fail;
}

int syscon_multprec_Laurent_size_limit ( int k, int *szl )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(607,&k,szl,c);

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
   fail = _ada_use_c2phc(67,&size,buffer,c);
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
   fail = _ada_use_c2phc(106,&size,buffer,c);
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
   fail = _ada_use_c2phc(107,&size,buffer,c);
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
   fail = _ada_use_c2phc(108,&size,buffer,c);
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

   fail = _ada_use_c2phc(74,a,b,c);

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

   fail = _ada_use_c2phc(558,a,b,c);

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

   fail = _ada_use_c2phc(568,a,b,c);

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

   fail = _ada_use_c2phc(578,a,b,c);

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
   fail = _ada_use_c2phc(128,&size,buffer,c);
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
   fail = _ada_use_c2phc(559,&size,buffer,c);
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
   fail = _ada_use_c2phc(569,&size,buffer,c);
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
   fail = _ada_use_c2phc(579,&size,buffer,c);
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
   fail = _ada_use_c2phc(147,a,b,c);
   return fail;
}

int syscon_create_Jacobian_evaluator ( void )
{
   int fail,*a,*b;
   double *c;
   fail = _ada_use_c2phc(148,a,b,c);
   return fail;
}

int syscon_number_of_standard_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(24,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_number_of_standard_Laurent_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(124,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_number_of_dobldobl_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(334,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_number_of_dobldobl_Laurent_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(554,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_number_of_quaddobl_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(384,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_number_of_quaddobl_Laurent_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(564,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_number_of_multprec_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(444,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_retrieve_standard_term ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc(25,a,exp,c);  
   return fail;
}

int syscon_retrieve_dobldobl_term ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc(335,a,exp,c);  
   return fail;
}

int syscon_retrieve_dobldobl_Laurent_term
 ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc(555,a,exp,c);  
   return fail;
}

int syscon_retrieve_quaddobl_term
 ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc(385,a,exp,c);  
   return fail;
}

int syscon_retrieve_quaddobl_Laurent_term
 ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc(565,a,exp,c);  
   return fail;
}

int syscon_add_standard_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc(26,a,exp,c); 
   return fail;
}

int syscon_add_dobldobl_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc(336,a,exp,c); 
   return fail;
}

int syscon_add_dobldobl_Laurent_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc(556,a,exp,c); 
   return fail;
}

int syscon_add_quaddobl_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc(386,a,exp,c); 
   return fail;
}

int syscon_add_quaddobl_Laurent_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc(566,a,exp,c); 
   return fail;
}

int syscon_total_degree ( int *d )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(28,d,b,c);
   return fail;
}

int syscon_clear_standard_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(27,a,b,c);
   return fail;
}

int syscon_clear_standard_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(127,a,b,c);
   return fail;
}

int syscon_clear_dobldobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(337,a,b,c);
   return fail;
}

int syscon_clear_dobldobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(557,a,b,c);
   return fail;
}

int syscon_clear_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(387,a,b,c);
   return fail;
}

int syscon_clear_quaddobl_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(567,a,b,c);
   return fail;
}

int syscon_clear_multprec_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(447,a,b,c);
   return fail;
}

int syscon_clear_multprec_Laurent_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(577,a,b,c);
   return fail;
}

int syscon_number_of_symbols ( int *n )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(293,n,b,c);
   return fail;
}

int syscon_write_symbols ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(294,a,b,c);
   return fail;
}

int syscon_string_of_symbols ( int *n, char *s )
{
   int i,b[*n],fail;
   double *c;

   fail = _ada_use_c2phc(295,n,b,c);

   for(i=0; i<*n; i++)
      s[i] = (char)b[i];
   s[*n] = '\0';

   return fail;
}

int syscon_clear_symbol_table ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(29,a,b,c);
   return fail;
}

int syscon_remove_symbol_from_table ( int i )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(291,&i,b,c);
   return fail;
}

int syscon_remove_symbol_name_from_table ( int nc, char *s )
{
   int i,fail;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc(296,&nc,b,c);
   return fail;
}

int syscon_sort_embed_symbols ( int *nzz )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(292,nzz,b,c);
   return fail;
}

int syscon_standard_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc(306,&k,b,c);
   return fail;
}

int syscon_standard_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc(309,&nc,b,c);
   return fail;
}

int syscon_dobldobl_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc(307,&k,b,c);
   return fail;
}

int syscon_dobldobl_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc(310,&nc,b,c);
   return fail;
}

int syscon_quaddobl_drop_variable_by_index ( int k )
{ 
   int fail,*b;
   double *c;
   fail = _ada_use_c2phc(308,&k,b,c);
   return fail;
}

int syscon_quaddobl_drop_variable_by_name ( int nc, char *s )
{ 
   int fail,i;
   int b[nc];
   double *c;
   for(i=0; i<nc; i++) b[i] = (int) s[i];
   fail = _ada_use_c2phc(311,&nc,b,c);
   return fail;
}
