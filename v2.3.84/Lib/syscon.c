/* This file "syscon.c" contains the definitions of the operations
 * declared in the file "syscon.h". */

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

int syscon_read_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(20,a,b,c);
   return fail;
}

int syscon_read_Laurent_system ( void )
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

int syscon_read_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(380,a,b,c);
   return fail;
}

int syscon_write_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(21,a,b,c);
   return fail;
}

int syscon_write_Laurent_system ( void )
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

int syscon_write_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(381,a,b,c);
   return fail;
}

int syscon_number_of_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(22,length,b,c);
   return fail;
}

int syscon_number_of_Laurentials ( int *length )
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

int syscon_number_of_quaddobl_polynomials ( int *length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(382,length,b,c);
   return fail;
}

int syscon_initialize_number ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(23,&length,b,c);
   return fail;
}

int syscon_initialize_number_of_Laurentials ( int length )
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

int syscon_initialize_number_of_quaddobl_polynomials ( int length )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(383,&length,b,c);
   return fail;
}

int syscon_degree_of_polynomial ( int k, int *d )
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

int syscon_store_polynomial ( int nc, int n, int k, char *p )
{
   int a[3],b[nc],i,fail;
   double *c;

   a[0] = nc;
   a[1] = n;
   a[2] = k;
   for(i=0; i<nc; i++) b[i] = (int) p[i];

   fail = _ada_use_c2phc(76,a,b,c);

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

int syscon_load_polynomial ( int k, int *nc, char *p )
{
   int fail,i;
   int b[25600];
   int size = k;
   double *c;

   fail = _ada_use_c2phc(67,&size,b,c);
   /* printf("number of characters : %d\n",size); */
   for(i=0; i<size; i++) p[i] = (char) b[i];
   p[size] = '\0';
   /* printf("the string : %s\n",p); */

   return fail;
}

int syscon_store_Laurential ( int nc, int n, int k, char *p )
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

int syscon_number_of_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(24,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_number_of_Laurent_terms ( int i, int *nt )
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

int syscon_number_of_quaddobl_terms ( int i, int *nt )
{
   int a[2],*b,fail;
   double *c;
   a[1] = i;
   fail = _ada_use_c2phc(384,a,b,c);
   *nt = a[0];
   return fail;
}

int syscon_retrieve_term ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc(25,a,exp,c);  
   return fail;
}

int syscon_retrieve_dobldobl_term
 ( int i, int j, int n, int *exp, double *c )
{
   int a[3],fail;
   a[0] = n;
   a[1] = i;
   a[2] = j;
   fail = _ada_use_c2phc(335,a,exp,c);  
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

int syscon_add_term ( int i, int n, int *exp, double *c )
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

int syscon_add_quaddobl_term ( int i, int n, int *exp, double *c )
{
   int a[2],fail;
   a[0] = n;
   a[1] = i;
   fail = _ada_use_c2phc(386,a,exp,c); 
   return fail;
}

int syscon_total_degree ( int *d )
{
   int *b,fail;
   double *c;
   fail = _ada_use_c2phc(28,d,b,c);
   return fail;
}

int syscon_clear_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(27,a,b,c);
   return fail;
}

int syscon_clear_Laurent_system ( void )
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

int syscon_clear_quaddobl_system ( void )
{
   int *a,*b,fail;
   double *c;
   fail = _ada_use_c2phc(387,a,b,c);
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
