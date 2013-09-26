/* simple test program in C on the Ada procedure use_solcon */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_solcon ( int task, int *a, int *b, double *c );
extern void adafinal();

typedef struct list LIST;
struct list
{
   double rt;   /* real part of complex continuation parameter t */
   double it;   /* imaginary part of t */
   int m;       /* multiplicity of the solution vector */
   double *v;   /* solution is sequence of real and imaginary parts */
   double err;  /* error on the correction term from Newton's method */
   double rco;  /* inverse of estimate of condition number */
   double res;  /* residual of the solution vector */
   LIST *next;  /* next solution in the list */
};

LIST *push ( LIST *l, int n, double rt, double it, int m, double *v,
             double err, double rco, double res );
/* pushes the solution to the front of the list l */

void write_solution ( int n, double rt, double it, int m, double *v,
                      double err, double rco, double res );
/* writes the n-dimensional solution to the screen */

void write_solutions ( LIST *l, int n );
/* writes the list of n-dimensional solution vectors */

void test_solution_retrieval ( int n );
/* prompts the user for some solution number,
 * retrieves and prints the solution vector,
 * n is the complex dimension of the solution */

void test_replace_solution ( int n );
/* prompts the user for some solution number,
 * and replaces the solution */

LIST *test_retrievals ( int len, int n );
/* writes the solution list of length len and complex dimension n
 * on standard output, returns the list of solutions */

void test_additions ( int n, LIST *l );
/* fills up the container with the solutions in the list l */

void test_standard_solution_strings ( int len, int nvar );
/* tests operations on solution strings of dimension nvar,
 * where len is the size of the solution container */

void test_dobldobl_solution_strings ( int len, int nvar );
/* tests operations on double double solution strings of dimension nvar,
 * where len is the size of the solution container */

void test_quaddobl_solution_strings ( int len, int nvar );
/* tests operations on double double solution strings of dimension nvar,
 * where len is the size of the solution container */

void test_standard_solution_container ( void );
/* tests main operations in the solutions container */

void test_dobldobl_solution_container ( void );
/* tests main operations in the solutions container
 * for complex double double coordinates  */

void test_quaddobl_solution_container ( void );
/* tests main operations in the solutions container
 * for complex double double coordinates  */

void test_incremental_read_and_write ( void );
/* tests incremental read and write of solutions */

int main(void)
{
   int option;
   char ch;

   printf("\nMENU to test operations in solutions container :\n");
   printf("  1. test main operations in the container;\n");
   printf("  2. incremental read/write of solutions from/to file;\n");
   printf("  3. test container for double double solutions;\n");
   printf("  4. test container for quad double solutions.\n");
   printf("Type 1, 2, 3, or 4 to make your choice : "); 
   scanf("%d",&option);
   scanf("%c",&ch); /* skip new line character */
   printf("\n");

   adainit();
   if(option == 1)
      test_standard_solution_container();
   else if(option == 2)
      test_incremental_read_and_write();
   else if(option == 3)
      test_dobldobl_solution_container();
   else if(option == 4)
      test_quaddobl_solution_container();
   else
      printf("%d is wrong choice.  Please try again...\n");
   adafinal();

   return 0;
}

void retrieve_and_write ( int n, int k )
{
   int i,m,fail;
   double sol[2*n+5];

   fail = _ada_use_solcon(4,&k,&m,sol);

   if (fail == 1)
      printf("some failure occurred...\n");
   else
   { 
      printf("The solution vector :\n");
      for(i=0; i<n; i++)
         printf(" %.15e  %.15e\n", sol[2+2*i], sol[2+2*i+1]);
   }
}

void test_solution_retrieval ( int n )
{
   int k;

   printf("Give k : "); scanf("%d", &k);

   retrieve_and_write(n,k);
}

void test_replace_solution ( int n )
{
   int i,k,m[2],fail;
   double sol[2*n+5];

   printf("Give k : "); scanf("%d", &k);
   retrieve_and_write(n,k);

   printf("Setting solution vector to 1,2,3,...\n");
   sol[0] = 1.111;
   sol[1] = 2.222;
   m[0] = n;
   m[1] = 12345;
   sol[2*n+2] = 3.333;
   sol[2*n+3] = 4.444;
   sol[2*n+4] = 5.555; 
   for(i=0; i<2*n; i++) sol[2+i] = i;

   fail = _ada_use_solcon(5,&k,m,sol);

   if (fail == 1)
      printf("some failure occurred...\n");
   else
   { 
      printf("The changed solution retrieved again :\n");
      retrieve_and_write(n,k);
   }
}

void write_solution ( int n, double rt, double it, int m, double *v,
                      double err, double rco, double res )
{
   int i;

   printf("t : %.15e  %.15e\n", rt, it);
   printf("m : %d\n", m);
   printf("The solution for t :\n");
   for(i=0; i<n; i++)
      printf(" %.15e  %.15e\n", v[2*i], v[2*i+1]);
   printf("== err : %.3e = rco : %.3e = res : %.3e ==\n",err,rco,res);
}

LIST *push ( LIST *l, int n, double rt, double it, int m, double *v,
             double err, double rco, double res )
{
   int i;
   LIST *nl = (LIST*)calloc(1,sizeof(LIST));

   nl->rt = rt;
   nl->it = it;
   nl->m = m;
   nl->v = calloc(2*n,sizeof(double));
   for(i=0; i<2*n; i++)
      nl->v[i] = v[i];
   nl->err = err;
   nl->rco = rco;
   nl->res = res;
   nl->next = l;

   return nl;
}

void write_solutions ( LIST *l, int n )
{
   LIST *p;
   int k;

   for(p=l,k=1; p!=NULL; p=p->next,k++)
   {
      printf("Solution %d :\n",k);
      write_solution(n,p->rt,p->it,p->m,p->v,p->err,p->rco,p->res);
   }
}

LIST *test_retrievals ( int len, int n )
{
   int i,k,m,fail;
   double sol[2*n+5],err,rco,res;
   LIST *l = NULL;

   for(k=1; k<=len; k++)
   {
      fail = _ada_use_solcon(4,&k,&m,sol);
      printf("Solution %d :\n", k);
      err = sol[2*n+2];
      rco = sol[2*n+3];
      res = sol[2*n+4];
      write_solution(n,sol[0],sol[1],m,sol+2,err,rco,res);
      l = push(l,n,sol[0],sol[1],m,sol+2,err,rco,res);
   }

   return l;
}

void test_additions ( int n, LIST *l )
{
   int i,k,m[2],fail;
   LIST *p;
   double sol[2*n+5];

   m[0] = n;
   for(p=l,k=1; p!=NULL; p=p->next,k++)
   {
      sol[0] = p->rt;
      sol[1] = p->it;
      m[1] = p->m;
      for(i=0;i<2*n;i++) sol[2+i] = p->v[i];
      sol[2*n+2] = p->err;
      sol[2*n+3] = p->rco;
      sol[2*n+4] = p->res;
      fail = _ada_use_solcon(6,&k,m,sol);
   }
}

void test_standard_solution_strings ( int len, int nvar )
{
   int k,n,fail;
   double *c;

   printf("\nTesting strings of solutions ... \n\n");

   printf("Give index of solution ( in [1,%d] ) : ",len);
   scanf("%d",&k);
   fail = _ada_use_solcon(30,&k,&n,c);
   printf("-> number of characters in solution %d : %d\n",k,n);
   {
      int i,a[2];
      int v[n];
      char s[n+1];

      a[0] = k; a[1] = n;
      fail = _ada_use_solcon(31,a,v,c);
      for(i=0; i<n; i++) s[i] = (char) v[i];
      s[n] = '\0';
      printf("the solution %d :\n%s\n",k,s);
      a[0] = nvar; a[1] = n;
      fail = _ada_use_solcon(38,a,v,c);
      printf("the container after appending solution %d :\n",k);
      fail = _ada_use_solcon(1,a,v,c);
   }
}

void test_dobldobl_solution_strings ( int len, int nvar )
{
   int k,n,fail;
   double *c;

   printf("\nTesting strings of solutions ... \n\n");

   printf("Give index of solution ( in [1,%d] ) : ",len);
   scanf("%d",&k);
   fail = _ada_use_solcon(70,&k,&n,c);
   printf("-> number of characters in solution %d : %d\n",k,n);
   {
      int i,a[2];
      int v[n];
      char s[n+1];

      a[0] = k; a[1] = n;
      fail = _ada_use_solcon(71,a,v,c);
      for(i=0; i<n; i++) s[i] = (char) v[i];
      s[n] = '\0';
      printf("the solution %d :\n%s\n",k,s);
      a[0] = nvar; a[1] = n;
      fail = _ada_use_solcon(78,a,v,c);
      printf("the container after appending solution %d :\n",k);
      fail = _ada_use_solcon(41,a,v,c);
   }
}

void test_quaddobl_solution_strings ( int len, int nvar )
{
   int k,n,fail;
   double *c;

   printf("\nTesting strings of solutions ... \n\n");

   printf("Give index of solution ( in [1,%d] ) : ",len);
   scanf("%d",&k);
   fail = _ada_use_solcon(110,&k,&n,c);
   printf("-> number of characters in solution %d : %d\n",k,n);
   {
      int i,a[2];
      int v[n];
      char s[n+1];

      a[0] = k; a[1] = n;
      fail = _ada_use_solcon(111,a,v,c);
      for(i=0; i<n; i++) s[i] = (char) v[i];
      s[n] = '\0';
      printf("the solution %d :\n%s\n",k,s);
      a[0] = nvar; a[1] = n;
      fail = _ada_use_solcon(118,a,v,c);
      printf("the container after appending solution %d :\n",k);
      fail = _ada_use_solcon(81,a,v,c);
   }
}

void test_standard_solution_container ( void )
{
   int k,len,n,fail;
   double *c;
   LIST *l;

   fail = _ada_use_solcon(0,&k,&k,c);          /* read solution list */
   printf("The solution list :\n");
   fail = _ada_use_solcon(1,&k,&k,c);          /* print list to screen */
   fail = _ada_use_solcon(2,&k,&len,c);        /* get length of list */
   printf("Number of solutions : %d\n",len);
   /* get the dimension of the solutions in the container */
   fail = _ada_use_solcon(3,&k,&n,c);          /* get the dimension */
   printf("Dimension of vectors : %d\n",n);
   test_standard_solution_strings(len,n);
   /* test retrievals of solutions */
   test_solution_retrieval(n);
   test_replace_solution(n);
   l = test_retrievals(len,n);
   printf("The solution list is reverse order : \n");
   write_solutions(l,n);
   fail = _ada_use_solcon(7,&k,&k,c);          /* clear the container */
   test_additions(n,l);                        /* reconstruct */
   printf("The reconstructed list in reverse order :\n");
   fail = _ada_use_solcon(1,&k,&k,c); 
}

void test_dobldobl_solution_container ( void )
{
   int k,len,n,fail;
   double *c;

   printf("Testing container for double double solutions ...\n\n");

   fail = _ada_use_solcon(40,&k,&k,c);          /* read solution list */
   printf("The solution list :\n");
   fail = _ada_use_solcon(41,&k,&k,c);          /* print list to screen */
   fail = _ada_use_solcon(42,&k,&len,c);        /* get length of list */
   printf("Number of solutions : %d\n",len);
   /* get the dimension of the solutions in the container */
   fail = _ada_use_solcon(43,&k,&n,c);          /* get the dimension */
   printf("Dimension of vectors : %d\n",n);
   test_dobldobl_solution_strings(len,n);
}

void test_quaddobl_solution_container ( void )
{
   int k,len,n,fail;
   double *c;

   printf("Testing container for quad double solutions ...\n\n");

   fail = _ada_use_solcon(80,&k,&k,c);          /* read solution list */
   printf("The solution list :\n");
   fail = _ada_use_solcon(81,&k,&k,c);          /* print list to screen */
   fail = _ada_use_solcon(82,&k,&len,c);        /* get length of list */
   printf("Number of solutions : %d\n",len);
   /* get the dimension of the solutions in the container */
   fail = _ada_use_solcon(83,&k,&n,c);          /* get the dimension */
   printf("Dimension of vectors : %d\n",n);
   test_quaddobl_solution_strings(len,n);
}

void test_incremental_read_and_write ( void )
{
   int fail,*a,*b,len,dim,k,m[2],i;
   double *c,err,rco,res;
   char ans;

   fail = _ada_use_solcon(10,a,b,c);
   printf("\n");
   fail = _ada_use_solcon(11,a,b,c);
   printf("\nShould the file be scanned for a banner ? (y/n) ");
   scanf("%c",&ans);
   if (ans == 'y') fail = _ada_use_solcon(12,a,b,c);
   fail = _ada_use_solcon(13,&len,&dim,c);
   printf("#solutions : %d  dimension : %d\n",len,dim);
   fail = _ada_use_solcon(14,&len,&dim,c);
   c = (double*)calloc(2*dim+5,sizeof(double));
   k = 0;
   for(i=0; i<len; i++)
   {
      fail = _ada_use_solcon(15,&dim,&m[0],c);
      printf("Solution %d :\n", i+1);
      err = c[2*dim+2]; rco = c[2*dim+3]; res = c[2*dim+4];
      write_solution(dim,c[0],c[1],m[0],c+2,err,rco,res);
      m[1] = m[0]; m[0] = dim;
      fail = _ada_use_solcon(16,&k,m,c);
      printf("number of solutions written to file : %d\n",k);
   }
   fail = _ada_use_solcon(17,a,b,c); /* close input file */
   fail = _ada_use_solcon(18,a,b,c); /* close output file */
}
