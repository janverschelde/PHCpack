/* simple test program in C on the operations in the cells container */

#include <stdio.h>
#include <stdlib.h>
#include "phcpack.h"
#include "solcon.h"
#include "syscon.h"
#include "celcon.h"
#include "intcelcon.h"

int compute_mixed_volume ( void );
/*
 * DESCRIPTION :
 *   prompts the user for a polynomial system and then computes
 *   its mixed volume. */

int write_lifted_supports ( int n );
/* 
 * DESCRIPTION :
 *   writes the lifted supports, the points are n-vectors to screen,
 *   returns 0 if no failure occurred, otherwise 1 is returned. */

int query_cell ( int n, int r );
/*
 * DESCRIPTION :
 *   prompts the user for a cell number and lists the cell,
 *   where r is the number of different supports;
 *   returns 0 if no failure occurred, otherwise 1 is returned. */

int show_mixture ( int dim, int *r );
/* 
 * DESCRIPTION :
 *   queries the cells container for the type of mixture,
 *   the number "dim" on input should be at least as large as the
 *   number of different supports, if the query results in no failure,
 *  then the mixture is shown and *r contains #different supports. */

void read_and_retrieve ( void );
/*
 * DESCRIPTION :
 *   reads a mixed-cell configuration, initializes the cells container
 *   and tests the retrieval operations. */

int read_supports ( int r, int n );
/*
 * DESCRIPTION :
 *   reads r different lifted supports of dimension n
 *   and stores the lifted supports in the cells container,
 *   returns 1 if a failure occurred, 0 otherwise. */

int read_mixed_cell ( int r, int n );
/*
 * DESCRIPTION :
 *   prompts the user for a mixed cell and stores it in the container,
 *   where r equals the number of different supports and n the dimension,
 *   returns 1 if a failure occurred, 0 otherwise. */

void read_cells_and_create_start_system ( void );
/*
 * DESCRIPTION :
 *   reads a mixed-cell configuration, initializes the cells container
 *   and creates a start system which is then written to screen. */

void solve_standard_start_system ( int len );
/*
 * DESCRIPTION :
 *   Solves a random coefficient system in standard double precions.
 *   The input parameter in the number of cells in the container. */

void solve_dobldobl_start_system ( int len );
/*
 * DESCRIPTION :
 *   Solves a random coefficient system in double double precions.
 *   The input parameter in the number of cells in the container. */

void solve_quaddobl_start_system ( int len );
/*
 * DESCRIPTION :
 *   Solves a random coefficient system in quad double precions.
 *   The input parameter in the number of cells in the container. */

void read_cells_and_solve_start_system ( void );
/* 
 * DESCRIPTION :
 *   reads a mixed-cell configuration and then applies the polyhedral
 *   homotopies to solve a random coefficient system. */

void read_and_construct ( void );
/*
 * DESCRIPTION :
 *   prompts the user to provide data for a mixed-cell configuration. */

int prompt_for_precision ( void );
/*
 * DESCRIPTION :
 *   prompts the user to select a precision level which is returned
 *   as an int: 0 for double, 1 for double double, and 2 for quad double */

void test_intcelcon ( void );
/*
 * DESCRIPTION :
 *   Interactive test on the operations in the containers for cells
 *   in a subdivision induced by an integer valued lifting function. */

void check_intcelcon ( int n, int r, int *mix );
/*
 * DESCRIPTION :
 *   Prompts the user for lifted points of size n,
 *   to fill up r supports with type of mixture in mix. */

void write_integer_lifted_supports ( int n, int r, int *mix );
/*
 * DESCRIPTION :
 *   Given r lifted supports of type of mixture in mix,
 *   the n coordinates of each point are written to screen. */

void write_integer_mixed_cells ( int n );
/*
 * DESCRIPTION :
 *   Writes the mixed cells in the regular subdivision
 *   induced by an integer lifting, where n is the length
 *   of the lifted points. */

int main ( void )
{
   char ans;

   printf("\nTesting the Operations in the Cells Container...\n\n");

   adainit();

   printf("Choose one of the following testing options :\n");
   printf("  0. read a polynomial system and compute its mixed volume;\n");
   printf("  1. read a mixed-cell configuration, test retrievals;\n");
   printf("  2. create a mixed-cell configuration interactively;\n");
   printf("  3. read a mixed-cell configuration and create start system;\n");
   printf("  4. polyhedral homotopies solve a random coefficient system;\n");
   printf("  5. test operations on integer cells container.\n");
   printf("Type 1, 2, 3, 4, or 5 to make your choice : "); 
   ans = getchar();

   if(ans == '0')
      compute_mixed_volume();
   else if(ans == '1')
   {
      scanf("%c",&ans);        /* skip end of line symbol */
      read_and_retrieve();
   }
   else if(ans == '2')
      read_and_construct();
   else if(ans == '3')
   {
      scanf("%c",&ans);
      read_cells_and_create_start_system();
   }
   else if(ans == '4')
   {
      scanf("%c",&ans);
      read_cells_and_solve_start_system();
   }
   else if(ans == '5')
   {
      scanf("%c",&ans);
      test_intcelcon();
   }
   else
      printf("\nOption %c unknown.  Please come again.\n\n", ans);

   adafinal();

   return 0;
}

int compute_mixed_volume( void )
{
   int fail,mv,len,dim,r,choice;

   printf("\nMENU for mixed volume calculator :\n");
   printf("  1. use MixedVol;\n");
   printf("  2. use DEMiCs.\n");
   printf("Type 1 or 2 to choose : ");
   scanf("%d",&choice);

   fail = syscon_read_standard_system();
   printf("\nThe system in the container : \n");
   fail = syscon_write_standard_system();
   if(choice == 1)
      fail = mixed_volume(&mv);
   else
      fail = mixed_volume_by_demics(&mv);
   printf("\nThe mixed volume : %d\n",mv);
   fail = celcon_number_of_cells(&len);
   printf("\nnumber of mixed cells : %d\n",len);
   fail = celcon_dimension_of_points(&dim);
   printf("dimension of the lifted points : %d\n",dim);
   fail = show_mixture(dim,&r);
}

int write_lifted_supports ( int n )
{
   int fail,a,b,r,i,j,k;
   int nl[n];
   double *c,pt[n];
 
   fail = celcon_length_of_supports(&r,nl);
   printf("\nlength of the lifted support lists :");
   for(i=0; i<r; i++) printf(" %d",nl[i]); printf("\n");
   printf("the lifted supports:\n");
   for(i=0; i<r; i++)
   {
      for(j=0; j<nl[i]; j++)
      {
        /* printf("retrieving point %d from list %d ...\n",j+1,i+1); */
         fail = celcon_get_lifted_point(n,i+1,j+1,pt);
         for(k=0; k<n-1; k++) printf(" %d",(int) pt[k]);
         printf(" %.15le\n", pt[n-1]);
      }
      printf("\n");
   }
   return fail;
}

int query_cell ( int n, int r )
{
   int k,fail,*b,i,j,mv,nl[r],kk;
   double *c,normal[n],pt[n];

   printf("\nGive a cell number : "); scanf("%d", &k);

   fail = celcon_get_inner_normal(n,k,normal);
   if (fail==1)
      printf("an error happened...\n");
   else
   {
      printf("inner normal for cell %d :\n", k);
      for(i=0;i<n;i++) printf("%.15le\n", normal[i]);
      fail = celcon_number_of_points_in_cell(k,r,nl);
      printf("number of points in supports :");
      for(i=0;i<r;i++) printf(" %d",nl[i]); printf("\n");
      printf("points in the supports :\n");
      for(i=0;i<r;i++)
      {
         for(j=0;j<nl[i];j++)
         {
            fail = celcon_get_point_in_cell(n,k,i+1,j+1,pt);
            for(kk=0; kk<n-1; kk++) printf(" %d",(int) pt[kk]);
            printf(" %.15le\n", pt[n-1]);
         }
      }
      fail = celcon_mixed_volume(k,&mv);
      printf("mixed volume : %d\n",mv);
      {
         int cl[1+r+2*n];
         double inner_normal[n];

         fail = celcon_retrieve_mixed_cell(n,r,k,cl,inner_normal);
         printf("the inner normal for cell %d : \n", k);
         for(i=0; i<n; i++) printf(" %.15le\n", inner_normal[i]);
         printf("total number of points in supports : %d\n",cl[0]);
         printf("number of points in supports : ");
         for(i=1; i<=r; i++) printf(" %d", cl[i]); printf("\n");
         printf("labels of points : ");
         kk=r;
         for(i=1; i<=r; i++)
         {
            for(j=0; j<cl[i]; j++) printf(" %d", cl[++kk]);
            printf(" | ");
         }
         printf("\n");
      }
   }
   return fail;
}

int show_mixture ( int dim, int *r )
{
   int fail,*mix,i;
   double *c;

   mix = (int*)calloc(dim,sizeof(int));
   fail = celcon_type_of_mixture(r,mix);

   if(fail == 1)
      printf("An error occurred, type of mixture not available.\n");
   else
   {
      printf("number of different supports : %d\n",*r);
      printf("type of mixture :");
      for(i=0; i<*r; i++) printf(" %d", mix[i]); printf("\n");
   }

   return fail;
}

void read_and_retrieve ( void )
{
   int n,fail,*d,r;
   double *c;
   int len,dim;
   char ans;

   fail = celcon_read_mixed_cell_configuration();
   printf("Do you wish to see the cells ? (y/n) ");
   ans = getchar();
   if(ans == 'y')
      fail = celcon_write_mixed_cell_configuration();

   fail = celcon_number_of_cells(&len);
   printf("\nnumber of mixed cells : %d\n",len);
   fail = celcon_dimension_of_points(&dim);
   printf("dimension of the lifted points : %d\n",dim);

   fail = show_mixture(dim,&r);

   fail = write_lifted_supports(dim);
   fail = query_cell(dim,r);
}

int read_supports ( int r, int n )
{
   int cnt[r],i,j,k,fail;
   double x[n];

   for(i=0; i<r; i++)
   {
      printf("  number of points in support %d : ",i+1);
      scanf("%d",&cnt[i]);
   }

   for(i=0; i<r; i++)
     for(j=0; j<cnt[i]; j++)
     {
        printf("Reading point %d of support %d...\n",j+1,i+1);
        printf("  give %d coordinates : ",n-1);
        for(k=0; k<n-1; k++) scanf("%lf", &x[k]);
        printf("  give lifting : "); scanf("%lf",&x[n-1]);
        /* printf("  the lifted point is ");
           for(k=0; k<n; k++) printf(" %lf",x[k]); printf("\n"); */
        fail = celcon_append_lifted_point(n,i+1,x);
        if (fail != 0) return fail;
     }

   return fail;
}

int read_mixed_cell ( int r, int n )
{
/* a mixed cell is represented by three vectors:
   1) a[0] = r, a[1] = n, and a[2] = length of vector b;
   2) b[0] = total number of points in the cell
      b[i] = gives the points in the i-th support;
      b[1+r+k] = the label for the k-th lifted point in the cell;
   3) c = an n-vector, is the inner normal */

   int a[3],*b,cnt[r],i,j,k,fail,sum;
   double c[n];

   printf("Give the inner normal (%d-vector) : ",n);
   for(i=0; i<n; i++) scanf("%lf",&c[i]);

   printf("Reading the number of points in each support...\n");
   a[0] = r; a[1] = n;
   sum = 0;
   for(i=0; i<r; i++)
   {
      printf("  #points in support %d : ",i+1);
      scanf("%d",&cnt[i]);
      sum += cnt[i];
   }
   a[2] = 1+r+sum;

   b = (int*)calloc(a[2],sizeof(int));
   b[0] = a[2];
   for(i=0; i<r; i++) b[i+1] = cnt[i];
   printf("Reading the labels for each support...\n");
   k = r+1;
   for(i=0; i<r; i++)
   {
      printf("  give %d labels for support %d : ",cnt[i],i+1);
      for(j=0; j<cnt[i]; j++) scanf("%d",&b[k++]);
   }
   fail = celcon_append_mixed_cell(n,r,b[0],b,c);

   return fail;
}

void read_and_construct ( void )
{
   int r,*mix,i,fail,rr,n;
   double *c;

   printf("\nGive the number of different supports : ");
   scanf("%d",&r);
   mix = (int*)calloc(r,sizeof(int));
   for(i=0; i<r; i++)
   {
      printf("  how many times does support %d occur ? ", i+1);
      scanf("%d",&mix[i]);
   }
   fail = celcon_set_type_of_mixture(r,mix);
   fail = show_mixture(r,&rr);

   printf("\nGive the dimension of the lifted points : ");
   scanf("%d",&n);

   fail = read_supports(r,n);
   fail = write_lifted_supports(n);
   fail = read_mixed_cell(r,n);
   fail = celcon_write_mixed_cell_configuration();
}

void read_cells_and_create_start_system ( void )
{
   int dim,n,fail,*d,r,len,k,nbsols,mv;
   double *c;
   char ans = 'y';

   fail = celcon_read_mixed_cell_configuration();
   printf("\nReading a system to initialize the symbol table...");
   fail = read_standard_target_system();
   fail = celcon_dimension_of_points(&dim);
   printf("dimension of the lifted points : %d\n",dim);
   fail = show_mixture(dim,&r);
   fail = celcon_standard_random_coefficient_system();
   printf("The random coefficient start system :\n");
   fail = celcon_write_standard_random_coefficient_system();
   fail = celcon_standard_polyhedral_homotopy();

   fail = celcon_number_of_cells(&len);
   while (ans == 'y')
   {
      printf("Give a number to a mixed cell (<= %d) : ", len);
      scanf("%d",&k);
      fail = celcon_solve_standard_start_system(k,&nbsols);
      printf(" -> found %d start solutions from cell %d\n",nbsols,k);
      fail = celcon_mixed_volume(k,&mv);
      if (nbsols == mv)
         printf("#start solutions equals mixed volume %d, ok\n",mv);
      else
         printf("#start solutions does not equal mixed volume %d!!!\n",mv);
      printf("Do you wish to test another cell (y/n) ");
      scanf("%c",&ans); /* skip new line symbol */
      ans = getchar();
   }
}

void solve_standard_start_system ( int len ) 
{
   int fail,nb,tnb,k,tmv,mv,i;

   printf("creating a random coefficient system ...\n");

   fail = celcon_standard_random_coefficient_system();
   fail = celcon_standard_polyhedral_homotopy();

   printf("solving the binomial start systems ...\n");
   tnb = 0; tmv = 0;
   for(k=1; k<=len; k++)
   {
      fail = celcon_solve_standard_start_system(k,&nb);
      printf(" -> found %d start solutions from cell %d\n",nb,k);
      fail = celcon_mixed_volume(k,&mv);
      if (nb == mv)
         printf("#start solutions equals mixed volume %d, OK\n",mv);
      else
         printf("#start solutions does not equal mixed volume %d!!!\n",mv);
      tnb += nb; tmv += mv;
   }
   if (tnb == tmv)
      printf("Total #start solutions : %d = %d mixed volume, OK.\n",tnb,tmv);
   else
      printf("Total #start solutions : %d /= %d mixed volume!!!\n",tnb,tmv);

   printf("tracking the solution paths ...\n");
   for(k=1; k<=len; k++)
   {
      fail = celcon_mixed_volume(k,&mv);
      for(i=1; i<=mv; i++)
      {
         printf("Tracking path %d corresponding to cell %d ...\n",i,k);
         fail = celcon_track_standard_solution_path(k,i,0);
      }
   }
   printf("copying the target solutions to the solution container ...\n");
   for(k=1; k<=len; k++)
   {
      fail = celcon_mixed_volume(k,&mv);
      for(i=1; i<=mv; i++)
         fail = celcon_copy_target_standard_solution_to_container(k,i);
   }
   printf("writing random coefficient system and its solutions to file ...\n");
   fail = celcon_write_standard_random_coefficient_system();
   fail = solcon_write_standard_solutions();
}

void solve_dobldobl_start_system ( int len ) 
{
   int fail,nb,tnb,k,tmv,mv,i;

   printf("solving a random coefficient system with double doubles ...\n");

   fail = celcon_dobldobl_random_coefficient_system();
   fail = celcon_dobldobl_polyhedral_homotopy();

   printf("solving the binomial start systems ...\n");
   tnb = 0; tmv = 0;
   for(k=1; k<=len; k++)
   {
      fail = celcon_solve_dobldobl_start_system(k,&nb);
      printf(" -> found %d start solutions from cell %d\n",nb,k);
      fail = celcon_mixed_volume(k,&mv);
      if (nb == mv)
         printf("#start solutions equals mixed volume %d, OK\n",mv);
      else
         printf("#start solutions does not equal mixed volume %d!!!\n",mv);
      tnb += nb; tmv += mv;
   }
   if (tnb == tmv)
      printf("Total #start solutions : %d = %d mixed volume, OK.\n",tnb,tmv);
   else
      printf("Total #start solutions : %d /= %d mixed volume!!!\n",tnb,tmv);

   printf("tracking the solution paths ...\n");
   for(k=1; k<=len; k++)
   {
      fail = celcon_mixed_volume(k,&mv);
      for(i=1; i<=mv; i++)
      {
         printf("Tracking path %d corresponding to cell %d ...\n",i,k);
         fail = celcon_track_dobldobl_solution_path(k,i,0);
      }
   }
   printf("copying the target solutions to the solution container ...\n");
   for(k=1; k<=len; k++)
   {
      fail = celcon_mixed_volume(k,&mv);
      for(i=1; i<=mv; i++)
         fail = celcon_copy_target_dobldobl_solution_to_container(k,i);
   }
   printf("writing random coefficient system and its solutions to file ...\n");
   fail = celcon_write_dobldobl_random_coefficient_system();
   fail = solcon_write_dobldobl_solutions();
}

void solve_quaddobl_start_system ( int len ) 
{
   int fail,nb,tnb,k,tmv,mv,i;

   printf("solving a random coefficient system with quad doubles ...\n");

   fail = celcon_quaddobl_random_coefficient_system();
   fail = celcon_quaddobl_polyhedral_homotopy();

   printf("solving the binomial start systems ...\n");
   tnb = 0; tmv = 0;
   for(k=1; k<=len; k++)
   {
      fail = celcon_solve_quaddobl_start_system(k,&nb);
      printf(" -> found %d start solutions from cell %d\n",nb,k);
      fail = celcon_mixed_volume(k,&mv);
      if (nb == mv)
         printf("#start solutions equals mixed volume %d, OK\n",mv);
      else
         printf("#start solutions does not equal mixed volume %d!!!\n",mv);
      tnb += nb; tmv += mv;
   }
   if (tnb == tmv)
      printf("Total #start solutions : %d = %d mixed volume, OK.\n",tnb,tmv);
   else
      printf("Total #start solutions : %d /= %d mixed volume!!!\n",tnb,tmv);

   printf("tracking the solution paths ...\n");
   for(k=1; k<=len; k++)
   {
      fail = celcon_mixed_volume(k,&mv);
      for(i=1; i<=mv; i++)
      {
         printf("Tracking path %d corresponding to cell %d ...\n",i,k);
         fail = celcon_track_quaddobl_solution_path(k,i,0);
      }
   }
   printf("copying the target solutions to the solution container ...\n");
   for(k=1; k<=len; k++)
   {
      fail = celcon_mixed_volume(k,&mv);
      for(i=1; i<=mv; i++)
         fail = celcon_copy_target_quaddobl_solution_to_container(k,i);
   }
   printf("writing random coefficient system and its solutions to file ...\n");
   fail = celcon_write_quaddobl_random_coefficient_system();
   fail = solcon_write_quaddobl_solutions();
}

void read_cells_and_solve_start_system ( void )
{
   int fail,dim,r,len,prcs;

   fail = celcon_read_mixed_cell_configuration();
   printf("\nReading a system to initialize the symbol table...");
   fail = read_standard_target_system();
   fail = define_output_file();
   fail = celcon_dimension_of_points(&dim);
   printf("dimension of the lifted points : %d\n",dim);
   fail = show_mixture(dim,&r);
   fail = celcon_number_of_cells(&len);
   printf("number of cells in the configuration : %d\n",len);

   prcs = prompt_for_precision();
   if(prcs == 0)
      solve_standard_start_system(len);
   else if(prcs == 1)
      solve_dobldobl_start_system(len);
   else if(prcs == 2)
      solve_quaddobl_start_system(len);
   else
      printf("invalid precision level\n");
}

int prompt_for_precision ( void )
{
   int result;

   printf("\nMENU to select precision : \n");
   printf("  0. standard double precision;\n");
   printf("  1. double double precision;\n");
   printf("  2. quad double precision;\n");
   printf("Type 0, 1, or 2 to select the precision : ");
   scanf("%d", &result);

   return result;
}

void test_intcelcon ( void )
{
   int n,r,i,fail;

   printf("\nTesting the integer cells container ...\n");
   printf("-> give the size of the lifted points : ");
   scanf("%d", &n);
   printf("-> give the number of distinct supports : ");
   scanf("%d", &r);
   {
      int smx;
      int *mix;

      mix = (int*) calloc(r, sizeof(int));

      if(r == n-1)
         for(i=0; i<r; i++) mix[i] = 1;
      else
      {
         for(i=0; i<r; i++)
         {
            printf("-> give frequency of support %d : ", i+1);
            scanf("%d",&mix[i]);
         }
      }
      printf("The mixture : ");
      for(i=0, smx=0; i<r; smx += mix[i++]) printf(" %d", mix[i]);
      printf(" = %d", smx);
      if(smx == n-1)
      {
         printf("  okay\n");
         check_intcelcon(n,r,mix);
         write_integer_lifted_supports(n,r,mix);
         printf("computing the mixed cells ...\n");
         fail = intcelcon_make_subdivision();
         write_integer_mixed_cells(n);
      }
      else
         printf(" != %d\n",n-1);
   }
} 

void check_intcelcon ( int n, int r, int *mix )
{
   int fail,retr,i,j,k,nb;
   int retmix[n-1];

   fail = intcelcon_initialize_supports(r);
   printf("in check_intcelcon, mix = ");
   for(i=0; i<r; i++) printf(" %d", mix[i]);
   printf("\n");

   fail = intcelcon_set_type_of_mixture(r,mix);
   fail = intcelcon_type_of_mixture(&retr,retmix);

   printf("The retrieved number of supports : %d.\n",retr);
   printf("The retrieved mixture : ");
   for(i=0; i<retr; i++) printf(" %d", retmix[i]);
   printf("\n");

   for(i=0; i<r; i++)
   {
      printf("How many lifted points for support %d ? ", i);
      scanf("%d", &nb);
      for(j=0; j<nb; j++)
      {
         int point[n];
         
         printf("Give point %d : ", j);
         for(k=0; k<n; k++) scanf("%d", &point[k]);
         fail = intcelcon_append_lifted_point(n,i+1,point);
      }
   }
}

void write_integer_lifted_supports ( int n, int r, int *mix )
{
   int fail,rr,i,j,k,m;
   int nb[r];

   fail = intcelcon_length_of_supports(&rr,nb);
   printf("Number of different supports retrieved : %d.\n", rr);
   printf("Length of each support :");
   for(i=0; i<r; i++) printf(" %d", nb[i]);
   printf("\n");

   for(i=1; i<=r; i++)
   {
      printf("Points in support %d :\n", i);
      for(j=1; j<=nb[i-1]; j++)
      {
         int point[n];

         printf("-> point %d of %d :", j, nb[i-1]);
         fail = intcelcon_get_lifted_point(n,i,j,point);
         for(k=0; k<n; k++) printf(" %d", point[k]);
         printf("\n");
      }
   }
}

void write_integer_mixed_cells ( int n )
{
   int nbcells,fail,i,j,mv,mixvol;
   int normal[n];

   fail = intcelcon_number_of_cells(&nbcells);
   printf("The number of mixed cells : %d\n",nbcells);
   mixvol = 0;
   for(i=1; i<=nbcells; i++)
   {
      printf("Inner normal for cell %d :",i);
      fail = intcelcon_get_inner_normal(n,i,normal);
      for(j=0; j<n; j++) printf(" %d",normal[j]);
      fail = intcelcon_mixed_volume(i,&mv);
      printf(", volume : %d \n",mv);
      mixvol = mixvol + mv;
   }
   printf("The mixed volume : %d.\n",mixvol);
}
