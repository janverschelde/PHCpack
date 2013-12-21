/* simple test program in C on the Ada procedure use_celcon */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_celcon ( int task, int *a, int *b, double *c );
extern void adafinal();

int write_lifted_supports ( int n );
/* writes the lifted supports, the points are n-vectors to screen,
   returns 0 if no failure occurred, otherwise 1 is returned */

int query_cell ( int n, int r );
/* prompts the user for a cell number and lists the cell,
   where r is the number of different supports;
   returns 0 if no failure occurred, otherwise 1 is returned */

int show_mixture ( int dim, int *r );
/* queries the cells container for the type of mixture,
   the number "dim" on input should be at least as large as the
   number of different supports, if the query results in no failure,
   then the mixture is shown and *r contains #different supports */

void read_and_retrieve ( void );
/* reads a mixed-cell configuration, initializes the cells container
   and tests the retrieval operations */

int read_supports ( int r, int n );
/* reads r different lifted supports of dimension n
   and stores the lifted supports in the cells container,
   returns 1 if a failure occurred, 0 otherwise */

int read_mixed_cell ( int r, int n );
/* prompts the user for a mixed cell and stores it in the container,
   where r equals the number of different supports and n the dimension,
   returns 1 if a failure occurred, 0 otherwise */

void read_and_construct ( void );
/* prompts the user to provide data for a mixed-cell configuration */

int main ( void )
{
   char ans;

   printf("\nTesting the Operations in the Cells Container...\n\n");

   adainit();

   printf("Choose one of the following testing options :\n");
   printf("  1. read a mixed-cell configuration from file and retrievals;\n");
   printf("  2. create a mixed-cell configuration interactively.\n");
   printf("Type 1 or 2 to make your choice : "); 
   ans = getchar();

   if(ans == '1')
   {
      scanf("%c",&ans);        /* skip end of line symbol */
      read_and_retrieve();
   }
   else if(ans == '2')
      read_and_construct();
   else
      printf("\nOption %c unknown.  Please come again.\n\n", ans);

   adafinal();

   return 0;
}

int write_lifted_supports ( int n )
{
   int fail,a,b,r,i,j,k;
   int nl[n];
   double *c,pt[n];
 
   fail = _ada_use_celcon(5,&r,nl,c);
   printf("\nlength of the lifted support lists :");
   for(i=0; i<r; i++) printf(" %d",nl[i]); printf("\n");
   printf("the lifted supports:\n");
   for(i=0; i<r; i++)
   {
      a = i+1;
      for(j=0; j<nl[i]; j++)
      {
         b = j+1;
        /* printf("retrieving point %d from list %d ...\n",b,a); */
         fail = _ada_use_celcon(6,&a,&b,pt);
         for(k=0; k<n-1; k++) printf(" %d",(int) pt[k]);
         printf(" %.15le\n", pt[n-1]);
      }
      printf("\n");
   }
   return fail;
}

int query_cell ( int n, int r )
{
   int k,fail,*b,i,j,mv,nl[r],ind[2],kk;
   double *c,normal[n],pt[n];

   printf("\nGive a cell number : "); scanf("%d", &k);

   fail = _ada_use_celcon(7,&k,b,normal);
   if (fail==1)
      printf("an error happened...\n");
   else
   {
      printf("inner normal for cell %d :\n", k);
      for(i=0;i<n;i++) printf("%.15le\n", normal[i]);
      fail = _ada_use_celcon(8,&k,nl,c);
      printf("number of points in supports :");
      for(i=0;i<r;i++) printf(" %d",nl[i]); printf("\n");
      printf("points in the supports :\n");
      for(i=0;i<r;i++)
      {
         for(j=0;j<nl[i];j++)
         {
            ind[0] = i+1;
            ind[1] = j+1;
            fail = _ada_use_celcon(9,&k,ind,pt);
            for(kk=0; kk<n-1; kk++) printf(" %d",(int) pt[kk]);
            printf(" %.15le\n", pt[n-1]);
         }
      }
      fail = _ada_use_celcon(10,&k,&mv,c);
      printf("mixed volume : %d\n",mv);
      {
         int cl[1+r+2*n];
         double inner_normal[n];

         fail = _ada_use_celcon(15,&k,cl,inner_normal);
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

   mix = calloc(dim,sizeof(int));
   fail = _ada_use_celcon(4,r,mix,c);   /* type of mixture */

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

   fail = _ada_use_celcon(0,&n,d,c);     /* read mixed-cell configuration */
   printf("Do you wish to see the cells ? (y/n) ");
   ans = getchar();
   if(ans == 'y')
      fail = _ada_use_celcon(1,&n,d,c);  /* write mixed-cell configuration */

   fail = _ada_use_celcon(2,&len,d,c);   /* number of mixed cells */
   printf("\nnumber of mixed cells : %d\n",len);
   fail = _ada_use_celcon(3,&dim,d,c);   /* dimension of lifted points */
   printf("dimension of the lifted points : %d\n",dim);

   fail = show_mixture(dim,&r);

   fail = write_lifted_supports(dim);
   fail = query_cell(dim,r);
}

int read_supports ( int r, int n )
{
   int cnt[r],i,j,k,fail,i1;
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
        i1 = i+1;
        fail = _ada_use_celcon(12,&i1,&n,x);
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
   2) b[0] = total number of points in the cell,
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

   b = calloc(a[2],sizeof(int));
   b[0] = sum;
   for(i=0; i<r; i++) b[i+1] = cnt[i];
   printf("Reading the labels for each support...\n");
   k = r+1;
   for(i=0; i<r; i++)
   {
      printf("  give %d labels for support %d : ",cnt[i],i+1);
      for(j=0; j<cnt[i]; j++) scanf("%d",&b[k++]);
   }
   fail = _ada_use_celcon(13,a,b,c);

   return fail;
}

void read_and_construct ( void )
{
   int r,*mix,i,fail,rr,n;
   double *c;

   printf("\nGive the number of different supports : ");
   scanf("%d",&r);
   mix = calloc(r,sizeof(int));
   for(i=0; i<r; i++)
   {
      printf("  how many times does support %d occur ? ", i+1);
      scanf("%d",&mix[i]);
   }
   fail = _ada_use_celcon(11,&r,mix,c);
   fail = show_mixture(r,&rr);

   printf("\nGive the dimension of the lifted points : ");
   scanf("%d",&n);

   fail = read_supports(r,n);
   fail = write_lifted_supports(n);
   fail = read_mixed_cell(r,n);
   fail = _ada_use_celcon(1,&n,&r,c);  /* write mixed-cell configuration */
}
