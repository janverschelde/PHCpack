/* file schubert.c contains definitions of the prototypes of schubert.h */

#include <stdio.h>
#include "schubert.h"

int Pieri_root_count ( int m, int p, int q, int *r )
{
   int fail;
   int mpq[3];
   double *c;

   mpq[0] = m; mpq[1] = p; mpq[2] = q;
   
   fail = _ada_use_c2phc4c(223,mpq,r,c,0);
 
   return fail;
}

int resolve_Schubert_conditions
 ( int n, int k, int c, int *brackets, int verbose, int *r )
{
   int fail;
   double rc;
   int dim[4];

   dim[0] = n;
   dim[1] = k;
   dim[2] = c;
   dim[3] = verbose;

   fail = _ada_use_c2phc4c(228,dim,brackets,&rc,0);

   (*r) = (int) rc;

   return fail;
}

int standard_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets,
   int verbose, int verify, int minrep, int tosquare,
   int nbchar, char *filename, int *r, double *flags )
{
   int fail,i;
   int size = 2*(c-2)*n*n+1;
   double rc[size]; /* filename on input, count & flags on return */
   int dim[8];

   dim[0] = n;
   dim[1] = k;
   dim[2] = c;
   dim[3] = verbose;
   dim[4] = verify;
   dim[5] = nbchar;
   dim[6] = minrep;
   dim[7] = tosquare;

   for(i=0; i<nbchar; i++) rc[i] = (double) filename[i];

   fail = _ada_use_c2phc4c(229,dim,brackets,rc,0);

   (*r) = (int) rc[0];

   for(i=1; i<size; i++) flags[i-1] = rc[i];

   return fail;
}

int dobldobl_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets,
   int verbose, int verify, int minrep, int tosquare,
   int nbchar, char *filename, int *r, double *flags )
{
   int fail,i;
   int size = 4*(c-2)*n*n+1;
   double rc[size]; /* filename on input, count & flags on return */
   int dim[8];

   dim[0] = n;
   dim[1] = k;
   dim[2] = c;
   dim[3] = verbose;
   dim[4] = verify;
   dim[5] = nbchar;
   dim[6] = minrep;
   dim[7] = tosquare;

   for(i=0; i<nbchar; i++) rc[i] = (double) filename[i];

   fail = _ada_use_c2phc4c(180,dim,brackets,rc,0);

   (*r) = (int) rc[0];

   for(i=1; i<size; i++) flags[i-1] = rc[i];

   return fail;
}

int quaddobl_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets,
   int verbose, int verify, int minrep, int tosquare,
   int nbchar, char *filename, int *r, double *flags )
{
   int fail,i;
   int size = 8*(c-2)*n*n+1;
   double rc[size]; /* filename on input, count & flags on return */
   int dim[8];

   dim[0] = n;
   dim[1] = k;
   dim[2] = c;
   dim[3] = verbose;
   dim[4] = verify;
   dim[5] = nbchar;
   dim[6] = minrep;
   dim[7] = tosquare;

   for(i=0; i<nbchar; i++) rc[i] = (double) filename[i];

   fail = _ada_use_c2phc4c(181,dim,brackets,rc,0);

   (*r) = (int) rc[0];

   for(i=1; i<size; i++) flags[i-1] = rc[i];

   return fail;
}

int localization_poset ( int m, int p, int q, int *nc, char *ps )
{
   int i,fail;
   int mpq[3];
   int buffer_size = 10240;  /* must be calculated!!! */
   int b[buffer_size];
   double *c;

   mpq[0] = m; mpq[1] = p; mpq[2] = q;
   
   fail = _ada_use_c2phc4c(224,mpq,b,c,0);
   *nc = mpq[0];
   for(i=0; i<*nc; i++) ps[i] = (char) b[i];
   ps[*nc] = '\0';
 
   return fail;
}

int scan_input_planes
 ( int n, int mdim, char *A, double *cff )
/*
 * Scans the string A for n input m-planes.
 * Every m-plane has as many as mdim double coefficients.
 * The total number of coefficients equals n*mdim.
 * All coefficients are stored in the array cff,
 * which must be allocated for n*mdim doubles. */
{
   int i;
   double c;
   char *ptr;

   ptr = A;
   for(i=0; i<n*mdim; i++)         /* scan the planes */
   {
      while(*ptr==' ') ptr++;      /* skip spaces */
      sscanf(&ptr[0],"%lf",&c);    /* read the number */
      cff[i] = c;                  /* store the number */
      if(i < n*mdim-1)
         while(*ptr!=' ') ptr++;   /* skip the number */
   }

   return 0;
}

int scan_complex_interpolation_points
 ( int n, char *points, double *cff )
/*
 * Scans the string points for n interpolation points.
 * Every point is a complex number.
 * All coefficients are stored in the array cff,
 * which must be allocated for 2*n doubles. */
{
   int i;
   double c;
   char *ptr;

   ptr = points;
   for(i=0; i<2*n; i++)           /* scan the doubles */
   {
      while(*ptr==' ') ptr++;     /* skip spaces */
      sscanf(&ptr[0],"%lf",&c);   /* read the number */
      cff[i] = c;                 /* store the number */
      if(i<2*n-1)
         while(*ptr!=' ') ptr++;  /* skip the number */
   }

   return 0;
}

int scan_real_interpolation_points
 ( int n, char *points, double *cff )
/*
 * Scans the string points for n interpolation points.
 * Every point is one double.
 * All coefficients are stored in the array cff,
 * which must be allocated for n doubles. */
{
   int i;
   double c;
   char *ptr;

   ptr = points;
   for(i=0; i<n; i++)             /* scan the doubles */
   {
      while(*ptr==' ') ptr++;     /* skip spaces */
      sscanf(&ptr[0],"%lf",&c);   /* read the number */
      cff[i] = c;                 /* store the number */
      if(i < n-1)
         while(*ptr!=' ') ptr++;  /* skip the number */
   }

   return 0;
}

int show_input_planes
 ( int n, int mdim, double *cff )
/*
 * Prints the coefficients of the input m-planes
 * to screen, for testing purposes.
 * The array cff contains n*mdim doubles. */
{
   int i,j;
   int cnt = 0;

   for(i=0; i<n; i++)
   {
      printf("input plane %d :\n",i);
      for(j=0; j<mdim; j++)
      {
         printf("  %.17e",cff[cnt]);
         cnt = cnt + 1; 
         if(cnt % 2 == 0) printf("\n");
      }
      printf("\n");
   }
   return 0;
}

int show_interpolation_points ( int n, double *cff )
/*
 * Prints the coefficients of the interpolation
 * points to screen, for testing purposes.
 * The array cff contains 2*n doubles. */
{
   int i;

   printf("the interpolation points :\n");
   for(i=0; i<2*n; i++)
   {
      printf("  %.17e",cff[i]);
      if((i+1) % 2 == 0) printf("\n");
   }

   return 0;
}

int Pieri_polynomial_system
 ( int m, int p, int q, int nc, char *A, int is_real )
{
   int i,fail;
   int n = m*p + q*(m+p);
   int mdim = m*(m+p);
   int rdim = n*mdim;
   int tdim = 2*rdim;
   int mpq[3];
   int *b;
   double rcf[rdim];
   double cff[tdim];

   if(is_real == 1)
   {
      fail = scan_input_planes(n,mdim,A,rcf);
      for(i=0; i<rdim; i++)
      {
         cff[2*i] = rcf[i];
         cff[2*i+1] = 0.0;
      }
   }
   else
      fail = scan_input_planes(n,2*mdim,A,cff);

   /* fail = show_input_planes(n,mdim,cff); */

   mpq[0] = m; mpq[1] = p; mpq[2] = q;
   fail = _ada_use_c2phc4c(227,mpq,b,cff,0);

   return fail;
}

int run_Pieri_homotopies
 ( int m, int p, int q, int nc, int *r, char *A, char *pts )
{
   int i,fail;
   int n = m*p + q*(m+p);
   int mdim = 2*m*(m+p);
   int tdim = n*mdim;
   int mpq[3];
   double cff[tdim+2*n];
   double cffpts[2*n];
   int cnt = 0;

   fail = scan_input_planes(n,mdim,A,cff);

   /* fail = show_input_planes(n,mdim,cff); */
   if(q > 0)
   {
      fail = scan_complex_interpolation_points(n,pts,cffpts);
      /* fail = show_interpolation_points(n,cffpts); */
      for(i=0; i<2*n; i++) cff[tdim+i] = cffpts[i];
   }
   mpq[0] = m; mpq[1] = p; mpq[2] = q;

   fail = _ada_use_c2phc4c(225,mpq,r,cff,0);

   return fail;
}

int pack_coefficients ( int n, double *c, int *np, char *p )
/*
 * Writes the n double coefficients in c into the string p,
 * using as many as *np characters. */
{
   int i,j;
   const int size = 30;
   char buf[size];
   int index = 0;

   for(i=0; i<n; i++)
   {
      for(j=0; j<size; j++) buf[j] = ' ';
      sprintf(buf,"%.17e",c[i]);
      for(j=0; j<size; j++)
      {
         if((buf[j] == '\0') || (buf[j] == ' ')) break;
         p[index++] = buf[j];
      }
      if(i < n-1)
         p[index++] = ' ';
      else
         p[index++] = '\0';
   }
   *np = index-1;
   /* printf("number of characters : %d\n",*np);
   printf("the packed coefficients:\n");
   for(i=0; i<(*np); i++) printf("%c",p[i]);
   printf("\n"); */

   return 0;
}

int real_osculating_planes
 ( int m, int p, int q, int *nc, char *s, char *planes )
{
   int i,fail;
   const int n = m*p + q*(m+p);
   const int d = m+p;
   const int e = n*m*d;
   double cffpts[n];
   double coeffs[e];
   int mpq[3];
   int *b;

   mpq[0] = m; mpq[1] = p; mpq[2] = q;

   fail = scan_real_interpolation_points(n,s,cffpts);

   for(i=0; i<n; i++) coeffs[i] = cffpts[i];
   fail = _ada_use_c2phc4c(226,mpq,b,coeffs,0);

   /* printf("the coefficients in schubert :\n"); 
   for(i=0; i<e; i++) printf("%.17lf\n",coeffs[i]); */

   fail = pack_coefficients(e,coeffs,nc,planes);

   /* printf("packed into one string : %s\n",planes); */

   return fail;
}
