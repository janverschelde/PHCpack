// the main program to run the polynomial evaluation and differentiation

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <gqd_type.h>
#include "complexD.h"
#include "gqd_qd_util.h"
#include "DefineType.h"
#include "complexH.h"
#include "ped_kernelsT.h"

using namespace std;

#define d  0 
#define dd 1
#define qd 2

#ifdef precision
#define the_p precision
#else
#define the_p 0
#endif

#if(the_p == 0)
typedef double realD;
typedef double realH;
#elif(the_p == 1)
typedef gdd_real realD;
typedef dd_real realH;
#else
typedef gqd_real realD;
typedef qd_real realH;
#endif

void print ( complexD<realD>* a, complexH<realH>* b, int dim, string var );
void print ( complexD<realD>* a, complexH<realH>** b, int dimX, int dimY,
             int stride, string var);

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<realD> *x_h, complexD<realD> *c_h,
   complexD<realD> *factors_h, complexD<realD> *polvalues_h );
/*
 * THE GPU is used to evaluate a system and its Jacobian matrix.
 *
 * ON ENTRY :
 *   BS           number of threads in a block;
 *   dim          dimension of the problem;
 *   NM           number of monomials;
 *   NV           number of variables;
 *   deg          highest degree of the variables;
 *   r            frequency of the runs,
 *   m            NM divided by dim
 *   ncoefs       number of coefficients
 *   pos          indicate the participating variables in the monomials;
 *   exp          the exponents of the participating variables; 
 *   c            coefficients of the monomials;
 *   x            point where to evaluate.
 *
 * ON RETURN :
 *   factors_h    the factors common to evaluation and differentiation;
 *   polvalues_h  are the values of polynomials and their derivatives. */

template <class realH>
void CPU_evaldiff
 ( int dim, int NM, int NV, int deg, int r, int m, int *pos, int *exp,
   complexH<realH> *c, complexH<realH> *x,
   complexH<realH> *factors_h, complexH<realH> **polvalues_h );
/*
 * The CPU is used to evaluate a system and its Jacobian matrix.
 *
 * ON ENTRY :
 *   dim          dimension of the problem;
 *   NM           number of monomials;
 *   NV           number of variables;
 *   deg          highest degree of the variables;
 *   r            frequency of the runs,
 *   m            NM divided by dim
 *   pos          indicate the participating variables in the monomials;
 *   exp          the exponents of the participating variables; 
 *   c            coefficients of the monomials;
 *   x            point where to evaluate.
 *
 * ON RETURN :
 *   factors_h    the factors common to evaluation and differentiation;
 *   polvalues_h  are the values of polynomials and their derivatives. */

int main ( int argc, char *argv[] )
{
   int BS,dim,NM,NV,deg,r,mode;
   if(parse_arguments(argc,argv,&BS,&dim,&NM,&NV,&deg,&r,&mode) == 1) return 1;

   int timevalue = 1287178355; // fixed seed instead of timevalue=time(NULL)
   const int n = 32;
   complexD<realD> *x = (complexD<realD>*)calloc(n,sizeof(complexD<realD>));
   complexD<realD> *y = (complexD<realD>*)calloc(n,sizeof(complexD<realD>));
   complexD<realD> *yD2 = (complexD<realD>*)calloc(n,sizeof(complexD<realD>));
  
   srand(timevalue);
   complexD<realD> *xp_d = new complexD<realD>[dim];
   complexH<realH> *xp_h = new complexH<realH>[dim];
   random_point(dim,xp_d,xp_h);
   // print(xp_d,xp_h,dim,"x");
   // for (int i=0;i<dim;i++) cout << "xp =" << xp_h[i];
   storage<realH,4> pt; pt.copyToStor(xp_h);
   // pt.print(); 
   // cout << "M_PI=" << M_PI << endl;

   int pos_arr_h_int[NM*NV];
   int exp_arr_h_int[NM*NV];
   char pos_arr_h_char[NM*NV];
   char exp_arr_h_char[NM*NV];
   int ncoefs = NM*(NV+1);
   complexD<realD> *c_d = new complexD<realD>[ncoefs];
   complexH<realH> *c_h = new complexH<realH>[ncoefs];
   generate_system(dim,NM,NV,deg,pos_arr_h_int,pos_arr_h_char,
                   exp_arr_h_int,exp_arr_h_char,c_d,c_h);
   if(mode == 2) write_system<realH>(dim,NM,NV,c_h,pos_arr_h_int,exp_arr_h_int);
   // allocate space for output
   int m = NM/dim;
   // cout << "m=" << m << endl;
   complexD<realD> *factors_d = new complexD<realD>[NM];
   complexH<realH> *factors_h = new complexH<realH>[NM];
   int ant = ((dim*dim+dim)/BS + 1)*BS;
   // cout << "ant=" << ant << endl;
   complexD<realD> *polvalues_d = new complexD<realD>[ant];
   complexH<realH> **polvalues_h = new complexH<realH>*[dim];
   for(int i=0; i<dim; i++) polvalues_h[i] = new complexH<realH>[dim];

   if(mode == 0 || mode == 2)
      GPU_evaldiff(BS,dim,NM,NV,deg,r,m,ncoefs,pos_arr_h_char,
                   exp_arr_h_char,xp_d,c_d,factors_d,polvalues_d);

   if(mode == 1 || mode == 2)
   {
      CPU_evaldiff<realH>
         (dim,NM,NV,deg,r,m,pos_arr_h_int,
          exp_arr_h_int,c_h,xp_h,factors_h,polvalues_h);
      if(mode==2) // compare results of GPU with CPU
      {
         // print(factors_h, factors_d,NM,"factors");
         // print(polvalues_d,polvalues_h,dim,dim,dim,"poly_derivatives"); 
         error_on_factors(NM,factors_d,factors_h);
         error_on_derivatives(dim,polvalues_h,polvalues_d);
      }
   }
   for(int i = 0; i<n; i++) x[i].initH(double(i+2),0.12345);      

   int fail; // = sqrt_by_Newton(n,x,y);
   int fail1; //= squareCall(n,y,yD2);

   complexH<realH>* xH;
   complexH<realH>* yH;
   complexH<realH>* yD2H;
   xH = new complexH<realH>[n];
   yH = new complexH<realH>[n];
   yD2H= new complexH<realH>[n];

   comp1_gqdArr2qdArr(x, xH, n);
   comp1_gqdArr2qdArr(y, yH, n);
   comp1_gqdArr2qdArr(yD2, yD2H, n);

   return 0;
}

void print ( complexD<realD>* a, complexH<realH>* b, int dim, string var )
{
   complexH<realH> temp;

   for (int i=0;i<dim;i++)
   {
      comp1_gqd2qd(&a[i],&temp);
      cout << "GPU: " << var << "["<< i << "] = " << temp;
      cout << "CPU: " << var << "["<< i << "] = " << b[i];
   }
}

void print
 ( complexD<realD>* a, complexH<realH>** b, int dimX, int dimY,
   int stride, string var)
{

   complexH<realH> temp;

   for(int i=0;i<dimX;i++)
      for(int j=0;j<dimY;j++)
      {
         comp1_gqd2qd(&a[stride+dimY*j+i],&temp);
         cout << "GPU: " 
              << var << "["<< i << "]" << "["<< j << "] = " << temp;
         cout << "CPU: "
              << var << "["<< i << "]" << "["<< j << "] = " << b[i][j];
      }
}

template <class realH>
void degrees ( complexH<realH> **deg, complexH<realH> *bm, int dim, int maxd )
{
   for(int i=0; i<dim; i++)
   {
      deg[i][0].init(1.0,0.0);
      for(int j=1; j<=maxd; j++) deg[i][j] = deg[i][j-1]*bm[i];
   }
}

template <class realH>
void comp_factors
 ( int na, int no_monomials, int nvarm,
   int **monvarindexes, int **nonz_exponentsx,
   complexH<realH> *xval, complexH<realH> **deg, complexH<realH> *roots )
{
   for(int i=0; i<no_monomials; i++)
   {
      roots[i].init(1.0,0.0);
      for(int j=0; j<nvarm; j++)
         roots[i]=roots[i]*deg[monvarindexes[i][j]][nonz_exponentsx[i][j]];
   }
}

template <class realH>
void speeldif_h
 ( int na, int NM, int nvarm, int nders, int **monvarindexes,
   complexH<realH> *xval, complexH<realH> *roots, complexH<realH> *coefs,
   complexH<realH> *monvalues, complexH<realH> **monderivatives )
{
   complexH<realH> *products1;
   complexH<realH> *products2;
   products1 = new complexH<realH>[nvarm];
   products2 = new complexH<realH>[nvarm];
   complexH<realH> *top_mon_derivatives;
   top_mon_derivatives = new complexH<realH>[nvarm];
   for(int i=0; i<NM; i++)
   {
      products1[0] = xval[monvarindexes[i][0]];
      products2[0] = xval[monvarindexes[i][nvarm-1]];
      for(int j=0; j<(nvarm-2); j++)
      {
       	 products1[j+1] = products1[j]*xval[monvarindexes[i][j+1]];
         products2[j+1] = products2[j]*xval[monvarindexes[i][nvarm-2-j]];
      }
      top_mon_derivatives[0]=products2[nvarm-2];
      top_mon_derivatives[nvarm-1]=products1[nvarm-2];
      for(int j=0; j<(nvarm-2); j++)
         top_mon_derivatives[j+1] = products1[j]*products2[nvarm-3-j];
      for(int j=0; j<nvarm; j++)
         monderivatives[i][j] = roots[i]*top_mon_derivatives[j];
      monvalues[i] = monderivatives[i][0]*xval[monvarindexes[i][0]];
      for(int j=0; j<nvarm; j++)
         monderivatives[i][j] = monderivatives[i][j]*coefs[NM*j+i];
      monvalues[i] = monvalues[i] * coefs[nders+i];
   }
}

template <class realH>
void sum_mons
 ( int dim, int m, int nvarm, int **monvarindexes,
   complexH<realH> *monvalues, complexH<realH> **monderivatives,
   complexH<realH> *polyvalues, complexH<realH> **polyderivatives )
{
   for(int i=0; i<dim; i++)
   {
      polyvalues[i].init(0.0,0.0);
      for(int j=0; j<dim; j++)
          polyderivatives[i][j].init(0.0,0.0);
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<m; j++)
      {
       	 int mon_number = m*i+j;
         polyvalues[i] = polyvalues[i] + monvalues[mon_number];
         for(int k=0; k<nvarm ; k++)
             polyderivatives[i][monvarindexes[mon_number][k]]
                = polyderivatives[i][monvarindexes[mon_number][k]]
                + monderivatives[mon_number][k];
      }
}

template <class realH>
void CPU_evaldiff
 ( int dim, int NM, int NV, int deg, int r, int m, int *pos, int *exp,
   complexH<realH> *c, complexH<realH> *x,
   complexH<realH> *factors_h, complexH<realH> **polvalues_h )
{
   complexH<realH> **vdegrees = new complexH<realH>*[dim];
   for(int i=0; i<dim; i++) vdegrees[i] = new complexH<realH>[deg+1];
   int **positions_h = new int*[NM];
   int **exponents_h = new int*[NM];
   for(int i=0; i<NM; i++)
   {
      positions_h[i] = new int[NV];
      exponents_h[i] = new int[NV];
   }
   for(int i=0; i<NM; i++)
      for(int j=0; j<NV; j++) positions_h[i][j] = pos[NV*i+j];
   for(int i=0; i<NM; i++)
      for(int j=0; j<NV; j++) exponents_h[i][j] = exp[NV*i+j];
   complexH<realH> *monvalues = new complexH<realH>[NM];
   complexH<realH> **derivatives = new complexH<realH>*[NM];
   for(int i=0; i<NM; i++) derivatives[i] = new complexH<realH>[NV];
   int nders = NM*NV;
   complexH<realH> *polyvalues = new complexH<realH>[dim];
   for(int j=0; j<r; j++)
   {
      degrees<realH>(vdegrees,x,dim,deg);
      comp_factors<realH>
         (dim,NM,NV,positions_h,exponents_h,x,vdegrees,factors_h);
      speeldif_h<realH>
         (dim,NM,NV,nders,positions_h,x,factors_h,c,monvalues,derivatives);
      sum_mons<realH>
         (dim,m,NV,positions_h,monvalues,derivatives,polyvalues,polvalues_h);
   }
}
