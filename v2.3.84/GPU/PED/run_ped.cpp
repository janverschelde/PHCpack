#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <gqd_type.h>
#include "complex.h"
#include "gqd_qd_util.h"
// #include <qd/qd_real.h>
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
typedef double T;
typedef double T1;
#elif(the_p == 1)
typedef gdd_real T;
typedef dd_real T1;
#else
typedef gqd_real T;
typedef qd_real T1;
#endif

// typedef dd_real T1;


void print(complex<T>* a, complexH<T1>* b, int dim, string var);
void print(complex<T>* a, complexH<T1>** b, int dimX, int dimY, int stride, string var);

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int mode, int m,
   int ncoefs, char *pos_arr_h_char, char *exp_arr_h_char, complex<T> *x_h,
   complex<T> *c_h, complex<T> *factors_h, complex<T> *polvalues_h );

void CPU_evaldiff
 ( int dim, int NM, int NV, int deg, int r, int mode, int m,
   int *pos_arr_h_int, int *exp_arr_h_int, complexH<T1> *c,
   complexH<T1> *x, complex<T> *factors_h, complex<T> *polvalues_h );

void write_system ( int dim, 
int NM, 
int NV, 
complexH<T1> *c, 
int *myp, 
int *e );
// Writes the system of dimension dim, with number of monomials NM
// number of variables NV with coefficients in c, positions of nonzero
// exponents in p, and values of the exponents in e.

int main ( int argc, char *argv[] )
{

// complexH<T1> ExCompH; ExCompH.init(2.3,3.4); cout << "ExCompH=" << ExCompH;

  int BS,dim,NM,NV,deg,r,mode;
  if(parse_arguments(argc,argv,&BS,&dim,&NM,&NV,&deg,&r,&mode) == 1) return 1;

  int timevalue = 1287178355; // fixed seed instead of timevalue=time(NULL)
  const int n = 32;
  complex<T> *x = (complex<T>*)calloc(n,sizeof(complex<T>));
  complex<T> *y = (complex<T>*)calloc(n,sizeof(complex<T>));
  complex<T> *yD2 = (complex<T>*)calloc(n,sizeof(complex<T>));
  
   srand(timevalue);
   complex<T> *xp_h = new complex<T>[dim];
   complexH<T1> *xp = new complexH<T1>[dim];
   random_point(dim,xp_h,xp);
   // print(xp_h,xp,dim,"x");
   // for (int i=0;i<dim;i++) cout << "xp=" <<xp[i];
   storage<T1,4> pt; pt.copyToStor(xp);
   // pt.print(); 
   // cout << "M_PI=" << M_PI << endl;

   int pos_arr_h_int[NM*NV];
   char pos_arr_h_char[NM*NV];
   int exp_arr_h_int[NM*NV];
   char exp_arr_h_char[NM*NV];
   int ncoefs = NM*(NV+1);
   complex<T> *c_h = new complex<T>[ncoefs];
   complexH<T1> *c = new complexH<T1>[ncoefs];
   generate_system(dim,NM,NV,deg,pos_arr_h_int,pos_arr_h_char,
                   exp_arr_h_int,exp_arr_h_char,c_h,c);
   if(mode == 2) write_system(dim,NM,NV,c,pos_arr_h_int,exp_arr_h_int);
// allocate space for output
   int m = NM/dim;
   // cout << "m=" << m << endl;
   complex<T> *factors_h = new complex<T>[NM];
   int ant = ((dim*dim+dim)/BS + 1)*BS;
   // cout << "ant=" << ant << endl;
   complex<T> *polvalues_h = new complex<T>[ant];
   if(mode == 0 || mode == 2)
      GPU_evaldiff(BS,dim,NM,NV,deg,r,mode,m,ncoefs,pos_arr_h_char,
                   exp_arr_h_char,xp_h,c_h,factors_h,polvalues_h);

   if(mode == 1 || mode == 2)
      CPU_evaldiff(dim,NM,NV,deg,r,mode,m,pos_arr_h_int,exp_arr_h_int,
                   c,xp,factors_h,polvalues_h);

  for(int i = 0; i<n; i++) x[i].initH(double(i+2),0.12345);      

  int fail; // = sqrt_by_Newton(n,x,y);
  int fail1; //= squareCall(n,y,yD2);

  complexH<T1>* xH; complexH<T1>* yH; complexH<T1>* yD2H;
  xH= new complexH<T1>[n]; yH= new complexH<T1>[n]; yD2H= new complexH<T1>[n];

   comp1_gqdArr2qdArr(x, xH, n);
   comp1_gqdArr2qdArr(y, yH, n);
   comp1_gqdArr2qdArr(yD2, yD2H, n);

/*
 if(fail == 0 && fail1==0)
    for(int i=0;i<n;i++)
      {
        cout << "X["<< i << "]: \t \t" << xH[i];
        cout << "sqrtX[" << i << "]: \t" <<  yH[i];
        cout << "(sqrtX)^2[" << i << "]: \t" << yD2H[i];
        cout << endl;
       }
*/
// cout << sizeof(double) << endl;

  return 0;
}

void degrees ( complexH<T1> **deg, complexH<T1> *bm, int dim, int maxd )
{
   for(int i=0; i<dim; i++)
   {
      deg[i][0].init(1.0,0.0);
      for(int j=1; j<=maxd; j++)
      {  
         deg[i][j] = deg[i][j-1]*bm[i];
         // cout << "deg="<<deg[i][j];
      }
   }
}

void comp_factors ( int na, int no_monomials, int nvarm,
   int **monvarindexes, int **nonz_exponentsx,
   complexH<T1> *xval, complexH<T1> **deg, complexH<T1> *roots )
{
   for(int i=0; i<no_monomials; i++)
   {
      roots[i].init(1.0,0.0);
      // cout << "seq.roots=" << roots[i];
      for(int j=0; j<nvarm; j++)
      {
         roots[i]=roots[i]*deg[monvarindexes[i][j]][nonz_exponentsx[i][j]];
         // cout << deg[monvarindexes[i][j]][nonz_exponentsx[i][j]];
      }
      // cout << "seq.roots=" << roots[i];
      
   }
}

void speeldif_h
 ( int na, int NM, int nvarm, int nders, int **monvarindexes,
   complexH<T1> *xval, complexH<T1> *roots, complexH<T1> *coefs,
   complexH<T1> *monvalues, complexH<T1> **monderivatives )
{
   complexH<T1> *products1;
   complexH<T1> *products2;
   products1 = new complexH<T1>[nvarm];
   products2 = new complexH<T1>[nvarm];
   complexH<T1> *top_mon_derivatives;
   top_mon_derivatives = new complexH<T1>[nvarm];
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

void sum_mons
 ( int dim, int m, int nvarm, int **monvarindexes,
   complexH<T1> *monvalues, complexH<T1> **monderivatives,
   complexH<T1> *polyvalues, complexH<T1> **polyderivatives )
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

void CPU_evaldiff
 ( int dim, int NM, int NV, int deg, int r, int mode, int m,
   int *pos_arr_h_int, int *exp_arr_h_int, complexH<T1> *c,
   complexH<T1> *x, complex<T> *factors_h, complex<T> *polvalues_h )
// The CPU is used to evaluate a system and its Jacobian matrix.
{
   complexH<T1> **vdegrees = new complexH<T1>*[dim];
   for(int i=0; i<dim; i++) vdegrees[i] = new complexH<T1>[deg+1];
   int **positions_h = new int*[NM];
   int **exponents_h = new int*[NM];
   for(int i=0; i<NM; i++)
   {
      positions_h[i] = new int[NV];
      exponents_h[i] = new int[NV];
   }
   for(int i=0; i<NM; i++)
      for(int j=0; j<NV; j++) positions_h[i][j] = pos_arr_h_int[NV*i+j];
   for(int i=0; i<NM; i++)
      for(int j=0; j<NV; j++) exponents_h[i][j] = exp_arr_h_int[NV*i+j];
   complexH<T1> *factors_s = new complexH<T1>[NM];
   complexH<T1> *monvalues = new complexH<T1>[NM];
   complexH<T1> **derivatives = new complexH<T1>*[NM];
   for(int i=0; i<NM; i++) derivatives[i] = new complexH<T1>[NV];
   int nders = NM*NV;
   complexH<T1> *polyvalues = new complexH<T1>[dim];
   complexH<T1> **polyderivatives = new complexH<T1>*[dim];
   for(int i=0; i<dim; i++) polyderivatives[i] = new complexH<T1>[dim];
   for(int j=0; j<r; j++)
   {
      degrees(vdegrees,x,dim,deg);
      comp_factors(dim,NM,NV,positions_h,exponents_h,x,vdegrees,factors_s);
      // for(int i=0;i<NM;i++)
      //      cout << "sec. factors after comp-g="<< factors_s[i];
      // print(factors_h, factors_s,NM,"factors");
      speeldif_h(dim,NM,NV,nders,positions_h,x,factors_s,c,
                  monvalues,derivatives);
      // for(int i=0;i<NM;i++)
      //     cout << "sec. factors after comp-g=" <<factors_s[i];
      sum_mons(dim,m,NV,positions_h,monvalues,derivatives,
               polyvalues,polyderivatives);
      // print(polvalues_h,polyderivatives,dim,dim,dim,"poly_derivatives"); 
   }
   if(mode==2) // compare results of GPU with CPU
   {
      error_on_factors(NM,factors_h,factors_s);
      error_on_derivatives(dim,polyderivatives,polvalues_h);
   }
}

void print(complex<T>* a, complexH<T1>* b, int dim, string var)
{
   complexH<T1> temp;

   for (int i=0;i<dim;i++)
   {
      comp1_gqd2qd(&a[i],&temp);
      cout << "GPU: " << var << "["<< i << "] = " << temp;
      cout << "CPU: " << var << "["<< i << "] = " << b[i];
   }
}

void print(complex<T>* a, complexH<T1>** b, int dimX, int dimY, int stride, string var)
{

complexH<T1> temp;

for (int i=0;i<dimX;i++)
 for (int j=0;j<dimY;j++)
  {
    comp1_gqd2qd(&a[stride+dimY*j+i],&temp);
    cout << "GPU: " << var << "["<< i << "]" << "["<< j << "] = " << temp;
    cout << "CPU: " << var << "["<< i << "]" << "["<< j << "] = " << b[i][j];
  }

}

void write_system ( int dim, int NM, int NV, complexH<T1> *c, int *myp, int *e )
{
   cout << "          dimension : " << dim << endl;
   cout << "number of monomials : " << NM << endl;
   cout << "number of variables : " << NV << endl;
   cout << "   the coefficients : " << endl;
   cout << scientific << setprecision(8);
   for(int i=0; i<NM; i++)
   {
      cout << "c[" << i << "] : ";
      cout << "(" <<  c[i] << ")";
      for(int j=0; j<NV; j++)
         cout << "*x[" << myp[NV*i+j] << "]^" << e[NV*i+j];
      cout << endl;
   }
}

