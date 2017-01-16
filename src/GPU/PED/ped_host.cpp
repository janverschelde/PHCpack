// contains the definition for CPU_evaldiff in ped_host.h

#include "ped_host.h"

#define d  0 
#define dd 1
#define qd 2

#ifdef precision
#define the_p precision
#else
#define the_p 0
#endif

#if(the_p == 0)
typedef double realH;
#elif(the_p == 1)
typedef dd_real realH;
#else
typedef qd_real realH;
#endif

void degrees ( complexH<realH> **deg, complexH<realH> *bm, int dim, int maxd );

void comp_factors
 ( int na, int no_monomials, int nvarm,
   int **monvarindexes, int **nonz_exponentsx,
   complexH<realH> *xval, complexH<realH> **deg, complexH<realH> *roots );

void speeldif_h
 ( int na, int NM, int nvarm, int nders, int **monvarindexes,
   complexH<realH> *xval, complexH<realH> *roots, complexH<realH> *coefs,
   complexH<realH> *monvalues, complexH<realH> **monderivatives );

void sum_mons
 ( int dim, int m, int nvarm, int **monvarindexes,
   complexH<realH> *monvalues, complexH<realH> **monderivatives,
   complexH<realH> *polyvalues, complexH<realH> **polyderivatives );

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
      degrees(vdegrees,x,dim,deg);
      comp_factors(dim,NM,NV,positions_h,exponents_h,x,vdegrees,factors_h);
      speeldif_h
         (dim,NM,NV,nders,positions_h,x,factors_h,c,monvalues,derivatives);
      sum_mons
         (dim,m,NV,positions_h,monvalues,derivatives,polyvalues,polvalues_h);
   }
}

void degrees ( complexH<realH> **deg, complexH<realH> *bm, int dim, int maxd )
{
   for(int i=0; i<dim; i++)
   {
      deg[i][0].init(1.0,0.0);
      for(int j=1; j<=maxd; j++) deg[i][j] = deg[i][j-1]*bm[i];
   }
}

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
