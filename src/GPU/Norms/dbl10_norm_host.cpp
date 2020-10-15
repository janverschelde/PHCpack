// The file dbl10_norm_host.cpp defines the code for the functions
// specified in dbl10_norm_host.h.

#include "dbl10_norm_host.h"
#include "deca_double_functions.h"

using namespace std;

void make_copy
 ( int dim, double *orgrtb, double *orgrix, double *orgrmi, double *orgrrg,
   double *orgrpk, double *orgltb, double *orglix, double *orglmi, 
   double *orglrg, double *orglpk, double *duprtb, double *duprix,
   double *duprmi, double *duprrg, double *duprpk, double *dupltb,
   double *duplix, double *duplmi, double *duplrg, double *duplpk )
{
   for(int i=0; i<dim; i++)
   {
      duprtb[i] = orgrtb[i];
      duprix[i] = orgrix[i];
      duprmi[i] = orgrmi[i];
      duprrg[i] = orgrrg[i];
      duprpk[i] = orgrpk[i];
      dupltb[i] = orgltb[i];
      duplix[i] = orglix[i];
      duplmi[i] = orglmi[i];
      duplrg[i] = orglrg[i];
      duplpk[i] = orglpk[i];
   }
}

void CPU_norm
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double *normrtb, double *normrix, double *normrmi,
   double *normrrg, double *normrpk, double *normltb, double *normlix,
   double *normlmi, double *normlrg, double *normlpk )
{
   double sumrtb = 0.0;
   double sumrix = 0.0;
   double sumrmi = 0.0;
   double sumrrg = 0.0;
   double sumrpk = 0.0;
   double sumltb = 0.0;
   double sumlix = 0.0;
   double sumlmi = 0.0;
   double sumlrg = 0.0;
   double sumlpk = 0.0;
   double prodrtb,prodrix,prodrmi,prodrrg,prodrpk;
   double prodltb,prodlix,prodlmi,prodlrg,prodlpk;

   for(int i=0; i<dim; i++) // sum = sum + v[i]*v[i]
   {
      daf_sqr(vrtb[i],vrix[i],vrmi[i],vrrg[i],vrpk[i],
              vltb[i],vlix[i],vlmi[i],vlrg[i],vlpk[i],
              &prodrtb,&prodrix,&prodrmi,&prodrrg,&prodrpk,
              &prodltb,&prodlix,&prodlmi,&prodlrg,&prodlpk);
      daf_inc(&sumrtb,&sumrix,&sumrmi,&sumrrg,&sumrpk,
              &sumltb,&sumlix,&sumlmi,&sumlrg,&sumlpk,
              prodrtb,prodrix,prodrmi,prodrrg,prodrpk,
              prodltb,prodlix,prodlmi,prodlrg,prodlpk);
   }
   daf_sqrt(sumrtb,sumrix,sumrmi,sumrrg,sumrpk,
            sumltb,sumlix,sumlmi,sumlrg,sumlpk,
            normrtb,normrix,normrmi,normrrg,normrpk,
            normltb,normlix,normlmi,normlrg,normlpk);
}

void CPU_normalize
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double normrtb, double normrix, double normrmi, double normrrg,
   double normrpk, double normltb, double normlix, double normlmi,
   double normlrg, double normlpk )
{
   for(int i=0; i<dim; i++)
      daf_div(vrtb[i],vrix[i],vrmi[i],vrrg[i],vrpk[i],
              vltb[i],vlix[i],vlmi[i],vlrg[i],vlpk[i],
              normrtb,normrix,normrmi,normrrg,normrpk,
              normltb,normlix,normlmi,normlrg,normlpk,
              &vrtb[i],&vrix[i],&vrmi[i],&vrrg[i],&vrpk[i],
              &vltb[i],&vlix[i],&vlmi[i],&vlrg[i],&vlpk[i]);
}
