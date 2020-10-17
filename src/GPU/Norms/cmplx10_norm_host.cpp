// The file cmplx10_norm_host.cpp defines the code for the functions
// specified in cmplx10_norm_host.h.

#include "deca_double_functions.h"
#include "cmplx10_norm_host.h"

void make_copy
 ( int dim,
   double *orgrertb, double *orgrerix, double *orgrermi, double *orgrerrg,
   double *orgrerpk, double *orgreltb, double *orgrelix, double *orgrelmi,
   double *orgrelrg, double *orgrelpk, double *orgimrtb, double *orgimrix,
   double *orgimrmi, double *orgimrrg, double *orgimrpk, double *orgimltb,
   double *orgimlix, double *orgimlmi, double *orgimlrg, double *orgimlpk,
   double *duprertb, double *duprerix, double *duprermi, double *duprerrg,
   double *duprerpk, double *dupreltb, double *duprelix, double *duprelmi,
   double *duprelrg, double *duprelpk, double *dupimrtb, double *dupimrix,
   double *dupimrmi, double *dupimrrg, double *dupimrpk, double *dupimltb,
   double *dupimlix, double *dupimlmi, double *dupimlrg, double *dupimlpk )
{
   for(int i=0; i<dim; i++)
   {
      duprertb[i] = orgrertb[i];
      duprerix[i] = orgrerix[i];
      duprermi[i] = orgrermi[i];
      duprerrg[i] = orgrerrg[i];
      duprerpk[i] = orgrerpk[i];
      dupreltb[i] = orgreltb[i];
      duprelix[i] = orgrelix[i];
      duprelmi[i] = orgrelmi[i];
      duprelrg[i] = orgrelrg[i];
      duprelpk[i] = orgrelpk[i];
      dupimrtb[i] = orgimrtb[i];
      dupimrix[i] = orgimrix[i];
      dupimrmi[i] = orgimrmi[i];
      dupimrrg[i] = orgimrrg[i];
      dupimrpk[i] = orgimrpk[i];
      dupimltb[i] = orgimltb[i];
      dupimlix[i] = orgimlix[i];
      dupimlmi[i] = orgimlmi[i];
      dupimlrg[i] = orgimlrg[i];
      dupimlpk[i] = orgimlpk[i];
   }
}

void CPU_norm
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk )
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

   for(int k=0; k<dim; k++) // sum = sum + vre[k]*vre[k] + vim[k]*vim[k]
   {
      daf_sqr(  vrertb[k],vrerix[k],vrermi[k],vrerrg[k],vrerpk[k],
                vreltb[k],vrelix[k],vrelmi[k],vrelrg[k],vrelpk[k],
              &prodrtb, &prodrix, &prodrmi, &prodrrg, &prodrpk,
              &prodltb, &prodlix, &prodlmi, &prodlrg, &prodlpk);
      daf_inc(&sumrtb,&sumrix,&sumrmi,&sumrrg,&sumrpk,
              &sumltb,&sumlix,&sumlmi,&sumlrg,&sumlpk,
              prodrtb,prodrix,prodrmi,prodrrg,prodrpk,
              prodltb,prodlix,prodlmi,prodlrg,prodlpk);
      daf_sqr(  vimrtb[k],vimrix[k],vimrmi[k],vimrrg[k],vimrpk[k],
                vimltb[k],vimlix[k],vimlmi[k],vimlrg[k],vimlpk[k],
              &prodrtb, &prodrix, &prodrmi, &prodrrg, &prodrpk,
              &prodltb, &prodlix, &prodlmi, &prodlrg, &prodlpk);
      daf_inc(&sumrtb,&sumrix,&sumrmi,&sumrrg,&sumrpk,
              &sumltb,&sumlix,&sumlmi,&sumlrg,&sumlpk,
              prodrtb,prodrix,prodrmi,prodrrg,prodrpk,
              prodltb,prodlix,prodlmi,prodlrg,prodlpk);
   }
   daf_sqrt( sumrtb, sumrix, sumrmi, sumrrg, sumrpk,
             sumltb, sumlix, sumlmi, sumlrg, sumlpk,
            normrtb,normrix,normrmi,normrrg,normrpk,
            normltb,normlix,normlmi,normlrg,normlpk);
}

void CPU_normalize
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk, double *vimrtb, double *vimrix,
   double *vimrmi, double *vimrrg, double *vimrpk, double *vimltb,
   double *vimlix, double *vimlmi, double *vimlrg, double *vimlpk, int dim,
   double normrtb, double normrix, double normrmi, double normrrg,
   double normrpk, double normltb, double normlix, double normlmi,
   double normlrg, double normlpk )
{
   for(int i=0; i<dim; i++)
   {
      daf_div( vrertb[i], vrerix[i], vrermi[i], vrerrg[i], vrerpk[i],
               vreltb[i], vrelix[i], vrelmi[i], vrelrg[i], vrelpk[i],
              normrtb,   normrix,   normrmi,   normrrg,   normrpk,
              normltb,   normlix,   normlmi,   normlrg,   normlpk,
              &vrertb[i],&vrerix[i],&vrermi[i],&vrerrg[i],&vrerpk[i],
              &vreltb[i],&vrelix[i],&vrelmi[i],&vrelrg[i],&vrelpk[i]);
      daf_div( vimrtb[i], vimrix[i], vimrmi[i], vimrrg[i], vimrpk[i],
               vimltb[i], vimlix[i], vimlmi[i], vimlrg[i], vimlpk[i],
              normrtb,   normrix,   normrmi,   normrrg,   normrpk,
              normltb,   normlix,   normlmi,   normlrg,   normlpk,
              &vimrtb[i],&vimrix[i],&vimrmi[i],&vimrrg[i],&vimrpk[i],
              &vimltb[i],&vimlix[i],&vimlmi[i],&vimlrg[i],&vimlpk[i]);
   }
}
