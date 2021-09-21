/* The file dbl8_baqr_host.cpp defines the functions specified in
 * the file dbl8_baqr_host.h. */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "octo_double_functions.h"
#include "dbl8_factorizations.h"
#include "dbl8_baqr_host.h"

using namespace std;

void CPU_dbl8_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double **Vhihihi, double **Vlohihi, double **Vhilohi, double **Vlolohi,
   double **Vhihilo, double **Vlohilo, double **Vhilolo, double **Vlololo,
   double **Whihihi, double **Wlohihi, double **Whilohi, double **Wlolohi,
   double **Whihilo, double **Wlohilo, double **Whilolo, double **Wlololo )
{
   double *vhihihi = new double[nrows];
   double *vlohihi = new double[nrows];
   double *vhilohi = new double[nrows];
   double *vlolohi = new double[nrows];
   double *vhihilo = new double[nrows];
   double *vlohilo = new double[nrows];
   double *vhilolo = new double[nrows];
   double *vlololo = new double[nrows];
   double *phihihi = new double[ncols];
   double *plohihi = new double[ncols];
   double *philohi = new double[ncols];
   double *plolohi = new double[ncols];
   double *phihilo = new double[ncols];
   double *plohilo = new double[ncols];
   double *philolo = new double[ncols];
   double *plololo = new double[ncols];
   double zihihihi,zilohihi,zihilohi,zilolohi;
   double zihihilo,zilohilo,zihilolo,zilololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<nrows; i++)           //  W[0][i] = -B[0]*V[0][i];
   {
      odf_mul(Bhihihi[0],    Blohihi[0],    Bhilohi[0],    Blolohi[0],
              Bhihilo[0],    Blohilo[0],    Bhilolo[0],    Blololo[0],
              Vhihihi[0][i], Vlohihi[0][i], Vhilohi[0][i], Vlolohi[0][i],
              Vhihilo[0][i], Vlohilo[0][i], Vhilolo[0][i], Vlololo[0][i],
             &Whihihi[0][i],&Wlohihi[0][i],&Whilohi[0][i],&Wlolohi[0][i],
             &Whihilo[0][i],&Wlohilo[0][i],&Whilolo[0][i],&Wlololo[0][i]);
      odf_minus(&Whihihi[0][i],&Wlohihi[0][i],&Whilohi[0][i],&Wlolohi[0][i],
                &Whihilo[0][i],&Wlohilo[0][i],&Whilolo[0][i],&Wlololo[0][i]);
   }

   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++)
      {
         vhihihi[i] = Vhihihi[j][i];
         vlohihi[i] = Vlohihi[j][i];
         vhilohi[i] = Vhilohi[j][i];
         vlolohi[i] = Vlolohi[j][i];
         vhihilo[i] = Vhihilo[j][i];
         vlohilo[i] = Vlohilo[j][i];
         vhilolo[i] = Vhilolo[j][i];
         vlololo[i] = Vlololo[j][i];
      }
      for(int k=0; k<j; k++)
      {
         phihihi[k] = 0.0;             // compute k-th component of Y^T*v
         plohihi[k] = 0.0;             // over all rows of k-th column of Y
         philohi[k] = 0.0;
         plolohi[k] = 0.0;
         phihilo[k] = 0.0; 
         plohilo[k] = 0.0;
         philolo[k] = 0.0;
         plololo[k] = 0.0;

         for(int i=0; i<nrows; i++)     // p[k] = p[k] + V[k][i]*v[i]
         {
            odf_mul(Vhihihi[k][i],Vlohihi[k][i],Vhilohi[k][i],Vlolohi[k][i],
                    Vhihilo[k][i],Vlohilo[k][i],Vhilolo[k][i],Vlololo[k][i],
                    vhihihi[i],   vlohihi[i],   vhilohi[i],   vlolohi[i],
                    vhihilo[i],   vlohilo[i],   vhilolo[i],   vlololo[i],
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&phihihi[k],&plohihi[k],&philohi[k],&plolohi[k],
                    &phihilo[k],&plohilo[k],&philolo[k],&plololo[k],
                   acchihihi,  acclohihi,  acchilohi,  acclolohi,
                   acchihilo,  acclohilo,  acchilolo,  acclololo);
         }
      }
      for(int i=0; i<nrows; i++)
      {
         zihihihi = 0.0;               // compute i-th component of W*p
         zilohihi = 0.0;               // over all rows of k-th column of W
         zihilohi = 0.0;
         zilolohi = 0.0;
         zihihilo = 0.0; 
         zilohilo = 0.0;
         zihilolo = 0.0;
         zilololo = 0.0;

         for(int k=0; k<j; k++)         // zi = zi + W[k][i]*p[k]
         {
            odf_mul(Whihihi[k][i],Wlohihi[k][i],Whilohi[k][i],Wlolohi[k][i],
                    Whihilo[k][i],Wlohilo[k][i],Whilolo[k][i],Wlololo[k][i],
                    phihihi[k],   plohihi[k],   philohi[k],   plolohi[k],
                    phihilo[k],   plohilo[k],   philolo[k],   plololo[k],
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&zihihihi,&zilohihi,&zihilohi,&zilolohi,
                    &zihihilo,&zilohilo,&zihilolo,&zilololo,
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
         odf_inc(&zihihihi, &zilohihi, &zihilohi, &zilolohi,
                 &zihihilo, &zilohilo, &zihilolo, &zilololo,
                   vhihihi[i],vlohihi[i],vhilohi[i],vlolohi[i],
                   vhihilo[i],vlohilo[i],vhilolo[i],vlololo[i]);
         // W[j][i] = -B[j]*zi;
         odf_mul(Bhihihi[j],    Blohihi[j],    Bhilohi[j],    Blolohi[j],
                 Bhihilo[j],    Blohilo[j],    Bhilolo[j],    Blololo[j],
                zihihihi,      zilohihi,      zihilohi,      zilolohi,
                zihihilo,      zilohilo,      zihilolo,      zilololo,
                &Whihihi[j][i],&Wlohihi[j][i],&Whilohi[j][i],&Wlolohi[j][i],
                &Whihilo[j][i],&Wlohilo[j][i],&Whilolo[j][i],&Wlololo[j][i]);
         odf_minus(&Whihihi[j][i],&Wlohihi[j][i],&Whilohi[j][i],&Wlolohi[j][i],
                   &Whihilo[j][i],&Wlohilo[j][i],&Whilolo[j][i],&Wlololo[j][i]);
      }
   }
   free(vhihihi); free(vlohihi); free(vhilohi); free(vlolohi);
   free(vhihilo); free(vlohilo); free(vhilolo); free(vlololo);
   free(phihihi); free(plohihi); free(philohi); free(plolohi);
   free(phihilo); free(plohilo); free(philolo); free(plololo);
}

void CPU_cmplx8_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double **Vrehihihi, double **Vrelohihi,
   double **Vrehilohi, double **Vrelolohi,
   double **Vrehihilo, double **Vrelohilo,
   double **Vrehilolo, double **Vrelololo,
   double **Vimhihihi, double **Vimlohihi,
   double **Vimhilohi, double **Vimlolohi,
   double **Vimhihilo, double **Vimlohilo,
   double **Vimhilolo, double **Vimlololo,
   double **Wrehihihi, double **Wrelohihi,
   double **Wrehilohi, double **Wrelolohi,
   double **Wrehihilo, double **Wrelohilo,
   double **Wrehilolo, double **Wrelololo,
   double **Wimhihihi, double **Wimlohihi,
   double **Wimhilohi, double **Wimlolohi,
   double **Wimhihilo, double **Wimlohilo,
   double **Wimhilolo, double **Wimlololo )
{
   double *vrehihihi = new double[nrows];
   double *vrelohihi = new double[nrows];
   double *vrehilohi = new double[nrows];
   double *vrelolohi = new double[nrows];
   double *vrehihilo = new double[nrows];
   double *vrelohilo = new double[nrows];
   double *vrehilolo = new double[nrows];
   double *vrelololo = new double[nrows];
   double *vimhihihi = new double[nrows];
   double *vimlohihi = new double[nrows];
   double *vimhilohi = new double[nrows];
   double *vimlolohi = new double[nrows];
   double *vimhihilo = new double[nrows];
   double *vimlohilo = new double[nrows];
   double *vimhilolo = new double[nrows];
   double *vimlololo = new double[nrows];
   double *prehihihi = new double[ncols];
   double *prelohihi = new double[ncols];
   double *prehilohi = new double[ncols];
   double *prelolohi = new double[ncols];
   double *prehihilo = new double[ncols];
   double *prelohilo = new double[ncols];
   double *prehilolo = new double[ncols];
   double *prelololo = new double[ncols];
   double *pimhihihi = new double[ncols];
   double *pimlohihi = new double[ncols];
   double *pimhilohi = new double[ncols];
   double *pimlolohi = new double[ncols];
   double *pimhihilo = new double[ncols];
   double *pimlohilo = new double[ncols];
   double *pimhilolo = new double[ncols];
   double *pimlololo = new double[ncols];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double zi_rehihihi,zi_relohihi,zi_rehilohi,zi_relolohi;
   double zi_rehihilo,zi_relohilo,zi_rehilolo,zi_relololo;
   double zi_imhihihi,zi_imlohihi,zi_imhilohi,zi_imlolohi;
   double zi_imhihilo,zi_imlohilo,zi_imhilolo,zi_imlololo;

   for(int i=0; i<nrows; i++) // W[0][i] = -B[0]*V[0][i]
   {
      odf_mul(Bhihihi[0],      Blohihi[0],      Bhilohi[0],      Blolohi[0],
              Bhihilo[0],      Blohilo[0],      Bhilolo[0],      Blololo[0],
            Vrehihihi[0][i], Vrelohihi[0][i], Vrehilohi[0][i], Vrelolohi[0][i],
            Vrehihilo[0][i], Vrelohilo[0][i], Vrehilolo[0][i], Vrelololo[0][i],
           &Wrehihihi[0][i],&Wrelohihi[0][i],
           &Wrehilohi[0][i],&Wrelolohi[0][i],
           &Wrehihilo[0][i],&Wrelohilo[0][i],
           &Wrehilolo[0][i],&Wrelololo[0][i]);
      odf_minus(&Wrehihihi[0][i],&Wrelohihi[0][i],
                &Wrehilohi[0][i],&Wrelolohi[0][i],
                &Wrehihilo[0][i],&Wrelohilo[0][i],
                &Wrehilolo[0][i],&Wrelololo[0][i]);
      odf_mul(Bhihihi[0],      Blohihi[0],      Bhilohi[0],      Blolohi[0],
              Bhihilo[0],      Blohilo[0],      Bhilolo[0],      Blololo[0],
            Vimhihihi[0][i], Vimlohihi[0][i], Vimhilohi[0][i], Vimlolohi[0][i],
            Vimhihilo[0][i], Vimlohilo[0][i], Vimhilolo[0][i], Vimlololo[0][i],
           &Wimhihihi[0][i],&Wimlohihi[0][i],
           &Wimhilohi[0][i],&Wimlolohi[0][i],
           &Wimhihilo[0][i],&Wimlohilo[0][i],
           &Wimhilolo[0][i],&Wimlololo[0][i]);
      odf_minus(&Wimhihihi[0][i],&Wimlohihi[0][i],
                &Wimhilohi[0][i],&Wimlolohi[0][i],
                &Wimhihilo[0][i],&Wimlohilo[0][i],
                &Wimhilolo[0][i],&Wimlololo[0][i]);
   }
   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++)
      {
         vrehihihi[i] = Vrehihihi[j][i]; vrelohihi[i] = Vrelohihi[j][i];
         vrehilohi[i] = Vrehilohi[j][i]; vrelolohi[i] = Vrelolohi[j][i];
         vrehihilo[i] = Vrehihilo[j][i]; vrelohilo[i] = Vrelohilo[j][i];
         vrehilolo[i] = Vrehilolo[j][i]; vrelololo[i] = Vrelololo[j][i];
         vimhihihi[i] = Vimhihihi[j][i]; vimlohihi[i] = Vimlohihi[j][i];
         vimhilohi[i] = Vimhilohi[j][i]; vimlolohi[i] = Vimlolohi[j][i];
         vimhihilo[i] = Vimhihilo[j][i]; vimlohilo[i] = Vimlohilo[j][i];
         vimhilolo[i] = Vimhilolo[j][i]; vimlololo[i] = Vimlololo[j][i];
      }
      for(int k=0; k<j; k++)
      {
         prehihihi[k] = 0.0;        // compute k-th component of Y^H*v
         prelohihi[k] = 0.0;
         prehilohi[k] = 0.0;
         prelolohi[k] = 0.0;
         prehihilo[k] = 0.0; 
         prelohilo[k] = 0.0;
         prehilolo[k] = 0.0;
         prelololo[k] = 0.0;
         pimhihihi[k] = 0.0;
         pimlohihi[k] = 0.0;
         pimhilohi[k] = 0.0;
         pimlolohi[k] = 0.0;
         pimhihilo[k] = 0.0;
         pimlohilo[k] = 0.0;
         pimhilolo[k] = 0.0;
         pimlololo[k] = 0.0;

         for(int i=0; i<nrows; i++)    // p[k] = p[k] + V[k][i]*v[i];
         {
            // accre =   Vre[k][i]*vre[i] + Vim[k][i]*vim[i];
            // pre[k] = pre[k] + accre;
            odf_mul(Vrehihihi[k][i],Vrelohihi[k][i],
                    Vrehilohi[k][i],Vrelolohi[k][i],
                    Vrehihilo[k][i],Vrelohilo[k][i],
                    Vrehilolo[k][i],Vrelololo[k][i],
                    vrehihihi[i],   vrelohihi[i], vrehilohi[i], vrelolohi[i],
                    vrehihilo[i],   vrelohilo[i], vrehilolo[i], vrelololo[i],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&prehihihi[k],&prelohihi[k],&prehilohi[k],&prelolohi[k],
                    &prehihilo[k],&prelohilo[k],&prehilolo[k],&prelololo[k],
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
            odf_mul(Vimhihihi[k][i],Vimlohihi[k][i],
                    Vimhilohi[k][i],Vimlolohi[k][i],
                    Vimhihilo[k][i],Vimlohilo[k][i],
                    Vimhilolo[k][i],Vimlololo[k][i],
                    vimhihihi[i],   vimlohihi[i], vimhilohi[i], vimlolohi[i],
                    vimhihilo[i],   vimlohilo[i], vimhilolo[i], vimlololo[i],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&prehihihi[k],&prelohihi[k],&prehilohi[k],&prelolohi[k],
                    &prehihilo[k],&prelohilo[k],&prehilolo[k],&prelololo[k],
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
            // accim = - Vim[k][i]*vre[i] + Vre[k][i]*vim[i];
            // pim[k] = pim[k] + accim;
            odf_mul(Vimhihihi[k][i],Vimlohihi[k][i],
                    Vimhilohi[k][i],Vimlolohi[k][i],
                    Vimhihilo[k][i],Vimlohilo[k][i],
                    Vimhilolo[k][i],Vimlololo[k][i],
                    vrehihihi[i],   vrelohihi[i], vrehilohi[i], vrelolohi[i],
                    vrehihilo[i],   vrelohilo[i], vrehilolo[i], vrelololo[i],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_dec(&pimhihihi[k],&pimlohihi[k],&pimhilohi[k],&pimlolohi[k],
                    &pimhihilo[k],&pimlohilo[k],&pimhilolo[k],&pimlololo[k],
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
            odf_mul(Vrehihihi[k][i],Vrelohihi[k][i],
                    Vrehilohi[k][i],Vrelolohi[k][i],
                    Vrehihilo[k][i],Vrelohilo[k][i],
                    Vrehilolo[k][i],Vrelololo[k][i],
                    vimhihihi[i],   vimlohihi[i], vimhilohi[i], vimlolohi[i],
                    vimhihilo[i],   vimlohilo[i], vimhilolo[i], vimlololo[i],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&pimhihihi[k],&pimlohihi[k],&pimhilohi[k],&pimlolohi[k],
                    &pimhihilo[k],&pimlohilo[k],&pimhilolo[k],&pimlololo[k],
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         }
      }
      for(int i=0; i<nrows; i++)
      {
         zi_rehihihi = 0.0;           // compute i-th component of W*p
         zi_relohihi = 0.0;
         zi_rehilohi = 0.0;
         zi_relolohi = 0.0;
         zi_rehihilo = 0.0; 
         zi_relohilo = 0.0;
         zi_rehilolo = 0.0;
         zi_relololo = 0.0;
         zi_imhihihi = 0.0;
         zi_imlohihi = 0.0;
         zi_imhilohi = 0.0;
         zi_imlolohi = 0.0;
         zi_imhihilo = 0.0;
         zi_imlohilo = 0.0;
         zi_imhilolo = 0.0;
         zi_imlololo = 0.0;
 
         for(int k=0; k<j; k++)         // zi = zi + W[k][i]*p[k];
         {
            // accre = Wre[k][i]*pre[k] - Wim[k][i]*pim[k];
            // zi_re = zi_re + accre;
            odf_mul(Wrehihihi[k][i],Wrelohihi[k][i],
                    Wrehilohi[k][i],Wrelolohi[k][i],
                    Wrehihilo[k][i],Wrelohilo[k][i],
                    Wrehilolo[k][i],Wrelololo[k][i],
                    prehihihi[k],   prelohihi[k], prehilohi[k], prelolohi[k],
                    prehihilo[k],   prelohilo[k], prehilolo[k], prelololo[k],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&zi_rehihihi,&zi_relohihi,&zi_rehilohi,&zi_relolohi,
                    &zi_rehihilo,&zi_relohilo,&zi_rehilolo,&zi_relololo,
                       acchihihi,   acclohihi,   acchilohi,   acclolohi,
                       acchihilo,   acclohilo,   acchilolo,   acclololo);
            odf_mul(Wimhihihi[k][i],Wimlohihi[k][i],
                    Wimhilohi[k][i],Wimlolohi[k][i],
                    Wimhihilo[k][i],Wimlohilo[k][i],
                    Wimhilolo[k][i],Wimlololo[k][i],
                    pimhihihi[k],   pimlohihi[k], pimhilohi[k], pimlolohi[k],
                    pimhihilo[k],   pimlohilo[k], pimhilolo[k], pimlololo[k],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_dec(&zi_rehihihi,&zi_relohihi,&zi_rehilohi,&zi_relolohi,
                    &zi_rehihilo,&zi_relohilo,&zi_rehilolo,&zi_relololo,
                       acchihihi,   acclohihi,   acchilohi,   acclolohi,
                       acchihilo,   acclohilo,   acchilolo,   acclololo);
            // accim = Wim[k][i]*pre[k] + Wre[k][i]*pim[k];
            // zi_im = zi_im + accim;
            odf_mul(Wimhihihi[k][i],Wimlohihi[k][i],
                    Wimhilohi[k][i],Wimlolohi[k][i],
                    Wimhihilo[k][i],Wimlohilo[k][i],
                    Wimhilolo[k][i],Wimlololo[k][i],
                    prehihihi[k],   prelohihi[k], prehilohi[k], prelolohi[k],
                    prehihilo[k],   prelohilo[k], prehilolo[k], prelololo[k],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&zi_imhihihi,&zi_imlohihi,&zi_imhilohi,&zi_imlolohi,
                    &zi_imhihilo,&zi_imlohilo,&zi_imhilolo,&zi_imlololo,
                       acchihihi,   acclohihi,   acchilohi,   acclolohi,
                       acchihilo,   acclohilo,   acchilolo,   acclololo);
            odf_mul(Wrehihihi[k][i],Wrelohihi[k][i],
                    Wrehilohi[k][i],Wrelolohi[k][i],
                    Wrehihilo[k][i],Wrelohilo[k][i],
                    Wrehilolo[k][i],Wrelololo[k][i],
                    pimhihihi[k],   pimlohihi[k], pimhilohi[k], pimlolohi[k],
                    pimhihilo[k],   pimlohilo[k], pimhilolo[k], pimlololo[k],
                   &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                   &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&zi_imhihihi,&zi_imlohihi,&zi_imhilohi,&zi_imlolohi,
                    &zi_imhihilo,&zi_imlohilo,&zi_imhilolo,&zi_imlololo,
                       acchihihi,   acclohihi,   acchilohi,   acclolohi,
                       acchihilo,   acclohilo,   acchilolo,   acclololo);
         }
         // zi_re = zi_re + vre[i];
         odf_inc(&zi_rehihihi,&zi_relohihi,&zi_rehilohi,&zi_relolohi,
                 &zi_rehihilo,&zi_relohilo,&zi_rehilolo,&zi_relololo,
                    vrehihihi[i],vrelohihi[i],vrehilohi[i],vrelolohi[i],
                    vrehihilo[i],vrelohilo[i],vrehilolo[i],vrelololo[i]);
         // zi_im = zi_im + vim[i];
         odf_inc(&zi_imhihihi,&zi_imlohihi,&zi_imhilohi,&zi_imlolohi,
                 &zi_imhihilo,&zi_imlohilo,&zi_imhilolo,&zi_imlololo,
                    vimhihihi[i],vimlohihi[i],vimhilohi[i],vimlolohi[i],
                    vimhihilo[i],vimlohilo[i],vimhilolo[i],vimlololo[i]);
         // Wre[j][i] = -B[j]*zi_re;
         odf_mul(Bhihihi[j],      Blohihi[j],    Bhilohi[j],    Blolohi[j],
                 Bhihilo[j],      Blohilo[j],    Bhilolo[j],    Blololo[j],
             zi_rehihihi,     zi_relohihi,   zi_rehilohi,   zi_relolohi,
             zi_rehihilo,     zi_relohilo,   zi_rehilolo,   zi_relololo,
              &Wrehihihi[j][i],&Wrelohihi[j][i],
              &Wrehilohi[j][i],&Wrelolohi[j][i],
              &Wrehihilo[j][i],&Wrelohilo[j][i],
              &Wrehilolo[j][i],&Wrelololo[j][i]);
         odf_minus(&Wrehihihi[j][i],&Wrelohihi[j][i],
                   &Wrehilohi[j][i],&Wrelolohi[j][i],
                   &Wrehihilo[j][i],&Wrelohilo[j][i],
                   &Wrehilolo[j][i],&Wrelololo[j][i]);
         // Wim[j][i] = -B[j]*zi_im;
         odf_mul(Bhihihi[j],      Blohihi[j],    Bhilohi[j],    Blolohi[j],
                 Bhihilo[j],      Blohilo[j],    Bhilolo[j],    Blololo[j],
             zi_imhihihi,     zi_imlohihi,   zi_imhilohi,   zi_imlolohi,
             zi_imhihilo,     zi_imlohilo,   zi_imhilolo,   zi_imlololo,
              &Wimhihihi[j][i],&Wimlohihi[j][i],
              &Wimhilohi[j][i],&Wimlolohi[j][i],
              &Wimhihilo[j][i],&Wimlohilo[j][i],
              &Wimhilolo[j][i],&Wimlololo[j][i]);
         odf_minus(&Wimhihihi[j][i],&Wimlohihi[j][i],
                   &Wimhilohi[j][i],&Wimlolohi[j][i],
                   &Wimhihilo[j][i],&Wimlohilo[j][i],
                   &Wimhilolo[j][i],&Wimlololo[j][i]);
      }
   }
   free(vrehihihi); free(vrelohihi); free(vrehilohi); free(vrelolohi);
   free(vrehihilo); free(vrelohilo); free(vrehilolo); free(vrelololo);
   free(vimhihihi); free(vimlohihi); free(vimhilohi); free(vimlolohi);
   free(vimhihilo); free(vimlohilo); free(vimhilolo); free(vimlololo);
   free(prehihihi); free(prelohihi); free(prehilohi); free(prelolohi);
   free(prehihilo); free(prelohilo); free(prehilolo); free(prelololo);
   free(pimhihihi); free(pimlohihi); free(pimhilohi); free(pimlolohi);
   free(pimhihilo); free(pimlohilo); free(pimhilolo); free(pimlololo);
}

void CPU_dbl8_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Chihihi, double **Clohihi, double **Chilohi, double **Clolohi,
   double **Chihilo, double **Clohilo, double **Chilolo, double **Clololo,
   double **Yhihihi, double **Ylohihi, double **Yhilohi, double **Ylolohi,
   double **Yhihilo, double **Ylohilo, double **Yhilolo, double **Ylololo,
   double **Whihihi, double **Wlohihi, double **Whilohi, double **Wlolohi,
   double **Whihilo, double **Wlohilo, double **Whilolo, double **Wlololo,
   bool verbose )
{
   const int rowoff = idx*szt;            // row offset for C
   const int rowdim = nrows - rowoff;     // number of rows in Y and W
   const int coloff = (idx+1)*szt;        // column offset for C
   const int coldim = ncols - coloff;     // number of columns in C

   if(verbose)
   {
      cout << "updating R ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff << endl;
   }
   double **YWThihihi = new double*[rowdim];
   double **YWTlohihi = new double*[rowdim];
   double **YWThilohi = new double*[rowdim];
   double **YWTlolohi = new double*[rowdim];
   double **YWThihilo = new double*[rowdim];
   double **YWTlohilo = new double*[rowdim];
   double **YWThilolo = new double*[rowdim];
   double **YWTlololo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      YWThihihi[i] = new double[rowdim];
      YWTlohihi[i] = new double[rowdim];
      YWThilohi[i] = new double[rowdim];
      YWTlolohi[i] = new double[rowdim];
      YWThihilo[i] = new double[rowdim];
      YWTlohilo[i] = new double[rowdim];
      YWThilolo[i] = new double[rowdim];
      YWTlololo[i] = new double[rowdim];
   }
   double **prdhihihi = new double*[rowdim];
   double **prdlohihi = new double*[rowdim];
   double **prdhilohi = new double*[rowdim];
   double **prdlolohi = new double*[rowdim];
   double **prdhihilo = new double*[rowdim];
   double **prdlohilo = new double*[rowdim];
   double **prdhilolo = new double*[rowdim];
   double **prdlololo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      prdhihihi[i] = new double[coldim];
      prdlohihi[i] = new double[coldim];
      prdhilohi[i] = new double[coldim];
      prdlolohi[i] = new double[coldim];
      prdhihilo[i] = new double[coldim];
      prdlohilo[i] = new double[coldim];
      prdhilolo[i] = new double[coldim];
      prdlololo[i] = new double[coldim];
   }
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<rowdim; i++)        // compute Y*W^T
      for(int j=0; j<rowdim; j++)
      {
         YWThihihi[i][j] = 0.0;    // row i of Y with column j of W^T
         YWTlohihi[i][j] = 0.0;
         YWThilohi[i][j] = 0.0;
         YWTlolohi[i][j] = 0.0;
         YWThihilo[i][j] = 0.0;
         YWTlohilo[i][j] = 0.0;
         YWThilolo[i][j] = 0.0;
         YWTlololo[i][j] = 0.0;

         for(int k=0; k<szt; k++)     // YWT[i][j] += Y[k][i]*W[k][j]
         {
            odf_mul(Yhihihi[k][i],Ylohihi[k][i],Yhilohi[k][i],Ylolohi[k][i],
                    Yhihilo[k][i],Ylohilo[k][i],Yhilolo[k][i],Ylololo[k][i],
                    Whihihi[k][j],Wlohihi[k][j],Whilohi[k][j],Wlolohi[k][j],
                    Whihilo[k][j],Wlohilo[k][j],Whilolo[k][j],Wlololo[k][j],
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&YWThihihi[i][j],&YWTlohihi[i][j],
                    &YWThilohi[i][j],&YWTlolohi[i][j],
                    &YWThihilo[i][j],&YWTlohilo[i][j],
                    &YWThilolo[i][j],&YWTlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWThihihi[i][j] << "  " << YWTlohihi[i][j] << endl
                 << "            "
                 << YWThilohi[i][j] << "  " << YWTlolohi[i][j] << endl
                 << "            "
                 << YWThihilo[i][j] << "  " << YWTlohilo[i][j] << endl
                 << "            "
                 << YWThilolo[i][j] << "  " << YWTlololo[i][j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         prdhihihi[i][j] = 0.0;
         prdlohihi[i][j] = 0.0;
         prdhilohi[i][j] = 0.0;
         prdlolohi[i][j] = 0.0; 
         prdhihilo[i][j] = 0.0;
         prdlohilo[i][j] = 0.0;
         prdhilolo[i][j] = 0.0;
         prdlololo[i][j] = 0.0; 
                                      // prd[i][j]
         for(int k=0; k<rowdim; k++)  // += YWT[i][k]*C[rowoff+k][coloff+j]
         {
            odf_mul(YWThihihi[i][k],YWTlohihi[i][k],
                    YWThilohi[i][k],YWTlolohi[i][k],
                    YWThihilo[i][k],YWTlohilo[i][k],
                    YWThilolo[i][k],YWTlololo[i][k],
                      Chihihi[rowoff+k][coloff+j],Clohihi[rowoff+k][coloff+j],
                      Chilohi[rowoff+k][coloff+j],Clolohi[rowoff+k][coloff+j],
                      Chihilo[rowoff+k][coloff+j],Clohilo[rowoff+k][coloff+j],
                      Chilolo[rowoff+k][coloff+j],Clololo[rowoff+k][coloff+j],
                   &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                   &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdhihihi[i][j],&prdlohihi[i][j],
                    &prdhilohi[i][j],&prdlolohi[i][j],
                    &prdhihilo[i][j],&prdlohilo[i][j],
                    &prdhilolo[i][j],&prdlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   for(int i=0; i<rowdim; i++)      // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)   // C[rowoff+i][coloff+j] += prd[i][j]
      {
         odf_inc(&Chihihi[rowoff+i][coloff+j],&Clohihi[rowoff+i][coloff+j],
                 &Chilohi[rowoff+i][coloff+j],&Clolohi[rowoff+i][coloff+j],
                 &Chihilo[rowoff+i][coloff+j],&Clohilo[rowoff+i][coloff+j],
                 &Chilolo[rowoff+i][coloff+j],&Clololo[rowoff+i][coloff+j],
                 prdhihihi[i][j],prdlohihi[i][j],
                 prdhilohi[i][j],prdlolohi[i][j],
                 prdhihilo[i][j],prdlohilo[i][j],
                 prdhilolo[i][j],prdlololo[i][j]);
      }

   for(int i=0; i<rowdim; i++)
   {
      free(YWThihihi[i]); free(YWTlohihi[i]);
      free(YWThilohi[i]); free(YWTlolohi[i]);
      free(YWThihilo[i]); free(YWTlohilo[i]);
      free(YWThilolo[i]); free(YWTlololo[i]);
      free(prdhihihi[i]); free(prdlohihi[i]);
      free(prdhilohi[i]); free(prdlolohi[i]);
      free(prdhihilo[i]); free(prdlohilo[i]);
      free(prdhilolo[i]); free(prdlololo[i]);
   }
   free(YWThihihi); free(YWTlohihi); free(YWThilohi); free(YWTlolohi);
   free(YWThihilo); free(YWTlohilo); free(YWThilolo); free(YWTlololo);
   free(prdhihihi); free(prdlohihi); free(prdhilohi); free(prdlolohi);
   free(prdhihilo); free(prdlohilo); free(prdhilolo); free(prdlololo);
}

void CPU_cmplx8_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Crehihihi, double **Crelohihi,
   double **Crehilohi, double **Crelolohi,
   double **Crehihilo, double **Crelohilo,
   double **Crehilolo, double **Crelololo,
   double **Cimhihihi, double **Cimlohihi,
   double **Cimhilohi, double **Cimlolohi,
   double **Cimhihilo, double **Cimlohilo,
   double **Cimhilolo, double **Cimlololo,
   double **Yrehihihi, double **Yrelohihi,
   double **Yrehilohi, double **Yrelolohi,
   double **Yrehihilo, double **Yrelohilo,
   double **Yrehilolo, double **Yrelololo,
   double **Yimhihihi, double **Yimlohihi,
   double **Yimhilohi, double **Yimlolohi,
   double **Yimhihilo, double **Yimlohilo,
   double **Yimhilolo, double **Yimlololo,
   double **Wrehihihi, double **Wrelohihi,
   double **Wrehilohi, double **Wrelolohi,
   double **Wrehihilo, double **Wrelohilo,
   double **Wrehilolo, double **Wrelololo,
   double **Wimhihihi, double **Wimlohihi,
   double **Wimhilohi, double **Wimlolohi,
   double **Wimhihilo, double **Wimlohilo,
   double **Wimhilolo, double **Wimlololo,
   bool verbose )
{
   const int rowoff = idx*szt;             // row offset for C
   const int rowdim = nrows - rowoff;      // number of rows in Y and W
   const int coloff = (idx+1)*szt;         // column offset for C
   const int coldim = ncols - coloff;      // number of columns in C
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   if(verbose)
   {
      cout << "updating R ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff << endl;
   }
   double **YWTrehihihi = new double*[rowdim];
   double **YWTrelohihi = new double*[rowdim];
   double **YWTrehilohi = new double*[rowdim];
   double **YWTrelolohi = new double*[rowdim];
   double **YWTrehihilo = new double*[rowdim];
   double **YWTrelohilo = new double*[rowdim];
   double **YWTrehilolo = new double*[rowdim];
   double **YWTrelololo = new double*[rowdim];
   double **YWTimhihihi = new double*[rowdim];
   double **YWTimlohihi = new double*[rowdim];
   double **YWTimhilohi = new double*[rowdim];
   double **YWTimlolohi = new double*[rowdim];
   double **YWTimhihilo = new double*[rowdim];
   double **YWTimlohilo = new double*[rowdim];
   double **YWTimhilolo = new double*[rowdim];
   double **YWTimlololo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      YWTrehihihi[i] = new double[rowdim];
      YWTrelohihi[i] = new double[rowdim];
      YWTrehilohi[i] = new double[rowdim];
      YWTrelolohi[i] = new double[rowdim];
      YWTrehihilo[i] = new double[rowdim];
      YWTrelohilo[i] = new double[rowdim];
      YWTrehilolo[i] = new double[rowdim];
      YWTrelololo[i] = new double[rowdim];
      YWTimhihihi[i] = new double[rowdim];
      YWTimlohihi[i] = new double[rowdim];
      YWTimhilohi[i] = new double[rowdim];
      YWTimlolohi[i] = new double[rowdim];
      YWTimhihilo[i] = new double[rowdim];
      YWTimlohilo[i] = new double[rowdim];
      YWTimhilolo[i] = new double[rowdim];
      YWTimlololo[i] = new double[rowdim];
   }
   double **prdrehihihi = new double*[rowdim];
   double **prdrelohihi = new double*[rowdim];
   double **prdrehilohi = new double*[rowdim];
   double **prdrelolohi = new double*[rowdim];
   double **prdrehihilo = new double*[rowdim];
   double **prdrelohilo = new double*[rowdim];
   double **prdrehilolo = new double*[rowdim];
   double **prdrelololo = new double*[rowdim];
   double **prdimhihihi = new double*[rowdim];
   double **prdimlohihi = new double*[rowdim];
   double **prdimhilohi = new double*[rowdim];
   double **prdimlolohi = new double*[rowdim];
   double **prdimhihilo = new double*[rowdim];
   double **prdimlohilo = new double*[rowdim];
   double **prdimhilolo = new double*[rowdim];
   double **prdimlololo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      prdrehihihi[i] = new double[coldim];
      prdrelohihi[i] = new double[coldim];
      prdrehilohi[i] = new double[coldim];
      prdrelolohi[i] = new double[coldim];
      prdrehihilo[i] = new double[coldim];
      prdrelohilo[i] = new double[coldim];
      prdrehilolo[i] = new double[coldim];
      prdrelololo[i] = new double[coldim];
      prdimhihihi[i] = new double[coldim];
      prdimlohihi[i] = new double[coldim];
      prdimhilohi[i] = new double[coldim];
      prdimlolohi[i] = new double[coldim];
      prdimhihilo[i] = new double[coldim];
      prdimlohilo[i] = new double[coldim];
      prdimhilolo[i] = new double[coldim];
      prdimlololo[i] = new double[coldim];
   }
   for(int i=0; i<rowdim; i++)      // compute Y*W^H
      for(int j=0; j<rowdim; j++)
      {
         YWTrehihihi[i][j] = 0.0;   // row i of Y with column j of W^H
         YWTrelohihi[i][j] = 0.0;
         YWTrehilohi[i][j] = 0.0;
         YWTrelolohi[i][j] = 0.0;
         YWTrehihilo[i][j] = 0.0; 
         YWTrelohilo[i][j] = 0.0;
         YWTrehilolo[i][j] = 0.0;
         YWTrelololo[i][j] = 0.0;
         YWTimhihihi[i][j] = 0.0;
         YWTimlohihi[i][j] = 0.0;
         YWTimhilohi[i][j] = 0.0;
         YWTimlolohi[i][j] = 0.0;
         YWTimhihilo[i][j] = 0.0;
         YWTimlohilo[i][j] = 0.0;
         YWTimhilolo[i][j] = 0.0;
         YWTimlololo[i][j] = 0.0;

         for(int k=0; k<szt; k++)   // YWT[i][j] = YWT[i][j] + Y[k][i]*W[k][j]
         {  
            // accre = Yre[k][i]*Wre[k][j] + Yim[k][i]*Wim[k][j];
            // YWTre[i][j] = YWTre[i][j] + accre;
            odf_mul(Yrehihihi[k][i],Yrelohihi[k][i],
                    Yrehilohi[k][i],Yrelolohi[k][i],
                    Yrehihilo[k][i],Yrelohilo[k][i],
                    Yrehilolo[k][i],Yrelololo[k][i],
                    Wrehihihi[k][j],Wrelohihi[k][j],
                    Wrehilohi[k][j],Wrelolohi[k][j],
                    Wrehihilo[k][j],Wrelohilo[k][j],
                    Wrehilolo[k][j],Wrelololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&YWTrehihihi[i][j],&YWTrelohihi[i][j],
                    &YWTrehilohi[i][j],&YWTrelolohi[i][j],
                    &YWTrehihilo[i][j],&YWTrelohilo[i][j],
                    &YWTrehilolo[i][j],&YWTrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(Yimhihihi[k][i],Yimlohihi[k][i],
                    Yimhilohi[k][i],Yimlolohi[k][i],
                    Yimhihilo[k][i],Yimlohilo[k][i],
                    Yimhilolo[k][i],Yimlololo[k][i],
                    Wimhihihi[k][j],Wimlohihi[k][j],
                    Wimhilohi[k][j],Wimlolohi[k][j],
                    Wimhihilo[k][j],Wimlohilo[k][j],
                    Wimhilolo[k][j],Wimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&YWTrehihihi[i][j],&YWTrelohihi[i][j],
                    &YWTrehilohi[i][j],&YWTrelolohi[i][j],
                    &YWTrehihilo[i][j],&YWTrelohilo[i][j],
                    &YWTrehilolo[i][j],&YWTrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            // accim = Yim[k][i]*Wre[k][j] - Yre[k][i]*Wim[k][j];
            // YWTim[i][j] = YWTim[i][j] + accim;
            odf_mul(Yimhihihi[k][i],Yimlohihi[k][i],
                    Yimhilohi[k][i],Yimlolohi[k][i],
                    Yimhihilo[k][i],Yimlohilo[k][i],
                    Yimhilolo[k][i],Yimlololo[k][i],
                    Wrehihihi[k][j],Wrelohihi[k][j],
                    Wrehilohi[k][j],Wrelolohi[k][j],
                    Wrehihilo[k][j],Wrelohilo[k][j],
                    Wrehilolo[k][j],Wrelololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&YWTimhihihi[i][j],&YWTimlohihi[i][j],
                    &YWTimhilohi[i][j],&YWTimlolohi[i][j],
                    &YWTimhihilo[i][j],&YWTimlohilo[i][j],
                    &YWTimhilolo[i][j],&YWTimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(Yrehihihi[k][i],Yrelohihi[k][i],
                    Yrehilohi[k][i],Yrelolohi[k][i],
                    Yrehihilo[k][i],Yrelohilo[k][i],
                    Yrehilolo[k][i],Yrelololo[k][i],
                    Wimhihihi[k][j],Wimlohihi[k][j],
                    Wimhilohi[k][j],Wimlolohi[k][j],
                    Wimhihilo[k][j],Wimlohilo[k][j],
                    Wimhilolo[k][j],Wimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_dec(&YWTimhihihi[i][j],&YWTimlohihi[i][j],
                    &YWTimhilohi[i][j],&YWTimlolohi[i][j],
                    &YWTimhihilo[i][j],&YWTimlohilo[i][j],
                    &YWTimhilolo[i][j],&YWTimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "YWT[" << i << "][" << j << "]re : "
                 << YWTrehihihi[i][j] << "  " << YWTrelohihi[i][j] << endl
                 << "              "
                 << YWTrehilohi[i][j] << "  " << YWTrelolohi[i][j] << endl
                 << "              "
                 << YWTrehihilo[i][j] << "  " << YWTrelohilo[i][j] << endl
                 << "              "
                 << YWTrehilolo[i][j] << "  " << YWTrelololo[i][j] << endl;
            cout << "YWT[" << i << "][" << j << "]im : "
                 << YWTimhihihi[i][j] << "  " << YWTimlohihi[i][j] << endl
                 << "              "
                 << YWTimhilohi[i][j] << "  " << YWTimlolohi[i][j] << endl
                 << "              "
                 << YWTimhihilo[i][j] << "  " << YWTimlohilo[i][j] << endl
                 << "              "
                 << YWTimhilolo[i][j] << "  " << YWTimlololo[i][j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^H)*C
      for(int j=0; j<coldim; j++)
      {
         prdrehihihi[i][j] = 0.0;
         prdrelohihi[i][j] = 0.0;
         prdrehilohi[i][j] = 0.0;
         prdrelolohi[i][j] = 0.0;
         prdrehihilo[i][j] = 0.0;
         prdrelohilo[i][j] = 0.0;
         prdrehilolo[i][j] = 0.0;
         prdrelololo[i][j] = 0.0;
         prdimhihihi[i][j] = 0.0;
         prdimlohihi[i][j] = 0.0;
         prdimhilohi[i][j] = 0.0;
         prdimlolohi[i][j] = 0.0;
         prdimhihilo[i][j] = 0.0;
         prdimlohilo[i][j] = 0.0;
         prdimhilolo[i][j] = 0.0;
         prdimlololo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++)
         {  // prd[i][j] = prd[i][j] + YWT[i][k]*C[rowoff+k][coloff+j];
            // accre = YWTre[i][k]*Cre[rowoff+k][coloff+j]
            //       - YWTim[i][k]*Cim[rowoff+k][coloff+j];
            // prdre[i][j] = prdre[i][j] + accre;
            odf_mul(YWTrehihihi[i][k],YWTrelohihi[i][k],
                    YWTrehilohi[i][k],YWTrelolohi[i][k],
                    YWTrehihilo[i][k],YWTrelohilo[i][k],
                    YWTrehilolo[i][k],YWTrelololo[i][k],
                    Crehihihi[rowoff+k][coloff+j],Crelohihi[rowoff+k][coloff+j],
                    Crehilohi[rowoff+k][coloff+j],Crelolohi[rowoff+k][coloff+j],
                    Crehihilo[rowoff+k][coloff+j],Crelohilo[rowoff+k][coloff+j],
                    Crehilolo[rowoff+k][coloff+j],Crelololo[rowoff+k][coloff+j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdrehihihi[i][j],&prdrelohihi[i][j],
                    &prdrehilohi[i][j],&prdrelolohi[i][j],
                    &prdrehihilo[i][j],&prdrelohilo[i][j],
                    &prdrehilolo[i][j],&prdrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(YWTimhihihi[i][k],YWTimlohihi[i][k],
                    YWTimhilohi[i][k],YWTimlolohi[i][k],
                    YWTimhihilo[i][k],YWTimlohilo[i][k],
                    YWTimhilolo[i][k],YWTimlololo[i][k],
                    Cimhihihi[rowoff+k][coloff+j],Cimlohihi[rowoff+k][coloff+j],
                    Cimhilohi[rowoff+k][coloff+j],Cimlolohi[rowoff+k][coloff+j],
                    Cimhihilo[rowoff+k][coloff+j],Cimlohilo[rowoff+k][coloff+j],
                    Cimhilolo[rowoff+k][coloff+j],Cimlololo[rowoff+k][coloff+j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_dec(&prdrehihihi[i][j],&prdrelohihi[i][j],
                    &prdrehilohi[i][j],&prdrelolohi[i][j],
                    &prdrehihilo[i][j],&prdrelohilo[i][j],
                    &prdrehilolo[i][j],&prdrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            // accim = YWTim[i][k]*Cre[rowoff+k][coloff+j]
            //       + YWTre[i][k]*Cim[rowoff+k][coloff+j];
            // prdim[i][j] = prdim[i][j] + accim;
            odf_mul(YWTimhihihi[i][k],YWTimlohihi[i][k],
                    YWTimhilohi[i][k],YWTimlolohi[i][k],
                    YWTimhihilo[i][k],YWTimlohilo[i][k],
                    YWTimhilolo[i][k],YWTimlololo[i][k],
                    Crehihihi[rowoff+k][coloff+j],Crelohihi[rowoff+k][coloff+j],
                    Crehilohi[rowoff+k][coloff+j],Crelolohi[rowoff+k][coloff+j],
                    Crehihilo[rowoff+k][coloff+j],Crelohilo[rowoff+k][coloff+j],
                    Crehilolo[rowoff+k][coloff+j],Crelololo[rowoff+k][coloff+j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdimhihihi[i][j],&prdimlohihi[i][j],
                    &prdimhilohi[i][j],&prdimlolohi[i][j],
                    &prdimhihilo[i][j],&prdimlohilo[i][j],
                    &prdimhilolo[i][j],&prdimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(YWTrehihihi[i][k],YWTrelohihi[i][k],
                    YWTrehilohi[i][k],YWTrelolohi[i][k],
                    YWTrehihilo[i][k],YWTrelohilo[i][k],
                    YWTrehilolo[i][k],YWTrelololo[i][k],
                    Cimhihihi[rowoff+k][coloff+j],Cimlohihi[rowoff+k][coloff+j],
                    Cimhilohi[rowoff+k][coloff+j],Cimlolohi[rowoff+k][coloff+j],
                    Cimhihilo[rowoff+k][coloff+j],Cimlohilo[rowoff+k][coloff+j],
                    Cimhilolo[rowoff+k][coloff+j],Cimlololo[rowoff+k][coloff+j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdimhihihi[i][j],&prdimlohihi[i][j],
                    &prdimhilohi[i][j],&prdimlolohi[i][j],
                    &prdimhihilo[i][j],&prdimlohilo[i][j],
                    &prdimhilolo[i][j],&prdimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   if(verbose)
   {
      cout << endl;
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
         {
            cout << "prd[" << i << "][" << j << "]re : "
                 << prdrehihihi[i][j] << "  " << prdrelohihi[i][j] << endl
                 << "              "
                 << prdrehilohi[i][j] << "  " << prdrelolohi[i][j] << endl
                 << "              "
                 << prdrehihilo[i][j] << "  " << prdrelohilo[i][j] << endl
                 << "              "
                 << prdrehilolo[i][j] << "  " << prdrelololo[i][j] << endl;
            cout << "prd[" << i << "][" << j << "]im : "
                 << prdimhihihi[i][j] << "  " << prdimlohihi[i][j] << endl
                 << "              "
                 << prdimhilohi[i][j] << "  " << prdimlolohi[i][j] << endl
                 << "              "
                 << prdimhihilo[i][j] << "  " << prdimlohilo[i][j] << endl
                 << "              "
                 << prdimhilolo[i][j] << "  " << prdimlololo[i][j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)        // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         odf_inc(&Crehihihi[rowoff+i][coloff+j],&Crelohihi[rowoff+i][coloff+j],
                 &Crehilohi[rowoff+i][coloff+j],&Crelolohi[rowoff+i][coloff+j],
                 &Crehihilo[rowoff+i][coloff+j],&Crelohilo[rowoff+i][coloff+j],
                 &Crehilolo[rowoff+i][coloff+j],&Crelololo[rowoff+i][coloff+j],
                 prdrehihihi[i][j],prdrelohihi[i][j],
                 prdrehilohi[i][j],prdrelolohi[i][j],
                 prdrehihilo[i][j],prdrelohilo[i][j],
                 prdrehilolo[i][j],prdrelololo[i][j]);
         odf_inc(&Cimhihihi[rowoff+i][coloff+j],&Cimlohihi[rowoff+i][coloff+j],
                 &Cimhilohi[rowoff+i][coloff+j],&Cimlolohi[rowoff+i][coloff+j],
                 &Cimhihilo[rowoff+i][coloff+j],&Cimlohilo[rowoff+i][coloff+j],
                 &Cimhilolo[rowoff+i][coloff+j],&Cimlololo[rowoff+i][coloff+j],
                 prdimhihihi[i][j],prdimlohihi[i][j],
                 prdimhilohi[i][j],prdimlolohi[i][j],
                 prdimhihilo[i][j],prdimlohilo[i][j],
                 prdimhilolo[i][j],prdimlololo[i][j]);
      }

   for(int i=0; i<rowdim; i++)
   {
      free(YWTrehihihi[i]); free(YWTrelohihi[i]);
      free(YWTrehilohi[i]); free(YWTrelolohi[i]);
      free(YWTrehihilo[i]); free(YWTrelohilo[i]);
      free(YWTrehilolo[i]); free(YWTrelololo[i]);
      free(YWTimhihihi[i]); free(YWTimlohihi[i]);
      free(YWTimhilohi[i]); free(YWTimlolohi[i]);
      free(YWTimhihilo[i]); free(YWTimlohilo[i]);
      free(YWTimhilolo[i]); free(YWTimlololo[i]);
      free(prdrehihihi[i]); free(prdrelohihi[i]);
      free(prdrehilohi[i]); free(prdrelolohi[i]);
      free(prdrehihilo[i]); free(prdrelohilo[i]);
      free(prdrehilolo[i]); free(prdrelololo[i]);
      free(prdimhihihi[i]); free(prdimlohihi[i]);
      free(prdimhilohi[i]); free(prdimlolohi[i]);
      free(prdimhihilo[i]); free(prdimlohilo[i]);
      free(prdimhilolo[i]); free(prdimlololo[i]);
   }
   free(prdrehihihi); free(prdrelohihi); free(prdrehilohi); free(prdrelolohi);
   free(prdrehihilo); free(prdrelohilo); free(prdrehilolo); free(prdrelololo);
   free(prdimhihihi); free(prdimlohihi); free(prdimhilohi); free(prdimlolohi);
   free(prdimhihilo); free(prdimlohilo); free(prdimhilolo); free(prdimlololo);
   free(YWTrehihihi); free(YWTrelohihi); free(YWTrehilohi); free(YWTrelolohi);
   free(YWTrehihilo); free(YWTrelohilo); free(YWTrehilolo); free(YWTrelololo);
   free(YWTimhihihi); free(YWTimlohihi); free(YWTimhilohi); free(YWTimlolohi);
   free(YWTimhihilo); free(YWTimlohilo); free(YWTimhilolo); free(YWTimlololo);
}

void CPU_dbl8_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Yhihihi, double **Ylohihi, double **Yhilohi, double **Ylolohi,
   double **Yhihilo, double **Ylohilo, double **Yhilolo, double **Ylololo,
   double **Whihihi, double **Wlohihi, double **Whilohi, double **Wlolohi,
   double **Whihilo, double **Wlohilo, double **Whilolo, double **Wlololo,
   bool verbose )
{
   const int coloff = idx*szt;        // column offset for Q
   const int rowdim = dim - coloff;
   // the number of rows in Y and W is the number of columns to update

   if(verbose)
   {
      cout << "updating Q ..." << endl;
      cout << "-> dim : " << dim << "  szt : " << szt << "  idx : " << idx
           << "  rowdim : " << rowdim << "  coloff : " << coloff << endl;
   }
   double **WYThihihi = new double*[rowdim];
   double **WYTlohihi = new double*[rowdim];
   double **WYThilohi = new double*[rowdim];
   double **WYTlolohi = new double*[rowdim];
   double **WYThihilo = new double*[rowdim];
   double **WYTlohilo = new double*[rowdim];
   double **WYThilolo = new double*[rowdim];
   double **WYTlololo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      WYThihihi[i] = new double[rowdim];
      WYTlohihi[i] = new double[rowdim];
      WYThilohi[i] = new double[rowdim];
      WYTlolohi[i] = new double[rowdim];
      WYThihilo[i] = new double[rowdim];
      WYTlohilo[i] = new double[rowdim];
      WYThilolo[i] = new double[rowdim];
      WYTlololo[i] = new double[rowdim];
   }
   double **prdhihihi = new double*[dim];
   double **prdlohihi = new double*[dim];
   double **prdhilohi = new double*[dim];
   double **prdlolohi = new double*[dim];
   double **prdhihilo = new double*[dim];
   double **prdlohilo = new double*[dim];
   double **prdhilolo = new double*[dim];
   double **prdlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdhihihi[i] = new double[rowdim];
      prdlohihi[i] = new double[rowdim];
      prdhilohi[i] = new double[rowdim];
      prdlolohi[i] = new double[rowdim];
      prdhihilo[i] = new double[rowdim];
      prdlohilo[i] = new double[rowdim];
      prdhilolo[i] = new double[rowdim];
      prdlololo[i] = new double[rowdim];
   }
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<rowdim; i++)     // compute W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         WYThihihi[i][j] = 0.0;     // row i of W with column j of Y^T
         WYTlohihi[i][j] = 0.0;
         WYThilohi[i][j] = 0.0; 
         WYTlolohi[i][j] = 0.0; 
         WYThihilo[i][j] = 0.0; 
         WYTlohilo[i][j] = 0.0;
         WYThilolo[i][j] = 0.0; 
         WYTlololo[i][j] = 0.0;     // take k-th column of W

         for(int k=0; k<szt; k++)  // WYT[i][j] = WYT[i][j] + W[k][i]*Y[k][j]
         {
            odf_mul(Whihihi[k][i],Wlohihi[k][i],Whilohi[k][i],Wlolohi[k][i],
                    Whihilo[k][i],Wlohilo[k][i],Whilolo[k][i],Wlololo[k][i],
                    Yhihihi[k][j],Ylohihi[k][j],Yhilohi[k][j],Ylolohi[k][j],
                    Yhihilo[k][j],Ylohilo[k][j],Yhilolo[k][j],Ylololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&WYThihihi[i][j],&WYTlohihi[i][j],
                    &WYThilohi[i][j],&WYTlolohi[i][j],
                    &WYThihilo[i][j],&WYTlohilo[i][j],
                    &WYThilolo[i][j],&WYTlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   if(verbose)
   {
      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "Y[" << i << "][" << j << "] : "
                 << Yhihihi[j][i] << "  " << Ylohihi[j][i] << endl
                 << "          "
                 << Yhilohi[j][i] << "  " << Ylolohi[j][i] << endl
                 << "          "
                 << Yhihilo[j][i] << "  " << Ylohilo[j][i] << endl
                 << "          "
                 << Yhilolo[j][i] << "  " << Ylololo[j][i] << endl;

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "W[" << i << "][" << j << "] : "
                 << Whihihi[j][i] << "  " << Wlohihi[j][i] << endl
                 << "          "
                 << Whilohi[j][i] << "  " << Wlolohi[j][i] << endl
                 << "          "
                 << Whihilo[j][i] << "  " << Wlohilo[j][i] << endl
                 << "          "
                 << Whilolo[j][i] << "  " << Wlololo[j][i] << endl;

      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihihi[i][j] << "  " << WYTlohihi[i][j] << endl
                 << "            "
                 << WYThilohi[i][j] << "  " << WYTlolohi[i][j] << endl
                 << "            "
                 << WYThihilo[i][j] << "  " << WYTlohilo[i][j] << endl
                 << "            "
                 << WYThilolo[i][j] << "  " << WYTlololo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)          // prd = Q*W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         prdhihihi[i][j] = 0.0;
         prdlohihi[i][j] = 0.0;
         prdhilohi[i][j] = 0.0;
         prdlolohi[i][j] = 0.0;
         prdhihilo[i][j] = 0.0;
         prdlohilo[i][j] = 0.0;
         prdhilolo[i][j] = 0.0;
         prdlololo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++) // prd[i][j] += Q[i][coloff+k]*WYT[k][j]
         {
            odf_mul(Qhihihi[i][coloff+k],Qlohihi[i][coloff+k],
                    Qhilohi[i][coloff+k],Qlolohi[i][coloff+k],
                    Qhihilo[i][coloff+k],Qlohilo[i][coloff+k],
                    Qhilolo[i][coloff+k],Qlololo[i][coloff+k],
                    WYThihihi[k][j],WYTlohihi[k][j],
                    WYThilohi[k][j],WYTlolohi[k][j],
                    WYThihilo[k][j],WYTlohilo[k][j],
                    WYThilolo[k][j],WYTlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdhihihi[i][j],&prdlohihi[i][j],
                    &prdhilohi[i][j],&prdlolohi[i][j],
                    &prdhihilo[i][j],&prdlohilo[i][j],
                    &prdhilolo[i][j],&prdlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "QWYT[" << i << "][" << j << "] : "
                 << prdhihihi[i][j] << "  " << prdlohihi[i][j] << endl
                 << "             "
                 << prdhilohi[i][j] << "  " << prdlolohi[i][j] << endl
                 << "             "
                 << prdhihilo[i][j] << "  " << prdlohilo[i][j] << endl
                 << "             "
                 << prdhilolo[i][j] << "  " << prdlololo[i][j] << endl;

      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Qhihihi[i][coloff+j] << "  "
                 << Qlohihi[i][coloff+j] << endl
                 << "          "
                 << Qhilohi[i][coloff+j] << "  "
                 << Qlolohi[i][coloff+j] << endl
                 << "          "
                 << Qhihilo[i][coloff+j] << "  "
                 << Qlohilo[i][coloff+j] << endl
                 << "          "
                 << Qhilolo[i][coloff+j] << "  "
                 << Qlololo[i][coloff+j] << endl;
   }
   for(int i=0; i<dim; i++)        // Q = Q + Q*W*Y^T
      for(int j=0; j<rowdim; j++)  // Q[i][coloff+j] += prd[i][j];
      {
         odf_inc(&Qhihihi[i][coloff+j],&Qlohihi[i][coloff+j],
                 &Qhilohi[i][coloff+j],&Qlolohi[i][coloff+j],
                 &Qhihilo[i][coloff+j],&Qlohilo[i][coloff+j],
                 &Qhilolo[i][coloff+j],&Qlololo[i][coloff+j],
                 prdhihihi[i][j],prdlohihi[i][j],
                 prdhilohi[i][j],prdlolohi[i][j],
                 prdhihilo[i][j],prdlohilo[i][j],
                 prdhilolo[i][j],prdlololo[i][j]);
      }

   if(verbose)
   {
      cout << "Q after the update with QWYT :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Qhihihi[i][coloff+j] << "  "
                 << Qlohihi[i][coloff+j] << endl
                 << "          "
                 << Qhilohi[i][coloff+j] << "  "
                 << Qlolohi[i][coloff+j] << endl
                 << "          "
                 << Qhihilo[i][coloff+j] << "  "
                 << Qlohilo[i][coloff+j] << endl
                 << "          "
                 << Qhilolo[i][coloff+j] << "  "
                 << Qlololo[i][coloff+j] << endl;
   }
   for(int i=0; i<rowdim; i++)
   {
      free(WYThihihi[i]); free(WYTlohihi[i]);
      free(WYThilohi[i]); free(WYTlolohi[i]);
      free(WYThihilo[i]); free(WYTlohilo[i]);
      free(WYThilolo[i]); free(WYTlololo[i]);
   }
   for(int i=0; i<dim; i++)
   {
      free(prdhihihi[i]); free(prdlohihi[i]);
      free(prdhilohi[i]); free(prdlolohi[i]);
      free(prdhihilo[i]); free(prdlohilo[i]);
      free(prdhilolo[i]); free(prdlololo[i]);
   }
   free(WYThihihi); free(WYTlohihi); free(WYThilohi); free(WYTlolohi);
   free(WYThihilo); free(WYTlohilo); free(WYThilolo); free(WYTlololo);
   free(prdhihihi); free(prdlohihi); free(prdhilohi); free(prdlolohi);
   free(prdhihilo); free(prdlohilo); free(prdhilolo); free(prdlololo);
}

void CPU_cmplx8_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Yrehihihi, double **Yrelohihi,
   double **Yrehilohi, double **Yrelolohi,
   double **Yrehihilo, double **Yrelohilo,
   double **Yrehilolo, double **Yrelololo,
   double **Yimhihihi, double **Yimlohihi,
   double **Yimhilohi, double **Yimlolohi,
   double **Yimhihilo, double **Yimlohilo,
   double **Yimhilolo, double **Yimlololo,
   double **Wrehihihi, double **Wrelohihi,
   double **Wrehilohi, double **Wrelolohi,
   double **Wrehihilo, double **Wrelohilo,
   double **Wrehilolo, double **Wrelololo,
   double **Wimhihihi, double **Wimlohihi,
   double **Wimhilohi, double **Wimlolohi,
   double **Wimhihilo, double **Wimlohilo,
   double **Wimhilolo, double **Wimlololo, bool verbose )
{
   const int coloff = idx*szt;        // column offset for Q
   const int rowdim = dim - coloff;
   // the number of rows in Y and W is the number of columns to update
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   if(verbose)
   {
      cout << "updating Q ..." << endl;
      cout << "-> dim : " << dim << "  szt : " << szt << "  idx : " << idx
           << "  rowdim : " << rowdim << "  coloff : " << coloff << endl;
   }
   double **WYTrehihihi = new double*[rowdim];
   double **WYTrelohihi = new double*[rowdim];
   double **WYTrehilohi = new double*[rowdim];
   double **WYTrelolohi = new double*[rowdim];
   double **WYTrehihilo = new double*[rowdim];
   double **WYTrelohilo = new double*[rowdim];
   double **WYTrehilolo = new double*[rowdim];
   double **WYTrelololo = new double*[rowdim];
   double **WYTimhihihi = new double*[rowdim];
   double **WYTimlohihi = new double*[rowdim];
   double **WYTimhilohi = new double*[rowdim];
   double **WYTimlolohi = new double*[rowdim];
   double **WYTimhihilo = new double*[rowdim];
   double **WYTimlohilo = new double*[rowdim];
   double **WYTimhilolo = new double*[rowdim];
   double **WYTimlololo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      WYTrehihihi[i] = new double[rowdim];
      WYTrelohihi[i] = new double[rowdim];
      WYTrehilohi[i] = new double[rowdim];
      WYTrelolohi[i] = new double[rowdim];
      WYTrehihilo[i] = new double[rowdim];
      WYTrelohilo[i] = new double[rowdim];
      WYTrehilolo[i] = new double[rowdim];
      WYTrelololo[i] = new double[rowdim];
      WYTimhihihi[i] = new double[rowdim];
      WYTimlohihi[i] = new double[rowdim];
      WYTimhilohi[i] = new double[rowdim];
      WYTimlolohi[i] = new double[rowdim];
      WYTimhihilo[i] = new double[rowdim];
      WYTimlohilo[i] = new double[rowdim];
      WYTimhilolo[i] = new double[rowdim];
      WYTimlololo[i] = new double[rowdim];
   }
   double **prdrehihihi = new double*[dim];
   double **prdrelohihi = new double*[dim];
   double **prdrehilohi = new double*[dim];
   double **prdrelolohi = new double*[dim];
   double **prdrehihilo = new double*[dim];
   double **prdrelohilo = new double*[dim];
   double **prdrehilolo = new double*[dim];
   double **prdrelololo = new double*[dim];
   double **prdimhihihi = new double*[dim];
   double **prdimlohihi = new double*[dim];
   double **prdimhilohi = new double*[dim];
   double **prdimlolohi = new double*[dim];
   double **prdimhihilo = new double*[dim];
   double **prdimlohilo = new double*[dim];
   double **prdimhilolo = new double*[dim];
   double **prdimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrehihihi[i] = new double[rowdim];
      prdrelohihi[i] = new double[rowdim];
      prdrehilohi[i] = new double[rowdim];
      prdrelolohi[i] = new double[rowdim];
      prdrehihilo[i] = new double[rowdim];
      prdrelohilo[i] = new double[rowdim];
      prdrehilolo[i] = new double[rowdim];
      prdrelololo[i] = new double[rowdim];
      prdimhihihi[i] = new double[rowdim];
      prdimlohihi[i] = new double[rowdim];
      prdimhilohi[i] = new double[rowdim];
      prdimlolohi[i] = new double[rowdim];
      prdimhihilo[i] = new double[rowdim];
      prdimlohilo[i] = new double[rowdim];
      prdimhilolo[i] = new double[rowdim];
      prdimlololo[i] = new double[rowdim];
   }
   for(int i=0; i<rowdim; i++)        // compute W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         WYTrehihihi[i][j] = 0.0;     // row i of W with column j of Y^H
         WYTrelohihi[i][j] = 0.0;
         WYTrehilohi[i][j] = 0.0;
         WYTrelolohi[i][j] = 0.0;
         WYTrehihilo[i][j] = 0.0;
         WYTrelohilo[i][j] = 0.0;
         WYTrehilolo[i][j] = 0.0;
         WYTrelololo[i][j] = 0.0;
         WYTimhihihi[i][j] = 0.0;
         WYTimlohihi[i][j] = 0.0;
         WYTimhilohi[i][j] = 0.0;
         WYTimlolohi[i][j] = 0.0;
         WYTimhihilo[i][j] = 0.0;
         WYTimlohilo[i][j] = 0.0;
         WYTimhilolo[i][j] = 0.0;
         WYTimlololo[i][j] = 0.0;

         for(int k=0; k<szt; k++)     // WYT[i][j] += W[k][i]*Y[k][j]
         {
            // accre = Wre[k][i]*Yre[k][j] + Wim[k][i]*Yim[k][j];
            // WYTre[i][j] = WYTre[i][j] + accre;
            odf_mul(Wrehihihi[k][i],Wrelohihi[k][i],
                    Wrehilohi[k][i],Wrelolohi[k][i],
                    Wrehihilo[k][i],Wrelohilo[k][i],
                    Wrehilolo[k][i],Wrelololo[k][i],
                    Yrehihihi[k][j],Yrelohihi[k][j],
                    Yrehilohi[k][j],Yrelolohi[k][j],
                    Yrehihilo[k][j],Yrelohilo[k][j],
                    Yrehilolo[k][j],Yrelololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&WYTrehihihi[i][j],&WYTrelohihi[i][j],
                    &WYTrehilohi[i][j],&WYTrelolohi[i][j],
                    &WYTrehihilo[i][j],&WYTrelohilo[i][j],
                    &WYTrehilolo[i][j],&WYTrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(Wimhihihi[k][i],Wimlohihi[k][i],
                    Wimhilohi[k][i],Wimlolohi[k][i],
                    Wimhihilo[k][i],Wimlohilo[k][i],
                    Wimhilolo[k][i],Wimlololo[k][i],
                    Yimhihihi[k][j],Yimlohihi[k][j],
                    Yimhilohi[k][j],Yimlolohi[k][j],
                    Yimhihilo[k][j],Yimlohilo[k][j],
                    Yimhilolo[k][j],Yimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&WYTrehihihi[i][j],&WYTrelohihi[i][j],
                    &WYTrehilohi[i][j],&WYTrelolohi[i][j],
                    &WYTrehihilo[i][j],&WYTrelohilo[i][j],
                    &WYTrehilolo[i][j],&WYTrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            // accim = Wim[k][i]*Yre[k][j] - Wre[k][i]*Yim[k][j];
            // WYTim[i][j] = WYTim[i][j] + accim;
            odf_mul(Wimhihihi[k][i],Wimlohihi[k][i],
                    Wimhilohi[k][i],Wimlolohi[k][i],
                    Wimhihilo[k][i],Wimlohilo[k][i],
                    Wimhilolo[k][i],Wimlololo[k][i],
                    Yrehihihi[k][j],Yrelohihi[k][j],
                    Yrehilohi[k][j],Yrelolohi[k][j],
                    Yrehihilo[k][j],Yrelohilo[k][j],
                    Yrehilolo[k][j],Yrelololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&WYTimhihihi[i][j],&WYTimlohihi[i][j],
                    &WYTimhilohi[i][j],&WYTimlolohi[i][j],
                    &WYTimhihilo[i][j],&WYTimlohilo[i][j],
                    &WYTimhilolo[i][j],&WYTimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(Wrehihihi[k][i],Wrelohihi[k][i],
                    Wrehilohi[k][i],Wrelolohi[k][i],
                    Wrehihilo[k][i],Wrelohilo[k][i],
                    Wrehilolo[k][i],Wrelololo[k][i],
                    Yimhihihi[k][j],Yimlohihi[k][j],
                    Yimhilohi[k][j],Yimlolohi[k][j],
                    Yimhihilo[k][j],Yimlohilo[k][j],
                    Yimhilolo[k][j],Yimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_dec(&WYTimhihihi[i][j],&WYTimlohihi[i][j],
                    &WYTimhilohi[i][j],&WYTimlolohi[i][j],
                    &WYTimhihilo[i][j],&WYTimlohilo[i][j],
                    &WYTimhilolo[i][j],&WYTimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   if(verbose)
   {
      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
         {
            cout << "Y[" << i << "][" << j << "]re : "
                 << Yrehihihi[j][i] << "  " << Yrelohihi[j][i] << endl
                 << "            "
                 << Yrehilohi[j][i] << "  " << Yrelolohi[j][i] << endl
                 << "            "
                 << Yrehihilo[j][i] << "  " << Yrelohilo[j][i] << endl
                 << "            "
                 << Yrehilolo[j][i] << "  " << Yrelololo[j][i] << endl;
            cout << "Y[" << i << "][" << j << "]im : "
                 << Yimhihihi[j][i] << "  " << Yimlohihi[j][i] << endl
                 << "            "
                 << Yimhilohi[j][i] << "  " << Yimlolohi[j][i] << endl
                 << "            "
                 << Yimhihilo[j][i] << "  " << Yimlohilo[j][i] << endl
                 << "            "
                 << Yimhilolo[j][i] << "  " << Yimlololo[j][i] << endl;
         }

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
         {
            cout << "W[" << i << "][" << j << "]re : "
                 << Wrehihihi[j][i] << "  " << Wrelohihi[j][i] << endl
                 << "            "
                 << Wrehilohi[j][i] << "  " << Wrelolohi[j][i] << endl
                 << "            "
                 << Wrehihilo[j][i] << "  " << Wrelohilo[j][i] << endl
                 << "            "
                 << Wrehilolo[j][i] << "  " << Wrelololo[j][i] << endl;
            cout << "W[" << i << "][" << j << "]im : "
                 << Wimhihihi[j][i] << "  " << Wimlohihi[j][i] << endl
                 << "            "
                 << Wimhilohi[j][i] << "  " << Wimlolohi[j][i] << endl
                 << "            "
                 << Wimhihilo[j][i] << "  " << Wimlohilo[j][i] << endl
                 << "            "
                 << Wimhilolo[j][i] << "  " << Wimlololo[j][i] << endl;
         }

      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "WYT[" << i << "][" << j << "]re : "
                 << WYTrehihihi[i][j] << "  " << WYTrelohihi[i][j] << endl
                 << "            "
                 << WYTrehilohi[i][j] << "  " << WYTrelolohi[i][j] << endl
                 << "            "
                 << WYTrehihilo[i][j] << "  " << WYTrelohilo[i][j] << endl
                 << "            "
                 << WYTrehilolo[i][j] << "  " << WYTrelololo[i][j] << endl;
            cout << "WYT[" << i << "][" << j << "]im : "
                 << WYTimhihihi[i][j] << "  " << WYTimlohihi[i][j] << endl
                 << "            "
                 << WYTimhilohi[i][j] << "  " << WYTimlolohi[i][j] << endl
                 << "            "
                 << WYTimhihilo[i][j] << "  " << WYTimlohilo[i][j] << endl
                 << "            "
                 << WYTimhilolo[i][j] << "  " << WYTimlololo[i][j] << endl;
         }
   }
   for(int i=0; i<dim; i++)           // prd = Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         prdrehihihi[i][j] = 0.0;
         prdrelohihi[i][j] = 0.0;
         prdrehilohi[i][j] = 0.0;
         prdrelolohi[i][j] = 0.0;
         prdrehihilo[i][j] = 0.0;
         prdrelohilo[i][j] = 0.0;
         prdrehilolo[i][j] = 0.0;
         prdrelololo[i][j] = 0.0;
         prdimhihihi[i][j] = 0.0;
         prdimlohihi[i][j] = 0.0;
         prdimhilohi[i][j] = 0.0;
         prdimlolohi[i][j] = 0.0;
         prdimhihilo[i][j] = 0.0;
         prdimlohilo[i][j] = 0.0;
         prdimhilolo[i][j] = 0.0;
         prdimlololo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++) // prd[i][j] += Q[i][coloff+k]*WYT[k][j]
         {
            // accre = Qre[i][coloff+k]*WYTre[k][j]
            //       - Qim[i][coloff+k]*WYTim[k][j];
            // prdre[i][j] = prdre[i][j] + accre;
            odf_mul(Qrehihihi[i][coloff+k],Qrelohihi[i][coloff+k],
                    Qrehilohi[i][coloff+k],Qrelolohi[i][coloff+k],
                    Qrehihilo[i][coloff+k],Qrelohilo[i][coloff+k],
                    Qrehilolo[i][coloff+k],Qrelololo[i][coloff+k],
                    WYTrehihihi[k][j],WYTrelohihi[k][j],
                    WYTrehilohi[k][j],WYTrelolohi[k][j],
                    WYTrehihilo[k][j],WYTrelohilo[k][j],
                    WYTrehilolo[k][j],WYTrelololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdrehihihi[i][j],&prdrelohihi[i][j],
                    &prdrehilohi[i][j],&prdrelolohi[i][j],
                    &prdrehihilo[i][j],&prdrelohilo[i][j],
                    &prdrehilolo[i][j],&prdrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(Qimhihihi[i][coloff+k],Qimlohihi[i][coloff+k],
                    Qimhilohi[i][coloff+k],Qimlolohi[i][coloff+k],
                    Qimhihilo[i][coloff+k],Qimlohilo[i][coloff+k],
                    Qimhilolo[i][coloff+k],Qimlololo[i][coloff+k],
                    WYTimhihihi[k][j],WYTimlohihi[k][j],
                    WYTimhilohi[k][j],WYTimlolohi[k][j],
                    WYTimhihilo[k][j],WYTimlohilo[k][j],
                    WYTimhilolo[k][j],WYTimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_dec(&prdrehihihi[i][j],&prdrelohihi[i][j],
                    &prdrehilohi[i][j],&prdrelolohi[i][j],
                    &prdrehihilo[i][j],&prdrelohilo[i][j],
                    &prdrehilolo[i][j],&prdrelololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            // accim = Qim[i][coloff+k]*WYTre[k][j]
            //       + Qre[i][coloff+k]*WYTim[k][j];
            // prdim[i][j] = prdim[i][j] + accim;
            odf_mul(Qimhihihi[i][coloff+k],Qimlohihi[i][coloff+k],
                    Qimhilohi[i][coloff+k],Qimlolohi[i][coloff+k],
                    Qimhihilo[i][coloff+k],Qimlohilo[i][coloff+k],
                    Qimhilolo[i][coloff+k],Qimlololo[i][coloff+k],
                    WYTrehihihi[k][j],WYTrelohihi[k][j],
                    WYTrehilohi[k][j],WYTrelolohi[k][j],
                    WYTrehihilo[k][j],WYTrelohilo[k][j],
                    WYTrehilolo[k][j],WYTrelololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdimhihihi[i][j],&prdimlohihi[i][j],
                    &prdimhilohi[i][j],&prdimlolohi[i][j],
                    &prdimhihilo[i][j],&prdimlohilo[i][j],
                    &prdimhilolo[i][j],&prdimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
            odf_mul(Qrehihihi[i][coloff+k],Qrelohihi[i][coloff+k],
                    Qrehilohi[i][coloff+k],Qrelolohi[i][coloff+k],
                    Qrehihilo[i][coloff+k],Qrelohilo[i][coloff+k],
                    Qrehilolo[i][coloff+k],Qrelololo[i][coloff+k],
                    WYTimhihihi[k][j],WYTimlohihi[k][j],
                    WYTimhilohi[k][j],WYTimlolohi[k][j],
                    WYTimhihilo[k][j],WYTimlohilo[k][j],
                    WYTimhilolo[k][j],WYTimlololo[k][j],
                    &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                    &acchihilo,&acclohilo,&acchilolo,&acclololo);
            odf_inc(&prdimhihihi[i][j],&prdimlohihi[i][j],
                    &prdimhilohi[i][j],&prdimlolohi[i][j],
                    &prdimhihilo[i][j],&prdimlohilo[i][j],
                    &prdimhilolo[i][j],&prdimlololo[i][j],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }

   if(verbose)
   {
      cout << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "QWYT[" << i << "][" << j << "]re : "
                 << prdrehihihi[i][j] << "  " << prdrelohihi[i][j] << endl
                 << "               "
                 << prdrehilohi[i][j] << "  " << prdrelolohi[i][j] << endl
                 << "               "
                 << prdrehihilo[i][j] << "  " << prdrelohilo[i][j] << endl
                 << "               "
                 << prdrehilolo[i][j] << "  " << prdrelololo[i][j] << endl;
            cout << "QWYT[" << i << "][" << j << "]im : "
                 << prdimhihihi[i][j] << "  " << prdimlohihi[i][j] << endl
                 << "               "
                 << prdimhilohi[i][j] << "  " << prdimlolohi[i][j] << endl
                 << "               "
                 << prdimhihilo[i][j] << "  " << prdimlohilo[i][j] << endl
                 << "               "
                 << prdimhilolo[i][j] << "  " << prdimlololo[i][j] << endl;
         }

      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "Q[" << i << "][" << coloff+j << "]re : "
                 << Qrehihihi[i][coloff+j] << "  "
                 << Qrelohihi[i][coloff+j] << endl
                 << "            "
                 << Qrehilohi[i][coloff+j] << "  "
                 << Qrelolohi[i][coloff+j] << endl
                 << "            "
                 << Qrehihilo[i][coloff+j] << "  "
                 << Qrelohilo[i][coloff+j] << endl
                 << "            "
                 << Qrehilolo[i][coloff+j] << "  "
                 << Qrelololo[i][coloff+j] << endl;
            cout << "Q[" << i << "][" << coloff+j << "]im : "
                 << Qimhihihi[i][coloff+j] << "  "
                 << Qimlohihi[i][coloff+j] << endl
                 << "            "
                 << Qimhilohi[i][coloff+j] << "  "
                 << Qimlolohi[i][coloff+j] << endl
                 << "            "
                 << Qimhihilo[i][coloff+j] << "  "
                 << Qimlohilo[i][coloff+j] << endl
                 << "            "
                 << Qimhilolo[i][coloff+j] << "  "
                 << Qimlololo[i][coloff+j] << endl;
         }
   }
   for(int i=0; i<dim; i++)           // Q = Q + Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         // Qre[i][coloff+j] = Qre[i][coloff+j] + prdre[i][j];
         odf_inc(&Qrehihihi[i][coloff+j],&Qrelohihi[i][coloff+j],
                 &Qrehilohi[i][coloff+j],&Qrelolohi[i][coloff+j],
                 &Qrehihilo[i][coloff+j],&Qrelohilo[i][coloff+j],
                 &Qrehilolo[i][coloff+j],&Qrelololo[i][coloff+j],
                 prdrehihihi[i][j],prdrelohihi[i][j],
                 prdrehilohi[i][j],prdrelolohi[i][j],
                 prdrehihilo[i][j],prdrelohilo[i][j],
                 prdrehilolo[i][j],prdrelololo[i][j]);
         // Qim[i][coloff+j] = Qim[i][coloff+j] + prdim[i][j];
         odf_inc(&Qimhihihi[i][coloff+j],&Qimlohihi[i][coloff+j],
                 &Qimhilohi[i][coloff+j],&Qimlolohi[i][coloff+j],
                 &Qimhihilo[i][coloff+j],&Qimlohilo[i][coloff+j],
                 &Qimhilolo[i][coloff+j],&Qimlololo[i][coloff+j],
                 prdimhihihi[i][j],prdimlohihi[i][j],
                 prdimhilohi[i][j],prdimlolohi[i][j],
                 prdimhihilo[i][j],prdimlohilo[i][j],
                 prdimhilolo[i][j],prdimlololo[i][j]);
      }

   if(verbose)
   {
      cout << "Q after the update with QWYT :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "Q[" << i << "][" << coloff+j << "]re : "
                 << Qrehihihi[i][coloff+j] << "  "
                 << Qrelohihi[i][coloff+j] << endl
                 << "            "
                 << Qrehilohi[i][coloff+j] << "  "
                 << Qrelolohi[i][coloff+j] << endl
                 << "            "
                 << Qrehihilo[i][coloff+j] << "  "
                 << Qrelohilo[i][coloff+j] << endl
                 << "            "
                 << Qrehilolo[i][coloff+j] << "  "
                 << Qrelololo[i][coloff+j] << endl;
            cout << "Q[" << i << "][" << coloff+j << "]im : "
                 << Qimhihihi[i][coloff+j] << "  "
                 << Qimlohihi[i][coloff+j] << endl
                 << "            "
                 << Qimhilohi[i][coloff+j] << "  "
                 << Qimlolohi[i][coloff+j] << endl
                 << "            "
                 << Qimhihilo[i][coloff+j] << "  "
                 << Qimlohilo[i][coloff+j] << endl
                 << "            "
                 << Qimhilolo[i][coloff+j] << "  "
                 << Qimlololo[i][coloff+j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)
   {
      free(WYTrehihihi[i]); free(WYTrelohihi[i]);
      free(WYTrehilohi[i]); free(WYTrelolohi[i]);
      free(WYTrehihilo[i]); free(WYTrelohilo[i]);
      free(WYTrehilolo[i]); free(WYTrelololo[i]);
      free(WYTimhihihi[i]); free(WYTimlohihi[i]);
      free(WYTimhilohi[i]); free(WYTimlolohi[i]);
      free(WYTimhihilo[i]); free(WYTimlohilo[i]);
      free(WYTimhilolo[i]); free(WYTimlololo[i]);
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrehihihi[i]); free(prdrelohihi[i]);
      free(prdrehilohi[i]); free(prdrelolohi[i]);
      free(prdrehihilo[i]); free(prdrelohilo[i]);
      free(prdrehilolo[i]); free(prdrelololo[i]);
      free(prdimhihihi[i]); free(prdimlohihi[i]);
      free(prdimhilohi[i]); free(prdimlolohi[i]);
      free(prdimhihilo[i]); free(prdimlohilo[i]);
      free(prdimhilolo[i]); free(prdimlololo[i]);
   }
   free(WYTrehihihi); free(WYTrelohihi);
   free(WYTrehilohi); free(WYTrelolohi);
   free(WYTrehihilo); free(WYTrelohilo);
   free(WYTrehilolo); free(WYTrelololo);
   free(WYTimhihihi); free(WYTimlohihi);
   free(WYTimhilohi); free(WYTimlolohi);
   free(WYTimhihilo); free(WYTimlohilo);
   free(WYTimhilolo); free(WYTimlololo);
   free(prdrehihihi); free(prdrelohihi);
   free(prdrehilohi); free(prdrelolohi);
   free(prdrehihilo); free(prdrelohilo);
   free(prdrehilolo); free(prdrelololo);
   free(prdimhihihi); free(prdimlohihi);
   free(prdimhilohi); free(prdimlolohi);
   free(prdimhihilo); free(prdimlohilo);
   free(prdimhilolo); free(prdimlololo);
}

void CPU_dbl8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *lapsed, bool verbose )
{
   double betahihihi,betalohihi,betahilohi,betalolohi;
   double betahihilo,betalohilo,betahilolo,betalololo;
   double *xhihihi = new double[nrows]; // input vector for house
   double *xlohihi = new double[nrows];
   double *xhilohi = new double[nrows];
   double *xlolohi = new double[nrows];
   double *xhihilo = new double[nrows]; 
   double *xlohilo = new double[nrows];
   double *xhilolo = new double[nrows];
   double *xlololo = new double[nrows];
   double *vhihihi = new double[nrows]; // Householder vector
   double *vlohihi = new double[nrows];
   double *vhilohi = new double[nrows];
   double *vlolohi = new double[nrows];
   double *vhihilo = new double[nrows];
   double *vlohilo = new double[nrows];
   double *vhilolo = new double[nrows];
   double *vlololo = new double[nrows];
   double *Bhihihi = new double[szt];   // the betas
   double *Blohihi = new double[szt];
   double *Bhilohi = new double[szt];
   double *Blolohi = new double[szt];
   double *Bhihilo = new double[szt]; 
   double *Blohilo = new double[szt];
   double *Bhilolo = new double[szt];
   double *Blololo = new double[szt];
   double **Yhihihi = new double*[szt]; // Householder vectors in one block
   double **Ylohihi = new double*[szt];
   double **Yhilohi = new double*[szt];
   double **Ylolohi = new double*[szt];
   double **Yhihilo = new double*[szt]; 
   double **Ylohilo = new double*[szt];
   double **Yhilolo = new double*[szt];
   double **Ylololo = new double*[szt];
   double **Whihihi = new double*[szt]; // columns of W
   double **Wlohihi = new double*[szt];
   double **Whilohi = new double*[szt];
   double **Wlolohi = new double*[szt];
   double **Whihilo = new double*[szt];
   double **Wlohilo = new double*[szt];
   double **Whilolo = new double*[szt];
   double **Wlololo = new double*[szt];

   for(int j=0; j<szt; j++)
   {
      Yhihihi[j] = new double[nrows];
      Ylohihi[j] = new double[nrows];
      Yhilohi[j] = new double[nrows];
      Ylolohi[j] = new double[nrows];
      Yhihilo[j] = new double[nrows];
      Ylohilo[j] = new double[nrows];
      Yhilolo[j] = new double[nrows];
      Ylololo[j] = new double[nrows];
      Whihihi[j] = new double[nrows];
      Wlohihi[j] = new double[nrows];
      Whilohi[j] = new double[nrows];
      Wlolohi[j] = new double[nrows];
      Whihilo[j] = new double[nrows];
      Wlohilo[j] = new double[nrows];
      Whilolo[j] = new double[nrows];
      Wlololo[j] = new double[nrows];
   }
   clock_t start = clock();

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qhihihi[i][j] = 0.0;
         Qlohihi[i][j] = 0.0;
         Qhilohi[i][j] = 0.0;
         Qlolohi[i][j] = 0.0;
         Qhihilo[i][j] = 0.0;
         Qlohilo[i][j] = 0.0;
         Qhilolo[i][j] = 0.0;
         Qlololo[i][j] = 0.0;
      }
      Qhihihi[i][i] = 1.0;

      for(int j=0; j<ncols; j++) 
      {
         Rhihihi[i][j] = Ahihihi[i][j];
         Rlohihi[i][j] = Alohihi[i][j];
         Rhilohi[i][j] = Ahilohi[i][j];
         Rlolohi[i][j] = Alolohi[i][j];
         Rhihilo[i][j] = Ahihilo[i][j];
         Rlohilo[i][j] = Alohilo[i][j];
         Rhilolo[i][j] = Ahilolo[i][j];
         Rlololo[i][j] = Alololo[i][j];
      }
   }
   int colidx,endcol,nrowscol;

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      for(int L=0; L<szt; L++)    // L runs over the columns in one block
      {
         colidx = k*szt + L;      // index of the current column
         endcol = (k+1)*szt;      // 1 + last column index in current block
         nrowscol = nrows-colidx; // number of rows in current column

         for(int i=colidx; i<nrows; i++)
         {
            xhihihi[i-colidx] = Rhihihi[i][colidx];
            xlohihi[i-colidx] = Rlohihi[i][colidx];
            xhilohi[i-colidx] = Rhilohi[i][colidx];
            xlolohi[i-colidx] = Rlolohi[i][colidx];
            xhihilo[i-colidx] = Rhihilo[i][colidx];
            xlohilo[i-colidx] = Rlohilo[i][colidx];
            xhilolo[i-colidx] = Rhilolo[i][colidx];
            xlololo[i-colidx] = Rlololo[i][colidx];
         }
         CPU_dbl8_factors_house
            (nrowscol,xhihihi,    xlohihi,    xhilohi,    xlolohi,
                      xhihilo,    xlohilo,    xhilolo,    xlololo,
                      vhihihi,    vlohihi,    vhilohi,    vlolohi,
                      vhihilo,    vlohilo,    vhilolo,    vlololo,
                  &betahihihi,&betalohihi,&betahilohi,&betalolohi,
                  &betahihilo,&betalohilo,&betahilolo,&betalololo);

         if(verbose)
         {
            cout << "beta[" << colidx << "] : "
                 << betahihihi << "  " << betalohihi << endl
                 << "          "
                 << betahilohi << "  " << betalolohi << endl
                 << "          "
                 << betahihilo << "  " << betalohilo << endl
                 << "          "
                 << betahilolo << "  " << betalololo << endl;
            for(int i=colidx; i<nrows; i++)
               cout << "v[" << i-colidx << "] : "
                    << vhihihi[i-colidx] << "  " << vlohihi[i-colidx] << endl
                    << "       "
                    << vhilohi[i-colidx] << "  " << vlolohi[i-colidx] << endl
                    << "       "
                    << vhihilo[i-colidx] << "  " << vlohilo[i-colidx] << endl
                    << "       "
                    << vhilolo[i-colidx] << "  " << vlololo[i-colidx] << endl;
            cout << "the R matrix :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : "
                       << Rhihihi[i][j] << "  " << Rlohihi[i][j] << endl
                       << "          "
                       << Rhilohi[i][j] << "  " << Rlolohi[i][j] << endl
                       << "          "
                       << Rhihilo[i][j] << "  " << Rlohilo[i][j] << endl
                       << "          "
                       << Rhilolo[i][j] << "  " << Rlololo[i][j] << endl;
         }
         CPU_dbl8_factors_leftRupdate
            (nrows,endcol,colidx,Rhihihi,   Rlohihi,   Rhilohi,   Rlolohi,
                                 Rhihilo,   Rlohilo,   Rhilolo,   Rlololo,
                                 vhihihi,   vlohihi,   vhilohi,   vlolohi,
                                 vhihilo,   vlohilo,   vhilolo,   vlololo,
                              betahihihi,betalohihi,betahilohi,betalolohi,
                              betahihilo,betalohilo,betahilolo,betalololo);

         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : "
                       << Rhihihi[i][j] << "  " << Rlohihi[i][j] << endl
                       << "           "
                       << Rhilohi[i][j] << "  " << Rlolohi[i][j] << endl
                       << "           "
                       << Rhihilo[i][j] << "  " << Rlohilo[i][j] << endl
                       << "           "
                       << Rhilolo[i][j] << "  " << Rlololo[i][j] << endl;
         }
         Bhihihi[L] = betahihihi;
         Blohihi[L] = betalohihi;
         Bhilohi[L] = betahilohi;
         Blolohi[L] = betalolohi;
         Bhihilo[L] = betahihilo;
         Blohilo[L] = betalohilo;
         Bhilolo[L] = betahilolo;
         Blololo[L] = betalololo;

         for(int i=0; i<L; i++)
         {
            Yhihihi[L][i] = 0.0;
            Ylohihi[L][i] = 0.0;
            Yhilohi[L][i] = 0.0;
            Ylolohi[L][i] = 0.0;
            Yhihilo[L][i] = 0.0;
            Ylohilo[L][i] = 0.0;
            Yhilolo[L][i] = 0.0;
            Ylololo[L][i] = 0.0;
         }
         for(int i=0; i<nrowscol; i++)
         {
            Yhihihi[L][i+L] = vhihihi[i];
            Ylohihi[L][i+L] = vlohihi[i];
            Yhilohi[L][i+L] = vhilohi[i];
            Ylolohi[L][i+L] = vlolohi[i];
            Yhihilo[L][i+L] = vhihilo[i];
            Ylohilo[L][i+L] = vlohilo[i];
            Yhilolo[L][i+L] = vhilolo[i];
            Ylololo[L][i+L] = vlololo[i];
         }
      }
      CPU_dbl8_blocked_VB_to_W
         (nrows-k*szt,szt,Bhihihi,Blohihi,Bhilohi,Blolohi,
                          Bhihilo,Blohilo,Bhilolo,Blololo,
                          Yhihihi,Ylohihi,Yhilohi,Ylolohi,
                          Yhihilo,Ylohilo,Yhilolo,Ylololo,
                          Whihihi,Wlohihi,Whilohi,Wlolohi,
                          Whihilo,Wlohilo,Whilolo,Wlololo);

      if(k<nbt-1)
      {
         CPU_dbl8_blocked_leftRupdate
            (nrows,ncols,szt,k,Rhihihi,Rlohihi,Rhilohi,Rlolohi,
                               Rhihilo,Rlohilo,Rhilolo,Rlololo,
                               Yhihihi,Ylohihi,Yhilohi,Ylolohi,
                               Yhihilo,Ylohilo,Yhilolo,Ylololo,
                               Whihihi,Wlohihi,Whilohi,Wlolohi,
                               Whihilo,Wlohilo,Whilolo,Wlololo,verbose);
      }
      CPU_dbl8_blocked_rightQupdate
         (nrows,szt,k,Qhihihi,Qlohihi,Qhilohi,Qlolohi,
                      Qhihilo,Qlohilo,Qhilolo,Qlololo,
                      Yhihihi,Ylohihi,Yhilohi,Ylolohi,
                      Yhihilo,Ylohilo,Yhilolo,Ylololo,
                      Whihihi,Wlohihi,Whilohi,Wlolohi,
                      Whihilo,Wlohilo,Whilolo,Wlololo,verbose);
   }
   clock_t end = clock();
   *lapsed = double(end - start)/CLOCKS_PER_SEC;

   free(xhihihi); free(xlohihi); free(xhilohi); free(xlolohi);
   free(xhihilo); free(xlohilo); free(xhilolo); free(xlololo);
   free(vhihihi); free(vlohihi); free(vhilohi); free(vlolohi); 
   free(vhihilo); free(vlohilo); free(vhilolo); free(vlololo); 
   free(Bhihihi); free(Blohihi); free(Bhilohi); free(Blolohi);
   free(Bhihilo); free(Blohilo); free(Bhilolo); free(Blololo);

   for(int j=0; j<szt; j++)
   {
      free(Yhihihi[j]); free(Ylohihi[j]); free(Yhilohi[j]); free(Ylolohi[j]);
      free(Yhihilo[j]); free(Ylohilo[j]); free(Yhilolo[j]); free(Ylololo[j]);
      free(Whihihi[j]); free(Wlohihi[j]); free(Whilohi[j]); free(Wlolohi[j]);
      free(Whihilo[j]); free(Wlohilo[j]); free(Whilolo[j]); free(Wlololo[j]);
   }
   free(Yhihihi); free(Ylohihi); free(Yhilohi); free(Ylolohi);
   free(Yhihilo); free(Ylohilo); free(Yhilolo); free(Ylololo);
   free(Whihihi); free(Wlohihi); free(Whilohi); free(Wlolohi);
   free(Whihilo); free(Wlohilo); free(Whilolo); free(Wlololo);
}

void CPU_cmplx8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo, double *lapsec, bool verbose )
{
   double betahihihi,betalohihi,betahilohi,betalolohi;
   double betahihilo,betalohilo,betahilolo,betalololo;
   double *xrehihihi = new double[nrows]; // real input vector for house
   double *xrelohihi = new double[nrows];
   double *xrehilohi = new double[nrows];
   double *xrelolohi = new double[nrows];
   double *xrehihilo = new double[nrows]; 
   double *xrelohilo = new double[nrows];
   double *xrehilolo = new double[nrows];
   double *xrelololo = new double[nrows];
   double *ximhihihi = new double[nrows]; // imaginary input vector for house
   double *ximlohihi = new double[nrows];
   double *ximhilohi = new double[nrows];
   double *ximlolohi = new double[nrows];
   double *ximhihilo = new double[nrows]; 
   double *ximlohilo = new double[nrows];
   double *ximhilolo = new double[nrows];
   double *ximlololo = new double[nrows];
   double *vrehihihi = new double[nrows]; // real parts of a Householder vector
   double *vrelohihi = new double[nrows];
   double *vrehilohi = new double[nrows];
   double *vrelolohi = new double[nrows];
   double *vrehihilo = new double[nrows]; 
   double *vrelohilo = new double[nrows];
   double *vrehilolo = new double[nrows];
   double *vrelololo = new double[nrows];
   double *vimhihihi = new double[nrows]; // imag parts of a Householder vector
   double *vimlohihi = new double[nrows];
   double *vimhilohi = new double[nrows];
   double *vimlolohi = new double[nrows];
   double *vimhihilo = new double[nrows];
   double *vimlohilo = new double[nrows];
   double *vimhilolo = new double[nrows];
   double *vimlololo = new double[nrows];
   double *Bhihihi = new double[szt];     // the betas
   double *Blohihi = new double[szt];
   double *Bhilohi = new double[szt];
   double *Blolohi = new double[szt];
   double *Bhihilo = new double[szt]; 
   double *Blohilo = new double[szt];
   double *Bhilolo = new double[szt];
   double *Blololo = new double[szt];
   double **Yrehihihi = new double*[szt]; // Householder vectors in one block
   double **Yrelohihi = new double*[szt];
   double **Yrehilohi = new double*[szt];
   double **Yrelolohi = new double*[szt];
   double **Yrehihilo = new double*[szt];
   double **Yrelohilo = new double*[szt];
   double **Yrehilolo = new double*[szt];
   double **Yrelololo = new double*[szt];
   double **Yimhihihi = new double*[szt]; // imag parts of Householder vectors
   double **Yimlohihi = new double*[szt];
   double **Yimhilohi = new double*[szt];
   double **Yimlolohi = new double*[szt];
   double **Yimhihilo = new double*[szt];
   double **Yimlohilo = new double*[szt];
   double **Yimhilolo = new double*[szt];
   double **Yimlololo = new double*[szt];
   double **Wrehihihi = new double*[szt]; // real parts of the columns of W
   double **Wrelohihi = new double*[szt];
   double **Wrehilohi = new double*[szt]; 
   double **Wrelolohi = new double*[szt]; 
   double **Wrehihilo = new double*[szt];
   double **Wrelohilo = new double*[szt];
   double **Wrehilolo = new double*[szt]; 
   double **Wrelololo = new double*[szt]; 
   double **Wimhihihi = new double*[szt]; // imaginary parts of W
   double **Wimlohihi = new double*[szt];
   double **Wimhilohi = new double*[szt];
   double **Wimlolohi = new double*[szt];
   double **Wimhihilo = new double*[szt]; 
   double **Wimlohilo = new double*[szt];
   double **Wimhilolo = new double*[szt];
   double **Wimlololo = new double*[szt];

   for(int j=0; j<szt; j++)
   {
      Yrehihihi[j] = new double[nrows];
      Yrelohihi[j] = new double[nrows];
      Yrehilohi[j] = new double[nrows];
      Yrelolohi[j] = new double[nrows];
      Yrehihilo[j] = new double[nrows];
      Yrelohilo[j] = new double[nrows];
      Yrehilolo[j] = new double[nrows];
      Yrelololo[j] = new double[nrows];
      Yimhihihi[j] = new double[nrows];
      Yimlohihi[j] = new double[nrows];
      Yimhilohi[j] = new double[nrows];
      Yimlolohi[j] = new double[nrows];
      Yimhihilo[j] = new double[nrows];
      Yimlohilo[j] = new double[nrows];
      Yimhilolo[j] = new double[nrows];
      Yimlololo[j] = new double[nrows];
      Wrehihihi[j] = new double[nrows];
      Wrelohihi[j] = new double[nrows];
      Wrehilohi[j] = new double[nrows];
      Wrelolohi[j] = new double[nrows];
      Wrehihilo[j] = new double[nrows];
      Wrelohilo[j] = new double[nrows];
      Wrehilolo[j] = new double[nrows];
      Wrelololo[j] = new double[nrows];
      Wimhihihi[j] = new double[nrows];
      Wimlohihi[j] = new double[nrows];
      Wimhilohi[j] = new double[nrows];
      Wimlolohi[j] = new double[nrows];
      Wimhihilo[j] = new double[nrows];
      Wimlohilo[j] = new double[nrows];
      Wimhilolo[j] = new double[nrows];
      Wimlololo[j] = new double[nrows];
   }
   clock_t start = clock();

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qrehihihi[i][j] = 0.0; Qrelohihi[i][j] = 0.0;
         Qrehilohi[i][j] = 0.0; Qrelolohi[i][j] = 0.0;
         Qrehihilo[i][j] = 0.0; Qrelohilo[i][j] = 0.0;
         Qrehilolo[i][j] = 0.0; Qrelololo[i][j] = 0.0;
         Qimhihihi[i][j] = 0.0; Qimlohihi[i][j] = 0.0;
         Qimhilohi[i][j] = 0.0; Qimlolohi[i][j] = 0.0;
         Qimhihilo[i][j] = 0.0; Qimlohilo[i][j] = 0.0;
         Qimhilolo[i][j] = 0.0; Qimlololo[i][j] = 0.0;
      }
      Qrehihihi[i][i] = 1.0;

      for(int j=0; j<ncols; j++)
      {
         Rrehihihi[i][j] = Arehihihi[i][j]; Rrelohihi[i][j] = Arelohihi[i][j];
         Rrehilohi[i][j] = Arehilohi[i][j]; Rrelolohi[i][j] = Arelolohi[i][j];
         Rrehihilo[i][j] = Arehihilo[i][j]; Rrelohilo[i][j] = Arelohilo[i][j];
         Rrehilolo[i][j] = Arehilolo[i][j]; Rrelololo[i][j] = Arelololo[i][j];
         Rimhihihi[i][j] = Aimhihihi[i][j]; Rimlohihi[i][j] = Aimlohihi[i][j];
         Rimhilohi[i][j] = Aimhilohi[i][j]; Rimlolohi[i][j] = Aimlolohi[i][j];
         Rimhihilo[i][j] = Aimhihilo[i][j]; Rimlohilo[i][j] = Aimlohilo[i][j];
         Rimhilolo[i][j] = Aimhilolo[i][j]; Rimlololo[i][j] = Aimlololo[i][j];
      }
   }
   int colidx,endcol,nrowscol;

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      for(int L=0; L<szt; L++)    // L runs over the columns in one block
      {
         colidx = k*szt + L;      // index of the current column
         endcol = (k+1)*szt;      // 1 + last column index in current block
         nrowscol = nrows-colidx; // number of rows in current column

         for(int i=colidx; i<nrows; i++)
         {
            xrehihihi[i-colidx] = Rrehihihi[i][colidx];
            xrelohihi[i-colidx] = Rrelohihi[i][colidx];
            xrehilohi[i-colidx] = Rrehilohi[i][colidx];
            xrelolohi[i-colidx] = Rrelolohi[i][colidx];
            xrehihilo[i-colidx] = Rrehihilo[i][colidx];
            xrelohilo[i-colidx] = Rrelohilo[i][colidx];
            xrehilolo[i-colidx] = Rrehilolo[i][colidx];
            xrelololo[i-colidx] = Rrelololo[i][colidx];
            ximhihihi[i-colidx] = Rimhihihi[i][colidx];
            ximlohihi[i-colidx] = Rimlohihi[i][colidx];
            ximhilohi[i-colidx] = Rimhilohi[i][colidx];
            ximlolohi[i-colidx] = Rimlolohi[i][colidx];
            ximhihilo[i-colidx] = Rimhihilo[i][colidx];
            ximlohilo[i-colidx] = Rimlohilo[i][colidx];
            ximhilolo[i-colidx] = Rimhilolo[i][colidx];
            ximlololo[i-colidx] = Rimlololo[i][colidx];
         }
         CPU_cmplx8_factors_house
            (nrowscol,xrehihihi,  xrelohihi,  xrehilohi,  xrelolohi,
                      xrehihilo,  xrelohilo,  xrehilolo,  xrelololo,
                      ximhihihi,  ximlohihi,  ximhilohi,  ximlolohi,
                      ximhihilo,  ximlohilo,  ximhilolo,  ximlololo,
                      vrehihihi,  vrelohihi,  vrehilohi,  vrelolohi,
                      vrehihilo,  vrelohilo,  vrehilolo,  vrelololo,
                      vimhihihi,  vimlohihi,  vimhilohi,  vimlolohi,
                      vimhihilo,  vimlohilo,  vimhilolo,  vimlololo,
                    &betahihihi,&betalohihi,&betahilohi,&betalolohi,
                    &betahihilo,&betalohilo,&betahilolo,&betalololo);
         if(verbose)
         {
            cout << "beta[" << colidx << "] : "
                 << betahihihi << "  " << betalohihi << endl
                 << "          "
                 << betahilohi << "  " << betalolohi << endl
                 << "          "
                 << betahihilo << "  " << betalohilo << endl
                 << "          "
                 << betahilolo << "  " << betalololo << endl;

            for(int i=colidx; i<nrows; i++)
            {
               cout << "v[" << i-colidx << "]re : "
                    << vrehihihi[i-colidx] << "  "
                    << vrelohihi[i-colidx] << endl
                    << "         "
                    << vrehilohi[i-colidx] << "  "
                    << vrelolohi[i-colidx] << endl
                    << "         "
                    << vrehihilo[i-colidx] << "  "
                    << vrelohilo[i-colidx] << endl
                    << "         "
                    << vrehilolo[i-colidx] << "  "
                    << vrelololo[i-colidx] << endl;
               cout << "v[" << i-colidx << "]im : "
                    << vimhihihi[i-colidx] << "  "
                    << vimlohihi[i-colidx] << endl
                    << "         "
                    << vimhilohi[i-colidx] << "  "
                    << vimlolohi[i-colidx] << endl
                    << "         "
                    << vimhihilo[i-colidx] << "  "
                    << vimlohilo[i-colidx] << endl
                    << "         "
                    << vimhilolo[i-colidx] << "  "
                    << vimlololo[i-colidx] << endl;
            }
         }
         CPU_cmplx8_factors_leftRupdate
            (nrows,endcol,colidx,Rrehihihi, Rrelohihi, Rrehilohi, Rrelolohi,
                                 Rrehihilo, Rrelohilo, Rrehilolo, Rrelololo,
                                 Rimhihihi, Rimlohihi, Rimhilohi, Rimlolohi,
                                 Rimhihilo, Rimlohilo, Rimhilolo, Rimlololo,
                                 vrehihihi, vrelohihi, vrehilohi, vrelolohi,
                                 vrehihilo, vrelohilo, vrehilolo, vrelololo,
                                 vimhihihi, vimlohihi, vimhilohi, vimlolohi,
                                 vimhihilo, vimlohilo, vimhilolo, vimlololo,
                                betahihihi,betalohihi,betahilohi,betalolohi,
                                betahihilo,betalohilo,betahilolo,betalololo);
         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
               {
                  cout << "R[" << i << "][" << j << "]re : "
                       << Rrehihihi[i][j] << "  " << Rrelohihi[i][j] << endl
                       << "            "
                       << Rrehilohi[i][j] << "  " << Rrelolohi[i][j] << endl
                       << "            "
                       << Rrehihilo[i][j] << "  " << Rrelohilo[i][j] << endl
                       << "            "
                       << Rrehilolo[i][j] << "  " << Rrelololo[i][j] << endl;
                  cout << "R[" << i << "][" << j << "]im : "
                       << Rimhihihi[i][j] << "  " << Rimlohihi[i][j] << endl
                       << "            "
                       << Rimhilohi[i][j] << "  " << Rimlolohi[i][j] << endl
                       << "            "
                       << Rimhihilo[i][j] << "  " << Rimlohilo[i][j] << endl
                       << "            "
                       << Rimhilolo[i][j] << "  " << Rimlololo[i][j] << endl;
               }
         }
         Bhihihi[L] = betahihihi;
         Blohihi[L] = betalohihi;
         Bhilohi[L] = betahilohi;
         Blolohi[L] = betalolohi;
         Bhihilo[L] = betahihilo;
         Blohilo[L] = betalohilo;
         Bhilolo[L] = betahilolo;
         Blololo[L] = betalololo;

         for(int i=0; i<L; i++)
         {
            Yrehihihi[L][i] = 0.0;
            Yrelohihi[L][i] = 0.0;
            Yrehilohi[L][i] = 0.0;
            Yrelolohi[L][i] = 0.0;
            Yrehihilo[L][i] = 0.0;
            Yrelohilo[L][i] = 0.0;
            Yrehilolo[L][i] = 0.0;
            Yrelololo[L][i] = 0.0;
            Yimhihihi[L][i] = 0.0;
            Yimlohihi[L][i] = 0.0;
            Yimhilohi[L][i] = 0.0;
            Yimlolohi[L][i] = 0.0;
            Yimhihilo[L][i] = 0.0;
            Yimlohilo[L][i] = 0.0;
            Yimhilolo[L][i] = 0.0;
            Yimlololo[L][i] = 0.0;
         }
         for(int i=0; i<nrowscol; i++)
         {
            Yrehihihi[L][i+L] = vrehihihi[i];
            Yrelohihi[L][i+L] = vrelohihi[i];
            Yrehilohi[L][i+L] = vrehilohi[i];
            Yrelolohi[L][i+L] = vrelolohi[i];
            Yrehihilo[L][i+L] = vrehihilo[i];
            Yrelohilo[L][i+L] = vrelohilo[i];
            Yrehilolo[L][i+L] = vrehilolo[i];
            Yrelololo[L][i+L] = vrelololo[i];
            Yimhihihi[L][i+L] = vimhihihi[i];
            Yimlohihi[L][i+L] = vimlohihi[i];
            Yimhilohi[L][i+L] = vimhilohi[i];
            Yimlolohi[L][i+L] = vimlolohi[i];
            Yimhihilo[L][i+L] = vimhihilo[i];
            Yimlohilo[L][i+L] = vimlohilo[i];
            Yimhilolo[L][i+L] = vimhilolo[i];
            Yimlololo[L][i+L] = vimlololo[i];
         }
      }
      CPU_cmplx8_blocked_VB_to_W
         (nrows-k*szt,szt,Bhihihi,  Blohihi,  Bhilohi,  Blolohi,
                          Bhihilo,  Blohilo,  Bhilolo,  Blololo,
                        Yrehihihi,Yrelohihi,Yrehilohi,Yrelolohi,
                        Yrehihilo,Yrelohilo,Yrehilolo,Yrelololo,
                        Yimhihihi,Yimlohihi,Yimhilohi,Yimlolohi,
                        Yimhihilo,Yimlohilo,Yimhilolo,Yimlololo,
                        Wrehihihi,Wrelohihi,Wrehilohi,Wrelolohi,
                        Wrehihilo,Wrelohilo,Wrehilolo,Wrelololo,
                        Wimhihihi,Wimlohihi,Wimhilohi,Wimlolohi,
                        Wimhihilo,Wimlohilo,Wimhilolo,Wimlololo);
      if(k<nbt-1)
      {
         CPU_cmplx8_blocked_leftRupdate
            (nrows,ncols,szt,k,
             Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
             Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
             Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
             Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,
             Yrehihihi,Yrelohihi,Yrehilohi,Yrelolohi,
             Yrehihilo,Yrelohilo,Yrehilolo,Yrelololo,
             Yimhihihi,Yimlohihi,Yimhilohi,Yimlolohi,
             Yimhihilo,Yimlohilo,Yimhilolo,Yimlololo,
             Wrehihihi,Wrelohihi,Wrehilohi,Wrelolohi,
             Wrehihilo,Wrelohilo,Wrehilolo,Wrelololo,
             Wimhihihi,Wimlohihi,Wimhilohi,Wimlolohi,
             Wimhihilo,Wimlohilo,Wimhilolo,Wimlololo,verbose);
      }
      CPU_cmplx8_blocked_rightQupdate
         (nrows,szt,k,Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
                      Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
                      Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
                      Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
                      Yrehihihi,Yrelohihi,Yrehilohi,Yrelolohi,
                      Yrehihilo,Yrelohilo,Yrehilolo,Yrelololo,
                      Yimhihihi,Yimlohihi,Yimhilohi,Yimlolohi,
                      Yimhihilo,Yimlohilo,Yimhilolo,Yimlololo,
                      Wrehihihi,Wrelohihi,Wrehilohi,Wrelolohi,
                      Wrehihilo,Wrelohilo,Wrehilolo,Wrelololo,
                      Wimhihihi,Wimlohihi,Wimhilohi,Wimlolohi,
                      Wimhihilo,Wimlohilo,Wimhilolo,Wimlololo,verbose);
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(xrehihihi); free(xrelohihi); free(xrehilohi); free(xrelolohi);
   free(xrehihilo); free(xrelohilo); free(xrehilolo); free(xrelololo);
   free(ximhihihi); free(ximlohihi); free(ximhilohi); free(ximlolohi);
   free(ximhihilo); free(ximlohilo); free(ximhilolo); free(ximlololo);
   free(vrehihihi); free(vrelohihi); free(vrehilohi); free(vrelolohi); 
   free(vrehihilo); free(vrelohilo); free(vrehilolo); free(vrelololo); 
   free(vimhihihi); free(vimlohihi); free(vimhilohi); free(vimlolohi);
   free(vimhihilo); free(vimlohilo); free(vimhilolo); free(vimlololo);
   free(Bhihihi); free(Blohihi); free(Bhilohi); free(Blolohi);
   free(Bhihilo); free(Blohilo); free(Bhilolo); free(Blololo);

   for(int j=0; j<szt; j++)
   {
      free(Yrehihihi[j]); free(Yrelohihi[j]);
      free(Yrehilohi[j]); free(Yrelolohi[j]);
      free(Yrehihilo[j]); free(Yrelohilo[j]);
      free(Yrehilolo[j]); free(Yrelololo[j]);
      free(Yimhihihi[j]); free(Yimlohihi[j]);
      free(Yimhilohi[j]); free(Yimlolohi[j]);
      free(Yimhihilo[j]); free(Yimlohilo[j]);
      free(Yimhilolo[j]); free(Yimlololo[j]);
      free(Wrehihihi[j]); free(Wrelohihi[j]);
      free(Wrehilohi[j]); free(Wrelolohi[j]);
      free(Wrehihilo[j]); free(Wrelohilo[j]);
      free(Wrehilolo[j]); free(Wrelololo[j]);
      free(Wimhihihi[j]); free(Wimlohihi[j]);
      free(Wimhilohi[j]); free(Wimlolohi[j]);
      free(Wimhihilo[j]); free(Wimlohilo[j]);
      free(Wimhilolo[j]); free(Wimlololo[j]);
   }
   free(Yrehihihi); free(Yrelohihi); free(Yrehilohi); free(Yrelolohi);
   free(Yrehihilo); free(Yrelohilo); free(Yrehilolo); free(Yrelololo);
   free(Yimhihihi); free(Yimlohihi); free(Yimhilohi); free(Yimlolohi);
   free(Yimhihilo); free(Yimlohilo); free(Yimhilolo); free(Yimlololo);
   free(Wrehihihi); free(Wrelohihi); free(Wrehilohi); free(Wrelolohi);
   free(Wrehihilo); free(Wrelohilo); free(Wrehilolo); free(Wrelololo);
   free(Wimhihihi); free(Wimlohihi); free(Wimhilohi); free(Wimlolohi);
   free(Wimhihilo); free(Wimlohilo); free(Wimhilolo); free(Wimlololo);
}
