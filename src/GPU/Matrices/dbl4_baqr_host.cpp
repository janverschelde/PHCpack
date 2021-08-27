/* The file dbl4_baqr_host.cpp defines the functions specified in
 * the file dbl4_baqr_host.h. */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "quad_double_functions.h"
#include "dbl4_factorizations.h"
#include "dbl4_baqr_host.h"

using namespace std;

void CPU_dbl4_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double **Vhihi, double **Vlohi, double **Vhilo, double **Vlolo,
   double **Whihi, double **Wlohi, double **Whilo, double **Wlolo )
{
   double *vhihi = new double[nrows];
   double *vlohi = new double[nrows];
   double *vhilo = new double[nrows];
   double *vlolo = new double[nrows];
   double *phihi = new double[ncols];
   double *plohi = new double[ncols];
   double *philo = new double[ncols];
   double *plolo = new double[ncols];
   double zihihi,zilohi,zihilo,zilolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<nrows; i++)           //  W[0][i] = -B[0]*V[0][i];
   {
      qdf_mul(Bhihi[0],    Blohi[0],    Bhilo[0],    Blolo[0],
              Vhihi[0][i], Vlohi[0][i], Vhilo[0][i], Vlolo[0][i],
             &Whihi[0][i],&Wlohi[0][i],&Whilo[0][i],&Wlolo[0][i]);
      qdf_minus(&Whihi[0][i],&Wlohi[0][i],&Whilo[0][i],&Wlolo[0][i]);
   }

   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++)
      {
         vhihi[i] = Vhihi[j][i];
         vlohi[i] = Vlohi[j][i];
         vhilo[i] = Vhilo[j][i];
         vlolo[i] = Vlolo[j][i];
      }
      for(int k=0; k<j; k++)
      {
         phihi[k] = 0.0;                // compute k-th component of Y^T*v
         plohi[k] = 0.0;                // over all rows of k-th column of Y
         philo[k] = 0.0;
         plolo[k] = 0.0;

         for(int i=0; i<nrows; i++)     // p[k] = p[k] + V[k][i]*v[i]
         {
            qdf_mul(Vhihi[k][i],Vlohi[k][i],Vhilo[k][i],Vlolo[k][i],
                    vhihi[i],   vlohi[i],   vhilo[i],   vlolo[i],
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
            qdf_inc(&phihi[k],&plohi[k],&philo[k],&plolo[k],
                   acchihi,  acclohi,  acchilo,  acclolo);
         }
      }
      for(int i=0; i<nrows; i++)
      {
         zihihi = 0.0;                  // compute i-th component of W*p
         zilohi = 0.0;                  // over all rows of k-th column of W
         zihilo = 0.0;
         zilolo = 0.0;

         for(int k=0; k<j; k++)         // zi = zi + W[k][i]*p[k]
         {
            qdf_mul(Whihi[k][i],Wlohi[k][i],Whilo[k][i],Wlolo[k][i],
                    phihi[k],   plohi[k],   philo[k],   plolo[k],
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
            qdf_inc(&zihihi,&zilohi,&zihilo,&zilolo,
                    acchihi,acclohi,acchilo,acclolo);
         }
         qdf_inc(&zihihi, &zilohi, &zihilo, &zilolo,
                   vhihi[i],vlohi[i],vhilo[i],vlolo[i]);
         // W[j][i] = -B[j]*zi;
         qdf_mul(Bhihi[j],    Blohi[j],    Bhilo[j],    Blolo[j],
                zihihi,      zilohi    ,  zihilo,      zilolo,
                &Whihi[j][i],&Wlohi[j][i],&Whilo[j][i],&Wlolo[j][i]);
         qdf_minus(&Whihi[j][i],&Wlohi[j][i],&Whilo[j][i],&Wlolo[j][i]);
      }
   }
   free(vhihi); free(vlohi); free(vhilo); free(vlolo);
   free(phihi); free(plohi); free(philo); free(plolo);
}

void CPU_cmplx4_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double **Vrehihi, double **Vrelohi, double **Vrehilo, double **Vrelolo,
   double **Vimhihi, double **Vimlohi, double **Vimhilo, double **Vimlolo,
   double **Wrehihi, double **Wrelohi, double **Wrehilo, double **Wrelolo,
   double **Wimhihi, double **Wimlohi, double **Wimhilo, double **Wimlolo )
{
   double *vrehihi = new double[nrows];
   double *vrelohi = new double[nrows];
   double *vrehilo = new double[nrows];
   double *vrelolo = new double[nrows];
   double *vimhihi = new double[nrows];
   double *vimlohi = new double[nrows];
   double *vimhilo = new double[nrows];
   double *vimlolo = new double[nrows];
   double *prehihi = new double[ncols];
   double *prelohi = new double[ncols];
   double *prehilo = new double[ncols];
   double *prelolo = new double[ncols];
   double *pimhihi = new double[ncols];
   double *pimlohi = new double[ncols];
   double *pimhilo = new double[ncols];
   double *pimlolo = new double[ncols];
   double acchihi,acclohi,acchilo,acclolo;
   double zi_rehihi,zi_relohi,zi_rehilo,zi_relolo;
   double zi_imhihi,zi_imlohi,zi_imhilo,zi_imlolo;

   for(int i=0; i<nrows; i++) // W[0][i] = -B[0]*V[0][i]
   {
      qdf_mul(Bhihi[0],      Blohi[0],      Bhilo[0],      Blolo[0],
            Vrehihi[0][i], Vrelohi[0][i], Vrehilo[0][i], Vrelolo[0][i],
           &Wrehihi[0][i],&Wrelohi[0][i],&Wrehilo[0][i],&Wrelolo[0][i]);
      qdf_minus(&Wrehihi[0][i],&Wrelohi[0][i],&Wrehilo[0][i],&Wrelolo[0][i]);
      qdf_mul(Bhihi[0],      Blohi[0],      Bhilo[0],      Blolo[0],
            Vimhihi[0][i], Vimlohi[0][i], Vimhilo[0][i], Vimlolo[0][i],
           &Wimhihi[0][i],&Wimlohi[0][i],&Wimhilo[0][i],&Wimlolo[0][i]);
      qdf_minus(&Wimhihi[0][i],&Wimlohi[0][i],&Wimhilo[0][i],&Wimlolo[0][i]);
   }
   for(int j=1; j<ncols; j++)           // compute column j of W
   {
      for(int i=0; i<nrows; i++)
      {
         vrehihi[i] = Vrehihi[j][i]; vrelohi[i] = Vrelohi[j][i];
         vrehilo[i] = Vrehilo[j][i]; vrelolo[i] = Vrelolo[j][i];
         vimhihi[i] = Vimhihi[j][i]; vimlohi[i] = Vimlohi[j][i];
         vimhilo[i] = Vimhilo[j][i]; vimlolo[i] = Vimlolo[j][i];
      }
      for(int k=0; k<j; k++)
      {
         prehihi[k] = 0.0;            // compute k-th component of Y^H*v
         prelohi[k] = 0.0;
         prehilo[k] = 0.0;
         prelolo[k] = 0.0;
         pimhihi[k] = 0.0;
         pimlohi[k] = 0.0;
         pimhilo[k] = 0.0;
         pimlolo[k] = 0.0;

         for(int i=0; i<nrows; i++)    // p[k] = p[k] + V[k][i]*v[i];
         {
            // accre =   Vre[k][i]*vre[i] + Vim[k][i]*vim[i];
            // pre[k] = pre[k] + accre;
            qdf_mul(Vrehihi[k][i],Vrelohi[k][i],Vrehilo[k][i],Vrelolo[k][i],
                    vrehihi[i],   vrelohi[i],   vrehilo[i],   vrelolo[i],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&prehihi[k],&prelohi[k],&prehilo[k],&prelolo[k],
                     acchihi,    acclohi,    acchilo,    acclolo);
            qdf_mul(Vimhihi[k][i],Vimlohi[k][i],Vimhilo[k][i],Vimlolo[k][i],
                    vimhihi[i],   vimlohi[i],   vimhilo[i],   vimlolo[i],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&prehihi[k],&prelohi[k],&prehilo[k],&prelolo[k],
                     acchihi,    acclohi,    acchilo,    acclolo);
            // accim = - Vim[k][i]*vre[i] + Vre[k][i]*vim[i];
            // pim[k] = pim[k] + accim;
            qdf_mul(Vimhihi[k][i],Vimlohi[k][i],Vimhilo[k][i],Vimlolo[k][i],
                    vrehihi[i],   vrelohi[i],   vrehilo[i],   vrelolo[i],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_dec(&pimhihi[k],&pimlohi[k],&pimhilo[k],&pimlolo[k],
                     acchihi,    acclohi,    acchilo,    acclolo);
            qdf_mul(Vrehihi[k][i],Vrelohi[k][i],Vrehilo[k][i],Vrelolo[k][i],
                    vimhihi[i],   vimlohi[i],   vimhilo[i],   vimlolo[i],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&pimhihi[k],&pimlohi[k],&pimhilo[k],&pimlolo[k],
                     acchihi,    acclohi,    acchilo,    acclolo);
         }
      }
      for(int i=0; i<nrows; i++)
      {
         zi_rehihi = 0.0;               // compute i-th component of W*p
         zi_relohi = 0.0;
         zi_rehilo = 0.0;
         zi_relolo = 0.0;
         zi_imhihi = 0.0;
         zi_imlohi = 0.0;
         zi_imhilo = 0.0;
         zi_imlolo = 0.0;
 
         for(int k=0; k<j; k++)         // zi = zi + W[k][i]*p[k];
         {
            // accre = Wre[k][i]*pre[k] - Wim[k][i]*pim[k];
            // zi_re = zi_re + accre;
            qdf_mul(Wrehihi[k][i],Wrelohi[k][i],Wrehilo[k][i],Wrelolo[k][i],
                    prehihi[k],   prelohi[k],   prehilo[k],   prelolo[k],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&zi_rehihi,&zi_relohi,&zi_rehilo,&zi_relolo,
                       acchihi,   acclohi,   acchilo,   acclolo);
            qdf_mul(Wimhihi[k][i],Wimlohi[k][i],Wimhilo[k][i],Wimlolo[k][i],
                    pimhihi[k],   pimlohi[k],   pimhilo[k],   pimlolo[k],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_dec(&zi_rehihi,&zi_relohi,&zi_rehilo,&zi_relolo,
                       acchihi,   acclohi,   acchilo,   acclolo);
            // accim = Wim[k][i]*pre[k] + Wre[k][i]*pim[k];
            // zi_im = zi_im + accim;
            qdf_mul(Wimhihi[k][i],Wimlohi[k][i],Wimhilo[k][i],Wimlolo[k][i],
                    prehihi[k],   prelohi[k],   prehilo[k],   prelolo[k],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&zi_imhihi,&zi_imlohi,&zi_imhilo,&zi_imlolo,
                       acchihi,   acclohi,   acchilo,   acclolo);
            qdf_mul(Wrehihi[k][i],Wrelohi[k][i],Wrehilo[k][i],Wrelolo[k][i],
                    pimhihi[k],   pimlohi[k],   pimhilo[k],   pimlolo[k],
                   &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&zi_imhihi,&zi_imlohi,&zi_imhilo,&zi_imlolo,
                       acchihi,   acclohi,   acchilo,   acclolo);
         }
         // zi_re = zi_re + vre[i];
         qdf_inc(&zi_rehihi,&zi_relohi,&zi_rehilo,&zi_relolo,
                    vrehihi[i],vrelohi[i],vrehilo[i],vrelolo[i]);
         // zi_im = zi_im + vim[i];
         qdf_inc(&zi_imhihi,&zi_imlohi,&zi_imhilo,&zi_imlolo,
                    vimhihi[i],vimlohi[i],vimhilo[i],vimlolo[i]);
         // Wre[j][i] = -B[j]*zi_re;
         qdf_mul(Bhihi[j],      Blohi[j],      Bhilo[j],      Blolo[j],
             zi_rehihi,     zi_relohi,     zi_rehilo,     zi_relolo,
              &Wrehihi[j][i],&Wrelohi[j][i],&Wrehilo[j][i],&Wrelolo[j][i]);
         qdf_minus(&Wrehihi[j][i],&Wrelohi[j][i],
                   &Wrehilo[j][i],&Wrelolo[j][i]);
         // Wim[j][i] = -B[j]*zi_im;
         qdf_mul(Bhihi[j],      Blohi[j],      Bhilo[j],      Blolo[j],
             zi_imhihi,     zi_imlohi,     zi_imhilo,     zi_imlolo,
              &Wimhihi[j][i],&Wimlohi[j][i],&Wimhilo[j][i],&Wimlolo[j][i]);
         qdf_minus(&Wimhihi[j][i],&Wimlohi[j][i],
                   &Wimhilo[j][i],&Wimlolo[j][i]);
      }
   }
   free(vrehihi); free(vrelohi); free(vrehilo); free(vrelolo);
   free(vimhihi); free(vimlohi); free(vimhilo); free(vimlolo);
   free(prehihi); free(prelohi); free(prehilo); free(prelolo);
   free(pimhihi); free(pimlohi); free(pimhilo); free(pimlolo);
}

void CPU_dbl4_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Chihi, double **Clohi, double **Chilo, double **Clolo,
   double **Yhihi, double **Ylohi, double **Yhilo, double **Ylolo,
   double **Whihi, double **Wlohi, double **Whilo, double **Wlolo,
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
   double **YWThihi = new double*[rowdim];
   double **YWTlohi = new double*[rowdim];
   double **YWThilo = new double*[rowdim];
   double **YWTlolo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      YWThihi[i] = new double[rowdim];
      YWTlohi[i] = new double[rowdim];
      YWThilo[i] = new double[rowdim];
      YWTlolo[i] = new double[rowdim];
   }
   double **prdhihi = new double*[rowdim];
   double **prdlohi = new double*[rowdim];
   double **prdhilo = new double*[rowdim];
   double **prdlolo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      prdhihi[i] = new double[coldim];
      prdlohi[i] = new double[coldim];
      prdhilo[i] = new double[coldim];
      prdlolo[i] = new double[coldim];
   }
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<rowdim; i++)        // compute Y*W^T
      for(int j=0; j<rowdim; j++)
      {
         YWThihi[i][j] = 0.0;         // row i of Y with column j of W^T
         YWTlohi[i][j] = 0.0;
         YWThilo[i][j] = 0.0;
         YWTlolo[i][j] = 0.0;

         for(int k=0; k<szt; k++)     // YWT[i][j] += Y[k][i]*W[k][j]
         {
            qdf_mul(Yhihi[k][i],Ylohi[k][i],Yhilo[k][i],Ylolo[k][i],
                    Whihi[k][j],Wlohi[k][j],Whilo[k][j],Wlolo[k][j],
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
            qdf_inc(&YWThihi[i][j],&YWTlohi[i][j],
                    &YWThilo[i][j],&YWTlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWThihi[i][j] << "  " << YWTlohi[i][j] << endl
                 << "            "
                 << YWThilo[i][j] << "  " << YWTlolo[i][j] << endl;
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         prdhihi[i][j] = 0.0;
         prdlohi[i][j] = 0.0;
         prdhilo[i][j] = 0.0;
         prdlolo[i][j] = 0.0; 
                                      // prd[i][j]
         for(int k=0; k<rowdim; k++)  // += YWT[i][k]*C[rowoff+k][coloff+j]
         {
            qdf_mul(YWThihi[i][k],YWTlohi[i][k],YWThilo[i][k],YWTlolo[i][k],
                      Chihi[rowoff+k][coloff+j],Clohi[rowoff+k][coloff+j],
                      Chilo[rowoff+k][coloff+j],Clolo[rowoff+k][coloff+j],
                   &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdhihi[i][j],&prdlohi[i][j],
                    &prdhilo[i][j],&prdlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   for(int i=0; i<rowdim; i++)      // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)   // C[rowoff+i][coloff+j] += prd[i][j]
      {
         qdf_inc(&Chihi[rowoff+i][coloff+j],&Clohi[rowoff+i][coloff+j],
                 &Chilo[rowoff+i][coloff+j],&Clolo[rowoff+i][coloff+j],
                 prdhihi[i][j],prdlohi[i][j],prdhilo[i][j],prdlolo[i][j]);
      }

   for(int i=0; i<rowdim; i++)
   {
      free(YWThihi[i]); free(YWTlohi[i]); free(YWThilo[i]); free(YWTlolo[i]);
      free(prdhihi[i]); free(prdlohi[i]); free(prdhilo[i]); free(prdlolo[i]);
   }
   free(YWThihi); free(YWTlohi); free(YWThilo); free(YWTlolo);
   free(prdhihi); free(prdlohi); free(prdhilo); free(prdlolo);
}

void CPU_cmplx4_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Crehihi, double **Crelohi, double **Crehilo, double **Crelolo,
   double **Cimhihi, double **Cimlohi, double **Cimhilo, double **Cimlolo,
   double **Yrehihi, double **Yrelohi, double **Yrehilo, double **Yrelolo,
   double **Yimhihi, double **Yimlohi, double **Yimhilo, double **Yimlolo,
   double **Wrehihi, double **Wrelohi, double **Wrehilo, double **Wrelolo,
   double **Wimhihi, double **Wimlohi, double **Wimhilo, double **Wimlolo,
   bool verbose )
{
   const int rowoff = idx*szt;             // row offset for C
   const int rowdim = nrows - rowoff;      // number of rows in Y and W
   const int coloff = (idx+1)*szt;         // column offset for C
   const int coldim = ncols - coloff;      // number of columns in C
   double acchihi,acclohi,acchilo,acclolo;

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
   double **YWTrehihi = new double*[rowdim];
   double **YWTrelohi = new double*[rowdim];
   double **YWTrehilo = new double*[rowdim];
   double **YWTrelolo = new double*[rowdim];
   double **YWTimhihi = new double*[rowdim];
   double **YWTimlohi = new double*[rowdim];
   double **YWTimhilo = new double*[rowdim];
   double **YWTimlolo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      YWTrehihi[i] = new double[rowdim];
      YWTrelohi[i] = new double[rowdim];
      YWTrehilo[i] = new double[rowdim];
      YWTrelolo[i] = new double[rowdim];
      YWTimhihi[i] = new double[rowdim];
      YWTimlohi[i] = new double[rowdim];
      YWTimhilo[i] = new double[rowdim];
      YWTimlolo[i] = new double[rowdim];
   }
   double **prdrehihi = new double*[rowdim];
   double **prdrelohi = new double*[rowdim];
   double **prdrehilo = new double*[rowdim];
   double **prdrelolo = new double*[rowdim];
   double **prdimhihi = new double*[rowdim];
   double **prdimlohi = new double*[rowdim];
   double **prdimhilo = new double*[rowdim];
   double **prdimlolo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      prdrehihi[i] = new double[coldim];
      prdrelohi[i] = new double[coldim];
      prdrehilo[i] = new double[coldim];
      prdrelolo[i] = new double[coldim];
      prdimhihi[i] = new double[coldim];
      prdimlohi[i] = new double[coldim];
      prdimhilo[i] = new double[coldim];
      prdimlolo[i] = new double[coldim];
   }
   for(int i=0; i<rowdim; i++)      // compute Y*W^H
      for(int j=0; j<rowdim; j++)
      {
         YWTrehihi[i][j] = 0.0;       // row i of Y with column j of W^H
         YWTrelohi[i][j] = 0.0;
         YWTrehilo[i][j] = 0.0;
         YWTrelolo[i][j] = 0.0;
         YWTimhihi[i][j] = 0.0;
         YWTimlohi[i][j] = 0.0;
         YWTimhilo[i][j] = 0.0;
         YWTimlolo[i][j] = 0.0;

         for(int k=0; k<szt; k++)   // YWT[i][j] = YWT[i][j] + Y[k][i]*W[k][j]
         {  
            // accre = Yre[k][i]*Wre[k][j] + Yim[k][i]*Wim[k][j];
            // YWTre[i][j] = YWTre[i][j] + accre;
            qdf_mul(Yrehihi[k][i],Yrelohi[k][i],Yrehilo[k][i],Yrelolo[k][i],
                    Wrehihi[k][j],Wrelohi[k][j],Wrehilo[k][j],Wrelolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&YWTrehihi[i][j],&YWTrelohi[i][j],
                    &YWTrehilo[i][j],&YWTrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(Yimhihi[k][i],Yimlohi[k][i],Yimhilo[k][i],Yimlolo[k][i],
                    Wimhihi[k][j],Wimlohi[k][j],Wimhilo[k][j],Wimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&YWTrehihi[i][j],&YWTrelohi[i][j],
                    &YWTrehilo[i][j],&YWTrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            // accim = Yim[k][i]*Wre[k][j] - Yre[k][i]*Wim[k][j];
            // YWTim[i][j] = YWTim[i][j] + accim;
            qdf_mul(Yimhihi[k][i],Yimlohi[k][i],Yimhilo[k][i],Yimlolo[k][i],
                    Wrehihi[k][j],Wrelohi[k][j],Wrehilo[k][j],Wrelolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&YWTimhihi[i][j],&YWTimlohi[i][j],
                    &YWTimhilo[i][j],&YWTimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(Yrehihi[k][i],Yrelohi[k][i],Yrehilo[k][i],Yrelolo[k][i],
                    Wimhihi[k][j],Wimlohi[k][j],Wimhilo[k][j],Wimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_dec(&YWTimhihi[i][j],&YWTimlohi[i][j],
                    &YWTimhilo[i][j],&YWTimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "YWT[" << i << "][" << j << "]re : "
                 << YWTrehihi[i][j] << "  " << YWTrelohi[i][j] << endl
                 << "              "
                 << YWTrehilo[i][j] << "  " << YWTrelolo[i][j] << endl;
            cout << "YWT[" << i << "][" << j << "]im : "
                 << YWTimhihi[i][j] << "  " << YWTimlohi[i][j] << endl
                 << "              "
                 << YWTimhilo[i][j] << "  " << YWTimlolo[i][j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)        // prd = (Y*W^H)*C
      for(int j=0; j<coldim; j++)
      {
         prdrehihi[i][j] = 0.0;
         prdrelohi[i][j] = 0.0;
         prdrehilo[i][j] = 0.0;
         prdrelolo[i][j] = 0.0;
         prdimhihi[i][j] = 0.0;
         prdimlohi[i][j] = 0.0;
         prdimhilo[i][j] = 0.0;
         prdimlolo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++)
         {  // prd[i][j] = prd[i][j] + YWT[i][k]*C[rowoff+k][coloff+j];
            // accre = YWTre[i][k]*Cre[rowoff+k][coloff+j]
            //       - YWTim[i][k]*Cim[rowoff+k][coloff+j];
            // prdre[i][j] = prdre[i][j] + accre;
            qdf_mul(YWTrehihi[i][k],YWTrelohi[i][k],
                    YWTrehilo[i][k],YWTrelolo[i][k],
                    Crehihi[rowoff+k][coloff+j],Crelohi[rowoff+k][coloff+j],
                    Crehilo[rowoff+k][coloff+j],Crelolo[rowoff+k][coloff+j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdrehihi[i][j],&prdrelohi[i][j],
                    &prdrehilo[i][j],&prdrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(YWTimhihi[i][k],YWTimlohi[i][k],
                    YWTimhilo[i][k],YWTimlolo[i][k],
                    Cimhihi[rowoff+k][coloff+j],Cimlohi[rowoff+k][coloff+j],
                    Cimhilo[rowoff+k][coloff+j],Cimlolo[rowoff+k][coloff+j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_dec(&prdrehihi[i][j],&prdrelohi[i][j],
                    &prdrehilo[i][j],&prdrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            // accim = YWTim[i][k]*Cre[rowoff+k][coloff+j]
            //       + YWTre[i][k]*Cim[rowoff+k][coloff+j];
            // prdim[i][j] = prdim[i][j] + accim;
            qdf_mul(YWTimhihi[i][k],YWTimlohi[i][k],
                    YWTimhilo[i][k],YWTimlolo[i][k],
                    Crehihi[rowoff+k][coloff+j],Crelohi[rowoff+k][coloff+j],
                    Crehilo[rowoff+k][coloff+j],Crelolo[rowoff+k][coloff+j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdimhihi[i][j],&prdimlohi[i][j],
                    &prdimhilo[i][j],&prdimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(YWTrehihi[i][k],YWTrelohi[i][k],
                    YWTrehilo[i][k],YWTrelolo[i][k],
                    Cimhihi[rowoff+k][coloff+j],Cimlohi[rowoff+k][coloff+j],
                    Cimhilo[rowoff+k][coloff+j],Cimlolo[rowoff+k][coloff+j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdimhihi[i][j],&prdimlohi[i][j],
                    &prdimhilo[i][j],&prdimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   if(verbose)
   {
      cout << endl;
      for(int i=0; i<rowdim; i++)
         for(int j=0; j<coldim; j++)
         {
            cout << "prd[" << i << "][" << j << "]re : "
                 << prdrehihi[i][j] << "  " << prdrelohi[i][j] << endl
                 << "              "
                 << prdrehilo[i][j] << "  " << prdrelolo[i][j] << endl;
            cout << "prd[" << i << "][" << j << "]im : "
                 << prdimhihi[i][j] << "  " << prdimlohi[i][j] << endl
                 << "              "
                 << prdimhilo[i][j] << "  " << prdimlolo[i][j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)        // C = C + (Y*W^T)*C
      for(int j=0; j<coldim; j++)
      {
         qdf_inc(&Crehihi[rowoff+i][coloff+j],&Crelohi[rowoff+i][coloff+j],
                 &Crehilo[rowoff+i][coloff+j],&Crelolo[rowoff+i][coloff+j],
                 prdrehihi[i][j],prdrelohi[i][j],
                 prdrehilo[i][j],prdrelolo[i][j]);
         qdf_inc(&Cimhihi[rowoff+i][coloff+j],&Cimlohi[rowoff+i][coloff+j],
                 &Cimhilo[rowoff+i][coloff+j],&Cimlolo[rowoff+i][coloff+j],
                 prdimhihi[i][j],prdimlohi[i][j],
                 prdimhilo[i][j],prdimlolo[i][j]);
      }

   for(int i=0; i<rowdim; i++)
   {
      free(YWTrehihi[i]); free(YWTrelohi[i]);
      free(YWTrehilo[i]); free(YWTrelolo[i]);
      free(YWTimhihi[i]); free(YWTimlohi[i]);
      free(YWTimhilo[i]); free(YWTimlolo[i]);
      free(prdrehihi[i]); free(prdrelohi[i]);
      free(prdrehilo[i]); free(prdrelolo[i]);
      free(prdimhihi[i]); free(prdimlohi[i]);
      free(prdimhilo[i]); free(prdimlolo[i]);
   }
   free(prdrehihi); free(prdrelohi); free(prdrehilo); free(prdrelolo);
   free(prdimhihi); free(prdimlohi); free(prdimhilo); free(prdimlolo);
   free(YWTrehihi); free(YWTrelohi); free(YWTrehilo); free(YWTrelolo);
   free(YWTimhihi); free(YWTimlohi); free(YWTimhilo); free(YWTimlolo);
}

void CPU_dbl4_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Yhihi, double **Ylohi, double **Yhilo, double **Ylolo,
   double **Whihi, double **Wlohi, double **Whilo, double **Wlolo,
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
   double **WYThihi = new double*[rowdim];
   double **WYTlohi = new double*[rowdim];
   double **WYThilo = new double*[rowdim];
   double **WYTlolo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      WYThihi[i] = new double[rowdim];
      WYTlohi[i] = new double[rowdim];
      WYThilo[i] = new double[rowdim];
      WYTlolo[i] = new double[rowdim];
   }
   double **prdhihi = new double*[dim];
   double **prdlohi = new double*[dim];
   double **prdhilo = new double*[dim];
   double **prdlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdhihi[i] = new double[rowdim];
      prdlohi[i] = new double[rowdim];
      prdhilo[i] = new double[rowdim];
      prdlolo[i] = new double[rowdim];
   }
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<rowdim; i++)     // compute W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         WYThihi[i][j] = 0.0;        // row i of W with column j of Y^T
         WYTlohi[i][j] = 0.0;
         WYThilo[i][j] = 0.0; 
         WYTlolo[i][j] = 0.0;        // take k-th column of W

         for(int k=0; k<szt; k++)  // WYT[i][j] = WYT[i][j] + W[k][i]*Y[k][j]
         {
            qdf_mul(Whihi[k][i],Wlohi[k][i],Whilo[k][i],Wlolo[k][i],
                    Yhihi[k][j],Ylohi[k][j],Yhilo[k][j],Ylolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&WYThihi[i][j],&WYTlohi[i][j],
                    &WYThilo[i][j],&WYTlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   if(verbose)
   {
      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "Y[" << i << "][" << j << "] : "
                 << Yhihi[j][i] << "  " << Ylohi[j][i] << endl
                 << "          "
                 << Yhilo[j][i] << "  " << Ylolo[j][i] << endl;

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
            cout << "W[" << i << "][" << j << "] : "
                 << Whihi[j][i] << "  " << Wlohi[j][i] << endl
                 << "          "
                 << Whilo[j][i] << "  " << Wlolo[j][i] << endl;

      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihi[i][j] << "  " << WYTlohi[i][j] << endl
                 << "            "
                 << WYThilo[i][j] << "  " << WYTlolo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)          // prd = Q*W*Y^T
      for(int j=0; j<rowdim; j++)
      {
         prdhihi[i][j] = 0.0;
         prdlohi[i][j] = 0.0;
         prdhilo[i][j] = 0.0;
         prdlolo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++) // prd[i][j] += Q[i][coloff+k]*WYT[k][j]
         {
            qdf_mul(Qhihi[i][coloff+k],Qlohi[i][coloff+k],
                    Qhilo[i][coloff+k],Qlolo[i][coloff+k],
                    WYThihi[k][j],WYTlohi[k][j],WYThilo[k][j],WYTlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdhihi[i][j],&prdlohi[i][j],
                    &prdhilo[i][j],&prdlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   if(verbose)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "QWYT[" << i << "][" << j << "] : "
                 << prdhihi[i][j] << "  " << prdlohi[i][j] << endl
                 << "             "
                 << prdhilo[i][j] << "  " << prdlolo[i][j] << endl;

      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Qhihi[i][coloff+j] << "  "
                 << Qlohi[i][coloff+j] << endl
                 << "          "
                 << Qhilo[i][coloff+j] << "  "
                 << Qlolo[i][coloff+j] << endl;
   }
   for(int i=0; i<dim; i++)        // Q = Q + Q*W*Y^T
      for(int j=0; j<rowdim; j++)  // Q[i][coloff+j] += prd[i][j];
      {
         qdf_inc(&Qhihi[i][coloff+j],&Qlohi[i][coloff+j],
                 &Qhilo[i][coloff+j],&Qlolo[i][coloff+j],
                 prdhihi[i][j],prdlohi[i][j],prdhilo[i][j],prdlolo[i][j]);
      }

   if(verbose)
   {
      cout << "Q after the update with QWYT :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
            cout << "Q[" << i << "][" << coloff+j << "] : "
                 << Qhihi[i][coloff+j] << "  "
                 << Qlohi[i][coloff+j] << endl
                 << "          "
                 << Qhilo[i][coloff+j] << "  "
                 << Qlolo[i][coloff+j] << endl;
   }
   for(int i=0; i<rowdim; i++)
   {
      free(WYThihi[i]); free(WYTlohi[i]);
      free(WYThilo[i]); free(WYTlolo[i]);
   }
   for(int i=0; i<dim; i++)
   {
      free(prdhihi[i]); free(prdlohi[i]);
      free(prdhilo[i]); free(prdlolo[i]);
   }
   free(WYThihi); free(WYTlohi); free(WYThilo); free(WYTlolo);
   free(prdhihi); free(prdlohi); free(prdhilo); free(prdlolo);
}

void CPU_cmplx4_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Yrehihi, double **Yrelohi, double **Yrehilo, double **Yrelolo,
   double **Yimhihi, double **Yimlohi, double **Yimhilo, double **Yimlolo,
   double **Wrehihi, double **Wrelohi, double **Wrehilo, double **Wrelolo,
   double **Wimhihi, double **Wimlohi, double **Wimhilo, double **Wimlolo,
   bool verbose )
{
   const int coloff = idx*szt;        // column offset for Q
   const int rowdim = dim - coloff;
   // the number of rows in Y and W is the number of columns to update
   double acchihi,acclohi,acchilo,acclolo;

   if(verbose)
   {
      cout << "updating Q ..." << endl;
      cout << "-> dim : " << dim << "  szt : " << szt << "  idx : " << idx
           << "  rowdim : " << rowdim << "  coloff : " << coloff << endl;
   }
   double **WYTrehihi = new double*[rowdim];
   double **WYTrelohi = new double*[rowdim];
   double **WYTrehilo = new double*[rowdim];
   double **WYTrelolo = new double*[rowdim];
   double **WYTimhihi = new double*[rowdim];
   double **WYTimlohi = new double*[rowdim];
   double **WYTimhilo = new double*[rowdim];
   double **WYTimlolo = new double*[rowdim];

   for(int i=0; i<rowdim; i++)
   {
      WYTrehihi[i] = new double[rowdim];
      WYTrelohi[i] = new double[rowdim];
      WYTrehilo[i] = new double[rowdim];
      WYTrelolo[i] = new double[rowdim];
      WYTimhihi[i] = new double[rowdim];
      WYTimlohi[i] = new double[rowdim];
      WYTimhilo[i] = new double[rowdim];
      WYTimlolo[i] = new double[rowdim];
   }
   double **prdrehihi = new double*[dim];
   double **prdrelohi = new double*[dim];
   double **prdrehilo = new double*[dim];
   double **prdrelolo = new double*[dim];
   double **prdimhihi = new double*[dim];
   double **prdimlohi = new double*[dim];
   double **prdimhilo = new double*[dim];
   double **prdimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrehihi[i] = new double[rowdim];
      prdrelohi[i] = new double[rowdim];
      prdrehilo[i] = new double[rowdim];
      prdrelolo[i] = new double[rowdim];
      prdimhihi[i] = new double[rowdim];
      prdimlohi[i] = new double[rowdim];
      prdimhilo[i] = new double[rowdim];
      prdimlolo[i] = new double[rowdim];
   }
   for(int i=0; i<rowdim; i++)        // compute W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         WYTrehihi[i][j] = 0.0;       // row i of W with column j of Y^H
         WYTrelohi[i][j] = 0.0;
         WYTrehilo[i][j] = 0.0;
         WYTrelolo[i][j] = 0.0;
         WYTimhihi[i][j] = 0.0;
         WYTimlohi[i][j] = 0.0;
         WYTimhilo[i][j] = 0.0;
         WYTimlolo[i][j] = 0.0;

         for(int k=0; k<szt; k++)     // WYT[i][j] += W[k][i]*Y[k][j]
         {
            // accre = Wre[k][i]*Yre[k][j] + Wim[k][i]*Yim[k][j];
            // WYTre[i][j] = WYTre[i][j] + accre;
            qdf_mul(Wrehihi[k][i],Wrelohi[k][i],Wrehilo[k][i],Wrelolo[k][i],
                    Yrehihi[k][j],Yrelohi[k][j],Yrehilo[k][j],Yrelolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&WYTrehihi[i][j],&WYTrelohi[i][j],
                    &WYTrehilo[i][j],&WYTrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(Wimhihi[k][i],Wimlohi[k][i],Wimhilo[k][i],Wimlolo[k][i],
                    Yimhihi[k][j],Yimlohi[k][j],Yimhilo[k][j],Yimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&WYTrehihi[i][j],&WYTrelohi[i][j],
                    &WYTrehilo[i][j],&WYTrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            // accim = Wim[k][i]*Yre[k][j] - Wre[k][i]*Yim[k][j];
            // WYTim[i][j] = WYTim[i][j] + accim;
            qdf_mul(Wimhihi[k][i],Wimlohi[k][i],Wimhilo[k][i],Wimlolo[k][i],
                    Yrehihi[k][j],Yrelohi[k][j],Yrehilo[k][j],Yrelolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&WYTimhihi[i][j],&WYTimlohi[i][j],
                    &WYTimhilo[i][j],&WYTimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(Wrehihi[k][i],Wrelohi[k][i],Wrehilo[k][i],Wrelolo[k][i],
                    Yimhihi[k][j],Yimlohi[k][j],Yimhilo[k][j],Yimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_dec(&WYTimhihi[i][j],&WYTimlohi[i][j],
                    &WYTimhilo[i][j],&WYTimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   if(verbose)
   {
      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
         {
            cout << "Y[" << i << "][" << j << "]re : "
                 << Yrehihi[j][i] << "  " << Yrelohi[j][i] << endl
                 << "            "
                 << Yrehilo[j][i] << "  " << Yrelolo[j][i] << endl;
            cout << "Y[" << i << "][" << j << "]im : "
                 << Yimhihi[j][i] << "  " << Yimlohi[j][i] << endl
                 << "            "
                 << Yimhilo[j][i] << "  " << Yimlolo[j][i] << endl;
         }

      for(int j=0; j<szt; j++)
         for(int i=0; i<rowdim; i++)
         {
            cout << "W[" << i << "][" << j << "]re : "
                 << Wrehihi[j][i] << "  " << Wrelohi[j][i] << endl
                 << "            "
                 << Wrehilo[j][i] << "  " << Wrelolo[j][i] << endl;
            cout << "W[" << i << "][" << j << "]im : "
                 << Wimhihi[j][i] << "  " << Wimlohi[j][i] << endl
                 << "            "
                 << Wimhilo[j][i] << "  " << Wimlolo[j][i] << endl;
         }

      for(int i=0; i<rowdim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "WYT[" << i << "][" << j << "]re : "
                 << WYTrehihi[i][j] << "  " << WYTrelohi[i][j] << endl
                 << "            "
                 << WYTrehilo[i][j] << "  " << WYTrelolo[i][j] << endl;
            cout << "WYT[" << i << "][" << j << "]im : "
                 << WYTimhihi[i][j] << "  " << WYTimlohi[i][j] << endl
                 << "            "
                 << WYTimhilo[i][j] << "  " << WYTimlolo[i][j] << endl;
         }
   }
   for(int i=0; i<dim; i++)           // prd = Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         prdrehihi[i][j] = 0.0;
         prdrelohi[i][j] = 0.0;
         prdrehilo[i][j] = 0.0;
         prdrelolo[i][j] = 0.0;
         prdimhihi[i][j] = 0.0;
         prdimlohi[i][j] = 0.0;
         prdimhilo[i][j] = 0.0;
         prdimlolo[i][j] = 0.0;

         for(int k=0; k<rowdim; k++) // prd[i][j] += Q[i][coloff+k]*WYT[k][j]
         {
            // accre = Qre[i][coloff+k]*WYTre[k][j]
            //       - Qim[i][coloff+k]*WYTim[k][j];
            // prdre[i][j] = prdre[i][j] + accre;
            qdf_mul(Qrehihi[i][coloff+k],Qrelohi[i][coloff+k],
                    Qrehilo[i][coloff+k],Qrelolo[i][coloff+k],
                    WYTrehihi[k][j],WYTrelohi[k][j],
                    WYTrehilo[k][j],WYTrelolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdrehihi[i][j],&prdrelohi[i][j],
                    &prdrehilo[i][j],&prdrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(Qimhihi[i][coloff+k],Qimlohi[i][coloff+k],
                    Qimhilo[i][coloff+k],Qimlolo[i][coloff+k],
                    WYTimhihi[k][j],WYTimlohi[k][j],
                    WYTimhilo[k][j],WYTimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_dec(&prdrehihi[i][j],&prdrelohi[i][j],
                    &prdrehilo[i][j],&prdrelolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            // accim = Qim[i][coloff+k]*WYTre[k][j]
            //       + Qre[i][coloff+k]*WYTim[k][j];
            // prdim[i][j] = prdim[i][j] + accim;
            qdf_mul(Qimhihi[i][coloff+k],Qimlohi[i][coloff+k],
                    Qimhilo[i][coloff+k],Qimlolo[i][coloff+k],
                    WYTrehihi[k][j],WYTrelohi[k][j],
                    WYTrehilo[k][j],WYTrelolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdimhihi[i][j],&prdimlohi[i][j],
                    &prdimhilo[i][j],&prdimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
            qdf_mul(Qrehihi[i][coloff+k],Qrelohi[i][coloff+k],
                    Qrehilo[i][coloff+k],Qrelolo[i][coloff+k],
                    WYTimhihi[k][j],WYTimlohi[k][j],
                    WYTimhilo[k][j],WYTimlolo[k][j],
                    &acchihi,&acclohi,&acchilo,&acclolo);
            qdf_inc(&prdimhihi[i][j],&prdimlohi[i][j],
                    &prdimhilo[i][j],&prdimlolo[i][j],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }

   if(verbose)
   {
      cout << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "QWYT[" << i << "][" << j << "]re : "
                 << prdrehihi[i][j] << "  " << prdrelohi[i][j] << endl
                 << "               "
                 << prdrehilo[i][j] << "  " << prdrelolo[i][j] << endl;
            cout << "QWYT[" << i << "][" << j << "]im : "
                 << prdimhihi[i][j] << "  " << prdimlohi[i][j] << endl
                 << "               "
                 << prdimhilo[i][j] << "  " << prdimlolo[i][j] << endl;
         }

      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "Q[" << i << "][" << coloff+j << "]re : "
                 << Qrehihi[i][coloff+j] << "  "
                 << Qrelohi[i][coloff+j] << endl
                 << "            "
                 << Qrehilo[i][coloff+j] << "  "
                 << Qrelolo[i][coloff+j] << endl;
            cout << "Q[" << i << "][" << coloff+j << "]im : "
                 << Qimhihi[i][coloff+j] << "  "
                 << Qimlohi[i][coloff+j] << endl
                 << "            "
                 << Qimhilo[i][coloff+j] << "  "
                 << Qimlolo[i][coloff+j] << endl;
         }
   }
   for(int i=0; i<dim; i++)           // Q = Q + Q*W*Y^H
      for(int j=0; j<rowdim; j++)
      {
         // Qre[i][coloff+j] = Qre[i][coloff+j] + prdre[i][j];
         qdf_inc(&Qrehihi[i][coloff+j],&Qrelohi[i][coloff+j],
                 &Qrehilo[i][coloff+j],&Qrelolo[i][coloff+j],
                 prdrehihi[i][j],prdrelohi[i][j],
                 prdrehilo[i][j],prdrelolo[i][j]);
         // Qim[i][coloff+j] = Qim[i][coloff+j] + prdim[i][j];
         qdf_inc(&Qimhihi[i][coloff+j],&Qimlohi[i][coloff+j],
                 &Qimhilo[i][coloff+j],&Qimlolo[i][coloff+j],
                 prdimhihi[i][j],prdimlohi[i][j],
                 prdimhilo[i][j],prdimlolo[i][j]);
      }

   if(verbose)
   {
      cout << "Q after the update with QWYT :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<rowdim; j++)
         {
            cout << "Q[" << i << "][" << coloff+j << "]re : "
                 << Qrehihi[i][coloff+j] << "  "
                 << Qrelohi[i][coloff+j] << endl
                 << "            "
                 << Qrehilo[i][coloff+j] << "  "
                 << Qrelolo[i][coloff+j] << endl;
            cout << "Q[" << i << "][" << coloff+j << "]im : "
                 << Qimhihi[i][coloff+j] << "  "
                 << Qimlohi[i][coloff+j] << endl
                 << "            "
                 << Qimhilo[i][coloff+j] << "  "
                 << Qimlolo[i][coloff+j] << endl;
         }
   }
   for(int i=0; i<rowdim; i++)
   {
      free(WYTrehihi[i]); free(WYTrelohi[i]);
      free(WYTrehilo[i]); free(WYTrelolo[i]);
      free(WYTimhihi[i]); free(WYTimlohi[i]);
      free(WYTimhilo[i]); free(WYTimlolo[i]);
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrehihi[i]); free(prdrelohi[i]);
      free(prdrehilo[i]); free(prdrelolo[i]);
      free(prdimhihi[i]); free(prdimlohi[i]);
      free(prdimhilo[i]); free(prdimlolo[i]);
   }
   free(WYTrehihi); free(WYTrelohi); free(WYTrehilo); free(WYTrelolo);
   free(WYTimhihi); free(WYTimlohi); free(WYTimhilo); free(WYTimlolo);
   free(prdrehihi); free(prdrelohi); free(prdrehilo); free(prdrelolo);
   free(prdimhihi); free(prdimlohi); free(prdimhilo); free(prdimlolo);
}

void CPU_dbl4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *lapsed, bool verbose )
{
   double betahihi,betalohi,betahilo,betalolo;
   double *xhihi = new double[nrows]; // input vector for house
   double *xlohi = new double[nrows];
   double *xhilo = new double[nrows];
   double *xlolo = new double[nrows];
   double *vhihi = new double[nrows]; // Householder vector
   double *vlohi = new double[nrows];
   double *vhilo = new double[nrows];
   double *vlolo = new double[nrows];
   double *Bhihi = new double[szt];   // the betas
   double *Blohi = new double[szt];
   double *Bhilo = new double[szt];
   double *Blolo = new double[szt];
   double **Yhihi = new double*[szt]; // Householder vectors in one block
   double **Ylohi = new double*[szt];
   double **Yhilo = new double*[szt];
   double **Ylolo = new double*[szt];
   double **Whihi = new double*[szt]; // columns of W
   double **Wlohi = new double*[szt];
   double **Whilo = new double*[szt];
   double **Wlolo = new double*[szt];

   for(int j=0; j<szt; j++)
   {
      Yhihi[j] = new double[nrows];
      Ylohi[j] = new double[nrows];
      Yhilo[j] = new double[nrows];
      Ylolo[j] = new double[nrows];
      Whihi[j] = new double[nrows];
      Wlohi[j] = new double[nrows];
      Whilo[j] = new double[nrows];
      Wlolo[j] = new double[nrows];
   }
   clock_t start = clock();

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qhihi[i][j] = 0.0;
         Qlohi[i][j] = 0.0;
         Qhilo[i][j] = 0.0;
         Qlolo[i][j] = 0.0;
      }
      Qhihi[i][i] = 1.0;

      for(int j=0; j<ncols; j++) 
      {
         Rhihi[i][j] = Ahihi[i][j];
         Rlohi[i][j] = Alohi[i][j];
         Rhilo[i][j] = Ahilo[i][j];
         Rlolo[i][j] = Alolo[i][j];
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
            xhihi[i-colidx] = Rhihi[i][colidx];
            xlohi[i-colidx] = Rlohi[i][colidx];
            xhilo[i-colidx] = Rhilo[i][colidx];
            xlolo[i-colidx] = Rlolo[i][colidx];
         }
         CPU_dbl4_factors_house
            (nrowscol,xhihi,    xlohi,    xhilo,    xlolo,
                      vhihi,    vlohi,    vhilo,    vlolo,
                  &betahihi,&betalohi,&betahilo,&betalolo);

         if(verbose)
         {
            cout << "beta[" << colidx << "] : "
                 << betahihi << "  " << betalohi << endl
                 << "          "
                 << betahilo << "  " << betalolo << endl;
            for(int i=colidx; i<nrows; i++)
               cout << "v[" << i-colidx << "] : "
                    << vhihi[i-colidx] << "  " << vlohi[i-colidx] << endl
                    << "       "
                    << vhilo[i-colidx] << "  " << vlolo[i-colidx] << endl;
            cout << "the R matrix :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : "
                       << Rhihi[i][j] << "  " << Rlohi[i][j] << endl
                       << "          "
                       << Rhilo[i][j] << "  " << Rlolo[i][j] << endl;
         }
         CPU_dbl4_factors_leftRupdate
            (nrows,endcol,colidx,Rhihi,   Rlohi,   Rhilo,   Rlolo,
                                 vhihi,   vlohi,   vhilo,   vlolo,
                              betahihi,betalohi,betahilo,betalolo);

         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
                  cout << "R[" << i << "][" << j << "] : "
                       << Rhihi[i][j] << "  " << Rlohi[i][j] << endl
                       << "           "
                       << Rhilo[i][j] << "  " << Rlolo[i][j] << endl;
         }
         Bhihi[L] = betahihi;
         Blohi[L] = betalohi;
         Bhilo[L] = betahilo;
         Blolo[L] = betalolo;

         for(int i=0; i<L; i++)
         {
            Yhihi[L][i] = 0.0;
            Ylohi[L][i] = 0.0;
            Yhilo[L][i] = 0.0;
            Ylolo[L][i] = 0.0;
         }
         for(int i=0; i<nrowscol; i++)
         {
            Yhihi[L][i+L] = vhihi[i];
            Ylohi[L][i+L] = vlohi[i];
            Yhilo[L][i+L] = vhilo[i];
            Ylolo[L][i+L] = vlolo[i];
         }
      }
      CPU_dbl4_blocked_VB_to_W
         (nrows-k*szt,szt,Bhihi,Blohi,Bhilo,Blolo,
                          Yhihi,Ylohi,Yhilo,Ylolo,
                          Whihi,Wlohi,Whilo,Wlolo);

      if(k<nbt-1)
      {
         CPU_dbl4_blocked_leftRupdate
            (nrows,ncols,szt,k,Rhihi,Rlohi,Rhilo,Rlolo,
                               Yhihi,Ylohi,Yhilo,Ylolo,
                               Whihi,Wlohi,Whilo,Wlolo,verbose);
      }
      CPU_dbl4_blocked_rightQupdate
         (nrows,szt,k,Qhihi,Qlohi,Qhilo,Qlolo,
                      Yhihi,Ylohi,Yhilo,Ylolo,
                      Whihi,Wlohi,Whilo,Wlolo,verbose);
   }
   clock_t end = clock();
   *lapsed = double(end - start)/CLOCKS_PER_SEC;

   free(xhihi); free(xlohi); free(xhilo); free(xlolo);
   free(vhihi); free(vlohi); free(vhilo); free(vlolo); 
   free(Bhihi); free(Blohi); free(Bhilo); free(Blolo);

   for(int j=0; j<szt; j++)
   {
      free(Yhihi[j]); free(Ylohi[j]); free(Yhilo[j]); free(Ylolo[j]);
      free(Whihi[j]); free(Wlohi[j]); free(Whilo[j]); free(Wlolo[j]);
   }
   free(Yhihi); free(Ylohi); free(Yhilo); free(Ylolo);
   free(Whihi); free(Wlohi); free(Whilo); free(Wlolo);
}

void CPU_cmplx4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *lapsec, bool verbose )
{
   double betahihi,betalohi,betahilo,betalolo;
   double *xrehihi = new double[nrows]; // real input vector for house
   double *xrelohi = new double[nrows];
   double *xrehilo = new double[nrows];
   double *xrelolo = new double[nrows];
   double *ximhihi = new double[nrows]; // imaginary input vector for house
   double *ximlohi = new double[nrows];
   double *ximhilo = new double[nrows];
   double *ximlolo = new double[nrows];
   double *vrehihi = new double[nrows]; // real parts of a Householder vector
   double *vrelohi = new double[nrows];
   double *vrehilo = new double[nrows];
   double *vrelolo = new double[nrows];
   double *vimhihi = new double[nrows]; // imag parts of a Householder vector
   double *vimlohi = new double[nrows];
   double *vimhilo = new double[nrows];
   double *vimlolo = new double[nrows];
   double *Bhihi = new double[szt];     // the betas
   double *Blohi = new double[szt];
   double *Bhilo = new double[szt];
   double *Blolo = new double[szt];
   double **Yrehihi = new double*[szt]; // Householder vectors in one block
   double **Yrelohi = new double*[szt];
   double **Yrehilo = new double*[szt];
   double **Yrelolo = new double*[szt];
   double **Yimhihi = new double*[szt]; // imag parts of Householder vectors
   double **Yimlohi = new double*[szt];
   double **Yimhilo = new double*[szt];
   double **Yimlolo = new double*[szt];
   double **Wrehihi = new double*[szt]; // real parts of the columns of W
   double **Wrelohi = new double*[szt];
   double **Wrehilo = new double*[szt]; 
   double **Wrelolo = new double*[szt]; 
   double **Wimhihi = new double*[szt]; // imaginary parts of the columns of W
   double **Wimlohi = new double*[szt];
   double **Wimhilo = new double*[szt];
   double **Wimlolo = new double*[szt];

   for(int j=0; j<szt; j++)
   {
      Yrehihi[j] = new double[nrows];
      Yrelohi[j] = new double[nrows];
      Yrehilo[j] = new double[nrows];
      Yrelolo[j] = new double[nrows];
      Yimhihi[j] = new double[nrows];
      Yimlohi[j] = new double[nrows];
      Yimhilo[j] = new double[nrows];
      Yimlolo[j] = new double[nrows];
      Wrehihi[j] = new double[nrows];
      Wrelohi[j] = new double[nrows];
      Wrehilo[j] = new double[nrows];
      Wrelolo[j] = new double[nrows];
      Wimhihi[j] = new double[nrows];
      Wimlohi[j] = new double[nrows];
      Wimhilo[j] = new double[nrows];
      Wimlolo[j] = new double[nrows];
   }
   clock_t start = clock();

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++)
      {
         Qrehihi[i][j] = 0.0; Qrelohi[i][j] = 0.0;
         Qrehilo[i][j] = 0.0; Qrelolo[i][j] = 0.0;
         Qimhihi[i][j] = 0.0; Qimlohi[i][j] = 0.0;
         Qimhilo[i][j] = 0.0; Qimlolo[i][j] = 0.0;
      }
      Qrehihi[i][i] = 1.0;

      for(int j=0; j<ncols; j++)
      {
         Rrehihi[i][j] = Arehihi[i][j]; Rrelohi[i][j] = Arelohi[i][j];
         Rrehilo[i][j] = Arehilo[i][j]; Rrelolo[i][j] = Arelolo[i][j];
         Rimhihi[i][j] = Aimhihi[i][j]; Rimlohi[i][j] = Aimlohi[i][j];
         Rimhilo[i][j] = Aimhilo[i][j]; Rimlolo[i][j] = Aimlolo[i][j];
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
            xrehihi[i-colidx] = Rrehihi[i][colidx];
            xrelohi[i-colidx] = Rrelohi[i][colidx];
            xrehilo[i-colidx] = Rrehilo[i][colidx];
            xrelolo[i-colidx] = Rrelolo[i][colidx];
            ximhihi[i-colidx] = Rimhihi[i][colidx];
            ximlohi[i-colidx] = Rimlohi[i][colidx];
            ximhilo[i-colidx] = Rimhilo[i][colidx];
            ximlolo[i-colidx] = Rimlolo[i][colidx];
         }
         CPU_cmplx4_factors_house
            (nrowscol,xrehihi,  xrelohi,  xrehilo,  xrelolo,
                      ximhihi,  ximlohi,  ximhilo,  ximlolo,
                      vrehihi,  vrelohi,  vrehilo,  vrelolo,
                      vimhihi,  vimlohi,  vimhilo,  vimlolo,
                    &betahihi,&betalohi,&betahilo,&betalolo);
         if(verbose)
         {
            cout << "beta[" << colidx << "] : "
                 << betahihi << "  " << betalohi << endl
                 << "          "
                 << betahilo << "  " << betalolo << endl;

            for(int i=colidx; i<nrows; i++)
            {
               cout << "v[" << i-colidx << "]re : "
                    << vrehihi[i-colidx] << "  " << vrelohi[i-colidx] << endl
                    << "         "
                    << vrehilo[i-colidx] << "  " << vrelolo[i-colidx] << endl;
               cout << "v[" << i-colidx << "]im : "
                    << vimhihi[i-colidx] << "  " << vimlohi[i-colidx] << endl
                    << "         "
                    << vimhilo[i-colidx] << "  " << vimlolo[i-colidx] << endl;
            }
         }
         CPU_cmplx4_factors_leftRupdate
            (nrows,endcol,colidx,Rrehihi, Rrelohi, Rrehilo, Rrelolo,
                                 Rimhihi, Rimlohi, Rimhilo, Rimlolo,
                                 vrehihi, vrelohi, vrehilo, vrelolo,
                                 vimhihi, vimlohi, vimhilo, vimlolo,
                                betahihi,betalohi,betahilo,betalolo);
         if(verbose)
         {
            cout << "the matrix after the update :" << endl;
            for(int i=0; i<nrows; i++)
               for(int j=0; j<ncols; j++)
               {
                  cout << "R[" << i << "][" << j << "]re : "
                       << Rrehihi[i][j] << "  " << Rrelohi[i][j] << endl
                       << "            "
                       << Rrehilo[i][j] << "  " << Rrelolo[i][j] << endl;
                  cout << "R[" << i << "][" << j << "]im : "
                       << Rimhihi[i][j] << "  " << Rimlohi[i][j] << endl
                       << "            "
                       << Rimhilo[i][j] << "  " << Rimlolo[i][j] << endl;
               }
         }
         Bhihi[L] = betahihi;
         Blohi[L] = betalohi;
         Bhilo[L] = betahilo;
         Blolo[L] = betalolo;

         for(int i=0; i<L; i++)
         {
            Yrehihi[L][i] = 0.0;
            Yrelohi[L][i] = 0.0;
            Yrehilo[L][i] = 0.0;
            Yrelolo[L][i] = 0.0;
            Yimhihi[L][i] = 0.0;
            Yimlohi[L][i] = 0.0;
            Yimhilo[L][i] = 0.0;
            Yimlolo[L][i] = 0.0;
         }
         for(int i=0; i<nrowscol; i++)
         {
            Yrehihi[L][i+L] = vrehihi[i];
            Yrelohi[L][i+L] = vrelohi[i];
            Yrehilo[L][i+L] = vrehilo[i];
            Yrelolo[L][i+L] = vrelolo[i];
            Yimhihi[L][i+L] = vimhihi[i];
            Yimlohi[L][i+L] = vimlohi[i];
            Yimhilo[L][i+L] = vimhilo[i];
            Yimlolo[L][i+L] = vimlolo[i];
         }
      }
      CPU_cmplx4_blocked_VB_to_W
         (nrows-k*szt,szt,Bhihi,  Blohi,  Bhilo,  Blolo,
                        Yrehihi,Yrelohi,Yrehilo,Yrelolo,
                        Yimhihi,Yimlohi,Yimhilo,Yimlolo,
                        Wrehihi,Wrelohi,Wrehilo,Wrelolo,
                        Wimhihi,Wimlohi,Wimhilo,Wimlolo);
      if(k<nbt-1)
      {
         CPU_cmplx4_blocked_leftRupdate
            (nrows,ncols,szt,k,Rrehihi,Rrelohi,Rrehilo,Rrelolo,
                               Rimhihi,Rimlohi,Rimhilo,Rimlolo,
                               Yrehihi,Yrelohi,Yrehilo,Yrelolo,
                               Yimhihi,Yimlohi,Yimhilo,Yimlolo,
                               Wrehihi,Wrelohi,Wrehilo,Wrelolo,
                               Wimhihi,Wimlohi,Wimhilo,Wimlolo,verbose);
      }
      CPU_cmplx4_blocked_rightQupdate
         (nrows,szt,k,Qrehihi,Qrelohi,Qrehilo,Qrelolo,
                      Qimhihi,Qimlohi,Qimhilo,Qimlolo,
                      Yrehihi,Yrelohi,Yrehilo,Yrelolo,
                      Yimhihi,Yimlohi,Yimhilo,Yimlolo,
                      Wrehihi,Wrelohi,Wrehilo,Wrelolo,
                      Wimhihi,Wimlohi,Wimhilo,Wimlolo,verbose);
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(xrehihi); free(xrelohi); free(xrehilo); free(xrelolo);
   free(ximhihi); free(ximlohi); free(ximhilo); free(ximlolo);
   free(vrehihi); free(vrelohi); free(vrehilo); free(vrelolo); 
   free(vimhihi); free(vimlohi); free(vimhilo); free(vimlolo);
   free(Bhihi); free(Blohi); free(Bhilo); free(Blolo);

   for(int j=0; j<szt; j++)
   {
      free(Yrehihi[j]); free(Yrelohi[j]); free(Yrehilo[j]); free(Yrelolo[j]);
      free(Yimhihi[j]); free(Yimlohi[j]); free(Yimhilo[j]); free(Yimlolo[j]);
      free(Wrehihi[j]); free(Wrelohi[j]); free(Wrehilo[j]); free(Wrelolo[j]);
      free(Wimhihi[j]); free(Wimlohi[j]); free(Wimhilo[j]); free(Wimlolo[j]);
   }
   free(Yrehihi); free(Yrelohi); free(Yrehilo); free(Yrelolo);
   free(Yimhihi); free(Yimlohi); free(Yimhilo); free(Yimlolo);
   free(Wrehihi); free(Wrelohi); free(Wrehilo); free(Wrelolo);
   free(Wimhihi); free(Wimlohi); free(Wimhilo); free(Wimlolo);
}
