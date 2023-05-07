// The file cmplx8_newton_method.cpp defines the functions with prototypes in
// the file cmplx8_newton_method.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "random_numbers.h"
#include "random_monomials.h"
#include "octo_double_functions.h"
#include "dbl8_factorizations.h"
#include "dbl8_monomial_systems.h"
#include "unimodular_matrices.h"
#include "dbl8_systems_host.h"
#include "dbl8_systems_kernels.h"
#include "dbl8_bals_host.h"
#include "dbl8_tail_kernels.h"
#include "dbl8_bals_kernels.h"
#include "dbl8_newton_testers.h"
#include "cmplx8_newton_method.h"
#include "cmplx8_path_tracker.h"

int cmplx8_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehihihi, double **mbrelohihi,
   double **mbrehilohi, double **mbrelolohi,
   double **mbrehihilo, double **mbrelohilo,
   double **mbrehilolo, double **mbrelololo,
   double **mbimhihihi, double **mbimlohihi,
   double **mbimhilohi, double **mbimlolohi,
   double **mbimhihilo, double **mbimlohilo,
   double **mbimhilolo, double **mbimlololo,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo,
   double **accrehihihi, double **accrelohihi,
   double **accrehilohi, double **accrelolohi,
   double **accrehihilo, double **accrelohilo,
   double **accrehilolo, double **accrelololo,
   double **accimhihihi, double **accimlohihi,
   double **accimhilohi, double **accimlolohi,
   double **accimhihilo, double **accimlohilo,
   double **accimhilolo, double **accimlololo,
   double **inputrehihihi_h, double **inputrelohihi_h,
   double **inputrehilohi_h, double **inputrelolohi_h,
   double **inputrehihilo_h, double **inputrelohilo_h,
   double **inputrehilolo_h, double **inputrelololo_h,
   double **inputimhihihi_h, double **inputimlohihi_h,
   double **inputimhilohi_h, double **inputimlolohi_h,
   double **inputimhihilo_h, double **inputimlohilo_h,
   double **inputimhilolo_h, double **inputimlololo_h,
   double **inputrehihihi_d, double **inputrelohihi_d,
   double **inputrehilohi_d, double **inputrelolohi_d,
   double **inputrehihilo_d, double **inputrelohilo_d,
   double **inputrehilolo_d, double **inputrelololo_d,
   double **inputimhihihi_d, double **inputimlohihi_d,
   double **inputimhilohi_d, double **inputimlolohi_d,
   double **inputimhihilo_d, double **inputimlohilo_d,
   double **inputimhilolo_d, double **inputimlololo_d,
   double ***outputrehihihi_h, double ***outputrelohihi_h,
   double ***outputrehilohi_h, double ***outputrelolohi_h,
   double ***outputrehihilo_h, double ***outputrelohilo_h,
   double ***outputrehilolo_h, double ***outputrelololo_h,
   double ***outputimhihihi_h, double ***outputimlohihi_h,
   double ***outputimhilohi_h, double ***outputimlolohi_h,
   double ***outputimhihilo_h, double ***outputimlohilo_h,
   double ***outputimhilolo_h, double ***outputimlololo_h,
   double ***outputrehihihi_d, double ***outputrelohihi_d,
   double ***outputrehilohi_d, double ***outputrelolohi_d,
   double ***outputrehihilo_d, double ***outputrelohilo_d,
   double ***outputrehilolo_d, double ***outputrelololo_d,
   double ***outputimhihihi_d, double ***outputimlohihi_d,
   double ***outputimhilohi_d, double ***outputimlolohi_d,
   double ***outputimhihilo_d, double ***outputimlohilo_d,
   double ***outputimhilolo_d, double ***outputimlololo_d,
   double **funvalrehihihi_h, double **funvalrelohihi_h,
   double **funvalrehilohi_h, double **funvalrelolohi_h,
   double **funvalrehihilo_h, double **funvalrelohilo_h,
   double **funvalrehilolo_h, double **funvalrelololo_h,
   double **funvalimhihihi_h, double **funvalimlohihi_h,
   double **funvalimhilohi_h, double **funvalimlolohi_h,
   double **funvalimhihilo_h, double **funvalimlohilo_h,
   double **funvalimhilolo_h, double **funvalimlololo_h,
   double **funvalrehihihi_d, double **funvalrelohihi_d,
   double **funvalrehilohi_d, double **funvalrelolohi_d,
   double **funvalrehihilo_d, double **funvalrelohilo_d,
   double **funvalrehilolo_d, double **funvalrelololo_d,
   double **funvalimhihihi_d, double **funvalimlohihi_d,
   double **funvalimhilohi_d, double **funvalimlolohi_d,
   double **funvalimhihilo_d, double **funvalimlohilo_d,
   double **funvalimhilolo_d, double **funvalimlololo_d,
   double ***jacvalrehihihi_h, double ***jacvalrelohihi_h,
   double ***jacvalrehilohi_h, double ***jacvalrelolohi_h,
   double ***jacvalrehihilo_h, double ***jacvalrelohilo_h,
   double ***jacvalrehilolo_h, double ***jacvalrelololo_h,
   double ***jacvalimhihihi_h, double ***jacvalimlohihi_h,
   double ***jacvalimhilohi_h, double ***jacvalimlolohi_h,
   double ***jacvalimhihilo_h, double ***jacvalimlohilo_h,
   double ***jacvalimhilolo_h, double ***jacvalimlololo_h,
   double ***jacvalrehihihi_d, double ***jacvalrelohihi_d,
   double ***jacvalrehilohi_d, double ***jacvalrelolohi_d,
   double ***jacvalrehihilo_d, double ***jacvalrelohilo_d,
   double ***jacvalrehilolo_d, double ***jacvalrelololo_d,
   double ***jacvalimhihihi_d, double ***jacvalimlohihi_d,
   double ***jacvalimhilohi_d, double ***jacvalimlolohi_d,
   double ***jacvalimhihilo_d, double ***jacvalimlohilo_d,
   double ***jacvalimhilolo_d, double ***jacvalimlololo_d,
   double **rhsrehihihi_h, double **rhsrelohihi_h,
   double **rhsrehilohi_h, double **rhsrelolohi_h,
   double **rhsrehihilo_h, double **rhsrelohilo_h,
   double **rhsrehilolo_h, double **rhsrelololo_h,
   double **rhsimhihihi_h, double **rhsimlohihi_h,
   double **rhsimhilohi_h, double **rhsimlolohi_h,
   double **rhsimhihilo_h, double **rhsimlohilo_h,
   double **rhsimhilolo_h, double **rhsimlololo_h,
   double **rhsrehihihi_d, double **rhsrelohihi_d, 
   double **rhsrehilohi_d, double **rhsrelolohi_d, 
   double **rhsrehihilo_d, double **rhsrelohilo_d, 
   double **rhsrehilolo_d, double **rhsrelololo_d, 
   double **rhsimhihihi_d, double **rhsimlohihi_d,
   double **rhsimhilohi_d, double **rhsimlolohi_d,
   double **rhsimhihilo_d, double **rhsimlohilo_d,
   double **rhsimhilolo_d, double **rhsimlololo_d,
   double **urhsrehihihi_h, double **urhsrelohihi_h,
   double **urhsrehilohi_h, double **urhsrelolohi_h,
   double **urhsrehihilo_h, double **urhsrelohilo_h,
   double **urhsrehilolo_h, double **urhsrelololo_h,
   double **urhsimhihihi_h, double **urhsimlohihi_h,
   double **urhsimhilohi_h, double **urhsimlolohi_h,
   double **urhsimhihilo_h, double **urhsimlohilo_h,
   double **urhsimhilolo_h, double **urhsimlololo_h,
   double **urhsrehihihi_d, double **urhsrelohihi_d,
   double **urhsrehilohi_d, double **urhsrelolohi_d,
   double **urhsrehihilo_d, double **urhsrelohilo_d,
   double **urhsrehilolo_d, double **urhsrelololo_d,
   double **urhsimhihihi_d, double **urhsimlohihi_d,
   double **urhsimhilohi_d, double **urhsimlolohi_d,
   double **urhsimhihilo_d, double **urhsimlohilo_d,
   double **urhsimhilolo_d, double **urhsimlololo_d,
   double **solrehihihi_h, double **solrelohihi_h,
   double **solrehilohi_h, double **solrelolohi_h,
   double **solrehihilo_h, double **solrelohilo_h,
   double **solrehilolo_h, double **solrelololo_h,
   double **solimhihihi_h, double **solimlohihi_h, 
   double **solimhilohi_h, double **solimlolohi_h, 
   double **solimhihilo_h, double **solimlohilo_h, 
   double **solimhilolo_h, double **solimlololo_h, 
   double **solrehihihi_d, double **solrelohihi_d, 
   double **solrehilohi_d, double **solrelolohi_d, 
   double **solrehihilo_d, double **solrelohilo_d, 
   double **solrehilolo_d, double **solrelololo_d, 
   double **solimhihihi_d, double **solimlohihi_d,
   double **solimhilohi_d, double **solimlolohi_d,
   double **solimhihilo_d, double **solimlohilo_d,
   double **solimhilolo_d, double **solimlololo_d,
   double **Qrehihihi_h, double **Qrelohihi_h,
   double **Qrehilohi_h, double **Qrelolohi_h,
   double **Qrehihilo_h, double **Qrelohilo_h,
   double **Qrehilolo_h, double **Qrelololo_h,
   double **Qimhihihi_h, double **Qimlohihi_h,
   double **Qimhilohi_h, double **Qimlolohi_h,
   double **Qimhihilo_h, double **Qimlohilo_h,
   double **Qimhilolo_h, double **Qimlololo_h,
   double **Qrehihihi_d, double **Qrelohihi_d,
   double **Qrehilohi_d, double **Qrelolohi_d,
   double **Qrehihilo_d, double **Qrelohilo_d,
   double **Qrehilolo_d, double **Qrelololo_d,
   double **Qimhihihi_d, double **Qimlohihi_d, 
   double **Qimhilohi_d, double **Qimlolohi_d, 
   double **Qimhihilo_d, double **Qimlohilo_d, 
   double **Qimhilolo_d, double **Qimlololo_d, 
   double **Rrehihihi_h, double **Rrelohihi_h,
   double **Rrehilohi_h, double **Rrelolohi_h,
   double **Rrehihilo_h, double **Rrelohilo_h,
   double **Rrehilolo_h, double **Rrelololo_h,
   double **Rimhihihi_h, double **Rimlohihi_h, 
   double **Rimhilohi_h, double **Rimlolohi_h, 
   double **Rimhihilo_h, double **Rimlohilo_h, 
   double **Rimhilolo_h, double **Rimlololo_h, 
   double **Rrehihihi_d, double **Rrelohihi_d,
   double **Rrehilohi_d, double **Rrelolohi_d,
   double **Rrehihilo_d, double **Rrelohilo_d,
   double **Rrehilolo_d, double **Rrelololo_d,
   double **Rimhihihi_d, double **Rimlohihi_d,
   double **Rimhilohi_d, double **Rimlolohi_d,
   double **Rimhihilo_d, double **Rimlohilo_d,
   double **Rimhilolo_d, double **Rimlololo_d,
   double *workvecrehihihi, double *workvecrelohihi,
   double *workvecrehilohi, double *workvecrelolohi,
   double *workvecrehihilo, double *workvecrelohilo,
   double *workvecrehilolo, double *workvecrelololo,
   double *workvecimhihihi, double *workvecimlohihi,
   double *workvecimhilohi, double *workvecimlolohi,
   double *workvecimhihilo, double *workvecimlohilo,
   double *workvecimhilolo, double *workvecimlololo,
   double **resvecrehihihi, double **resvecrelohihi,
   double **resvecrehilohi, double **resvecrelolohi,
   double **resvecrehihilo, double **resvecrelohilo,
   double **resvecrehilolo, double **resvecrelololo,
   double **resvecimhihihi, double **resvecimlohihi,
   double **resvecimhilohi, double **resvecimlolohi,
   double **resvecimhihilo, double **resvecimlohilo,
   double **resvecimhilolo, double **resvecimlololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;
   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool zeroQ_h = true;
   bool zeroQ_d = true;
   bool noqr_h = false;
   bool noqr_d = false;
   int tailidx_h = 1;
   int tailidx_d = 1;
   int wrkdeg = 0; // working degree of precision
   int stepcnt = 0;

   double totcnvlapsedms = 0.0;
   double totqrlapsedms = 0.0;
   double totqtblapsedms = 0.0;
   double totbslapsedms = 0.0;
   double totupdlapsedms = 0.0;
   double totreslapsedms = 0.0;

   struct timeval begintime,endtime; // wall clock time of computations
   gettimeofday(&begintime,0);

   double dpr = 1.0; // dummy value

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step
              << " at degree " << wrkdeg << " ***" << endl;

      cmplx8_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,
          &tailidx_h,&tailidx_d,nvr,idx,exp,nbrfac,expfac,
          mbrehihihi,mbrelohihi,mbrehilohi,mbrelolohi,
          mbrehihilo,mbrelohilo,mbrehilolo,mbrelololo,
          mbimhihihi,mbimlohihi,mbimhilohi,mbimlolohi,
          mbimhihilo,mbimlohilo,mbimhilolo,mbimlololo,dpr,
          cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
          cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
          cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
          cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
          accrehihihi,accrelohihi,accrehilohi,accrelolohi,
          accrehihilo,accrelohilo,accrehilolo,accrelololo,
          accimhihihi,accimlohihi,accimhilohi,accimlolohi,
          accimhihilo,accimlohilo,accimhilolo,accimlololo,
          inputrehihihi_h,inputrelohihi_h,inputrehilohi_h,inputrelolohi_h,
          inputrehihilo_h,inputrelohilo_h,inputrehilolo_h,inputrelololo_h,
          inputimhihihi_h,inputimlohihi_h,inputimhilohi_h,inputimlolohi_h,
          inputimhihilo_h,inputimlohilo_h,inputimhilolo_h,inputimlololo_h,
          inputrehihihi_d,inputrelohihi_d,inputrehilohi_d,inputrelolohi_d,
          inputrehihilo_d,inputrelohilo_d,inputrehilolo_d,inputrelololo_d,
          inputimhihihi_d,inputimlohihi_d,inputimhilohi_d,inputimlolohi_d,
          inputimhihilo_d,inputimlohilo_d,inputimhilolo_d,inputimlololo_d,
          outputrehihihi_h,outputrelohihi_h,outputrehilohi_h,outputrelolohi_h,
          outputrehihilo_h,outputrelohilo_h,outputrehilolo_h,outputrelololo_h,
          outputimhihihi_h,outputimlohihi_h,outputimhilohi_h,outputimlolohi_h,
          outputimhihilo_h,outputimlohilo_h,outputimhilolo_h,outputimlololo_h,
          outputrehihihi_d,outputrelohihi_d,outputrehilohi_d,outputrelolohi_d,
          outputrehihilo_d,outputrelohilo_d,outputrehilolo_d,outputrelololo_d,
          outputimhihihi_d,outputimlohihi_d,outputimhilohi_d,outputimlolohi_d,
          outputimhihilo_d,outputimlohilo_d,outputimhilolo_d,outputimlololo_d,
          funvalrehihihi_h,funvalrelohihi_h,funvalrehilohi_h,funvalrelolohi_h,
          funvalrehihilo_h,funvalrelohilo_h,funvalrehilolo_h,funvalrelololo_h,
          funvalimhihihi_h,funvalimlohihi_h,funvalimhilohi_h,funvalimlolohi_h,
          funvalimhihilo_h,funvalimlohilo_h,funvalimhilolo_h,funvalimlololo_h,
          funvalrehihihi_d,funvalrelohihi_d,funvalrehilohi_d,funvalrelolohi_d,
          funvalrehihilo_d,funvalrelohilo_d,funvalrehilolo_d,funvalrelololo_d,
          funvalimhihihi_d,funvalimlohihi_d,funvalimhilohi_d,funvalimlolohi_d,
          funvalimhihilo_d,funvalimlohilo_d,funvalimhilolo_d,funvalimlololo_d,
          jacvalrehihihi_h,jacvalrelohihi_h,jacvalrehilohi_h,jacvalrelolohi_h,
          jacvalrehihilo_h,jacvalrelohilo_h,jacvalrehilolo_h,jacvalrelololo_h,
          jacvalimhihihi_h,jacvalimlohihi_h,jacvalimhilohi_h,jacvalimlolohi_h,
          jacvalimhihilo_h,jacvalimlohilo_h,jacvalimhilolo_h,jacvalimlololo_h,
          jacvalrehihihi_d,jacvalrelohihi_d,jacvalrehilohi_d,jacvalrelolohi_d,
          jacvalrehihilo_d,jacvalrelohilo_d,jacvalrehilolo_d,jacvalrelololo_d,
          jacvalimhihihi_d,jacvalimlohihi_d,jacvalimhilohi_d,jacvalimlolohi_d,
          jacvalimhihilo_d,jacvalimlohilo_d,jacvalimhilolo_d,jacvalimlololo_d,
          rhsrehihihi_h,rhsrelohihi_h,rhsrehilohi_h,rhsrelolohi_h,
          rhsrehihilo_h,rhsrelohilo_h,rhsrehilolo_h,rhsrelololo_h,
          rhsimhihihi_h,rhsimlohihi_h,rhsimhilohi_h,rhsimlolohi_h,
          rhsimhihilo_h,rhsimlohilo_h,rhsimhilolo_h,rhsimlololo_h,
          rhsrehihihi_d,rhsrelohihi_d,rhsrehilohi_d,rhsrelolohi_d,
          rhsrehihilo_d,rhsrelohilo_d,rhsrehilolo_d,rhsrelololo_d,
          rhsimhihihi_d,rhsimlohihi_d,rhsimhilohi_d,rhsimlolohi_d,
          rhsimhihilo_d,rhsimlohilo_d,rhsimhilolo_d,rhsimlololo_d,
          urhsrehihihi_h,urhsrelohihi_h,urhsrehilohi_h,urhsrelolohi_h,
          urhsrehihilo_h,urhsrelohilo_h,urhsrehilolo_h,urhsrelololo_h,
          urhsimhihihi_h,urhsimlohihi_h,urhsimhilohi_h,urhsimlolohi_h,
          urhsimhihilo_h,urhsimlohilo_h,urhsimhilolo_h,urhsimlololo_h,
          urhsrehihihi_d,urhsrelohihi_d,urhsrehilohi_d,urhsrelolohi_d,
          urhsrehihilo_d,urhsrelohilo_d,urhsrehilolo_d,urhsrelololo_d,
          urhsimhihihi_d,urhsimlohihi_d,urhsimhilohi_d,urhsimlolohi_d,
          urhsimhihilo_d,urhsimlohilo_d,urhsimhilolo_d,urhsimlololo_d,
          solrehihihi_h,solrelohihi_h,solrehilohi_h,solrelolohi_h,
          solrehihilo_h,solrelohilo_h,solrehilolo_h,solrelololo_h,
          solimhihihi_h,solimlohihi_h,solimhilohi_h,solimlolohi_h,
          solimhihilo_h,solimlohilo_h,solimhilolo_h,solimlololo_h,
          solrehihihi_d,solrelohihi_d,solrehilohi_d,solrelolohi_d,
          solrehihilo_d,solrelohilo_d,solrehilolo_d,solrelololo_d,
          solimhihihi_d,solimlohihi_d,solimhilohi_d,solimlolohi_d,
          solimhihilo_d,solimlohilo_d,solimhilolo_d,solimlololo_d,
          Qrehihihi_h,Qrelohihi_h,Qrehilohi_h,Qrelolohi_h,
          Qrehihilo_h,Qrelohilo_h,Qrehilolo_h,Qrelololo_h,
          Qimhihihi_h,Qimlohihi_h,Qimhilohi_h,Qimlolohi_h,
          Qimhihilo_h,Qimlohilo_h,Qimhilolo_h,Qimlololo_h,
          Qrehihihi_d,Qrelohihi_d,Qrehilohi_d,Qrelolohi_d,
          Qrehihilo_d,Qrelohilo_d,Qrehilolo_d,Qrelololo_d,
          Qimhihihi_d,Qimlohihi_d,Qimhilohi_d,Qimlolohi_d,
          Qimhihilo_d,Qimlohilo_d,Qimhilolo_d,Qimlololo_d,
          Rrehihihi_h,Rrelohihi_h,Rrehilohi_h,Rrelolohi_h,
          Rrehihilo_h,Rrelohilo_h,Rrehilolo_h,Rrelololo_h,
          Rimhihihi_h,Rimlohihi_h,Rimhilohi_h,Rimlolohi_h,
          Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,
          Rrehihihi_d,Rrelohihi_d,Rrehilohi_d,Rrelolohi_d,
          Rrehihilo_d,Rrelohilo_d,Rrehilolo_d,Rrelololo_d,
          Rimhihihi_d,Rimlohihi_d,Rimhilohi_d,Rimlolohi_d,
          Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,
          workvecrehihihi,workvecrelohihi,workvecrehilohi,workvecrelolohi,
          workvecrehihilo,workvecrelohilo,workvecrehilolo,workvecrelololo,
          workvecimhihihi,workvecimlohihi,workvecimhilohi,workvecimlolohi,
          workvecimhihilo,workvecimlohilo,workvecimhilolo,workvecimlololo,
          resvecrehihihi,resvecrelohihi,resvecrehilohi,resvecrelolohi,
          resvecrehihilo,resvecrelohilo,resvecrehilolo,resvecrelololo,
          resvecimhihihi,resvecimlohihi,resvecimhilohi,resvecimlolohi,
          resvecimhihilo,resvecimlohilo,resvecimhilolo,resvecimlololo,
          resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
          resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,
          &zeroQ_h,&noqr_h,&zeroQ_d,&noqr_d,
          &upidx_h,&bsidx_h,&upidx_d,&bsidx_d,
          &totcnvlapsedms,&totqrlapsedms,&totqtblapsedms,&totbslapsedms,
          &totupdlapsedms,&totreslapsedms,vrblvl,mode);

      stepcnt = stepcnt + 1;

      if(vrblvl > 0)
         cout << "up_h : " << upidx_h << "  bs_h : " << bsidx_h
              << "  tail_h : " << tailidx_h
              << "  up_d : " << upidx_d << "  bs_d : " << bsidx_d
              << "  tail_d : " << tailidx_d
              << "  wdeg : " << wrkdeg << endl;

      if((mode == 1) || (mode == 2)) if(tailidx_h >= deg) break;
      if((mode == 0) || (mode == 2)) if(tailidx_d >= deg) break;

      wrkdeg = wrkdeg + 1 + wrkdeg/2;
      if(wrkdeg > deg) wrkdeg = deg;
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   double walltimesec = seconds + microseconds*1.0e-6;

   double errsum = 0.0;

   cout << scientific << setprecision(16); // just in case vrblvl == 0
   cout << "The solution series : " << endl;
   for(int j=0; j<degp1; j++)
   {
      cout << "coefficient of degree " << j << " :" << endl;
      for(int i=0; i<dim; i++)
      {
         if(mode == 0)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputrehihihi_d[i][j] << "  "
                           << inputrelohihi_d[i][j] << endl << "  "
                           << inputrehilohi_d[i][j] << "  "
                           << inputrelolohi_d[i][j] << endl << "  "
                           << inputrehihilo_d[i][j] << "  "
                           << inputrelohilo_d[i][j] << endl << "  "
                           << inputrehilolo_d[i][j] << "  "
                           << inputrelololo_d[i][j] << endl << "  "
                           << inputimhihihi_d[i][j] << "  "
                           << inputimlohihi_d[i][j] << endl << "  "
                           << inputimhilohi_d[i][j] << "  "
                           << inputimlolohi_d[i][j] << endl << "  "
                           << inputimhihilo_d[i][j] << "  "
                           << inputimlohilo_d[i][j] << endl << "  "
                           << inputimhilolo_d[i][j] << "  "
                           << inputimlololo_d[i][j] << endl;
         }
         if(mode == 1)
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputrehihihi_h[i][j] << "  "
                           << inputrelohihi_h[i][j] << endl << "  "
                           << inputrehilohi_h[i][j] << "  "
                           << inputrelolohi_h[i][j] << endl << "  "
                           << inputrehihilo_h[i][j] << "  "
                           << inputrelohilo_h[i][j] << endl << "  "
                           << inputrehilolo_h[i][j] << "  "
                           << inputrelololo_h[i][j] << endl << "  "
                           << inputimhihihi_h[i][j] << "  "
                           << inputimlohihi_h[i][j] << endl << "  "
                           << inputimhilohi_h[i][j] << "  "
                           << inputimlolohi_h[i][j] << endl << "  "
                           << inputimhihilo_h[i][j] << "  "
                           << inputimlohilo_h[i][j] << endl << "  "
                           << inputimhilolo_h[i][j] << "  "
                           << inputimlololo_h[i][j] << endl;
         }
         if(mode == 2)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputrehihihi_d[i][j] << "  "
                           << inputrelohihi_d[i][j] << endl << "  "
                           << inputrehilohi_d[i][j] << "  "
                           << inputrelolohi_d[i][j] << endl << "  "
                           << inputrehihilo_d[i][j] << "  "
                           << inputrelohilo_d[i][j] << endl << "  "
                           << inputrehilolo_d[i][j] << "  "
                           << inputrelololo_d[i][j] << endl << "  "
                           << inputimhihihi_d[i][j] << "  "
                           << inputimlohihi_d[i][j] << endl << "  "
                           << inputimhilohi_d[i][j] << "  "
                           << inputimlolohi_d[i][j] << endl << "  "
                           << inputimhihilo_d[i][j] << "  "
                           << inputimlohilo_d[i][j] << endl << "  "
                           << inputimhilolo_d[i][j] << "  "
                           << inputimlololo_d[i][j] << endl;
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputrehihihi_h[i][j] << "  "
                           << inputrelohihi_h[i][j] << endl << "  "
                           << inputrehilohi_h[i][j] << "  "
                           << inputrelolohi_h[i][j] << endl << "  "
                           << inputrehihilo_h[i][j] << "  "
                           << inputrelohilo_h[i][j] << endl << "  "
                           << inputrehilolo_h[i][j] << "  "
                           << inputrelololo_h[i][j] << endl << "  "
                           << inputimhihihi_h[i][j] << "  "
                           << inputimlohihi_h[i][j] << endl << "  "
                           << inputimhilohi_h[i][j] << "  "
                           << inputimlolohi_h[i][j] << endl << "  "
                           << inputimhihilo_h[i][j] << "  "
                           << inputimlohilo_h[i][j] << endl << "  "
                           << inputimhilolo_h[i][j] << "  "
                           << inputimlololo_h[i][j] << endl;
            errsum += abs(inputrehihihi_h[i][j] - inputrehihihi_d[i][j])
                    + abs(inputrelohihi_h[i][j] - inputrelohihi_d[i][j])
                    + abs(inputrehilohi_h[i][j] - inputrehilohi_d[i][j])
                    + abs(inputrelolohi_h[i][j] - inputrelolohi_d[i][j])
                    + abs(inputrehihilo_h[i][j] - inputrehihilo_d[i][j])
                    + abs(inputrelohilo_h[i][j] - inputrelohilo_d[i][j])
                    + abs(inputrehilolo_h[i][j] - inputrehilolo_d[i][j])
                    + abs(inputrelololo_h[i][j] - inputrelololo_d[i][j])
                    + abs(inputimhihihi_h[i][j] - inputimhihihi_d[i][j])
                    + abs(inputimlohihi_h[i][j] - inputimlohihi_d[i][j])
                    + abs(inputimhilohi_h[i][j] - inputimhilohi_d[i][j])
                    + abs(inputimlolohi_h[i][j] - inputimlolohi_d[i][j])
                    + abs(inputimhihilo_h[i][j] - inputimhihilo_d[i][j])
                    + abs(inputimlohilo_h[i][j] - inputimlohilo_d[i][j])
                    + abs(inputimhilolo_h[i][j] - inputimhilolo_d[i][j])
                    + abs(inputimlololo_h[i][j] - inputimlololo_d[i][j]);
         }
      }
   }
   if(mode == 2) cout << "error : " << errsum << endl;

   cout << "Wall clock time on all " << stepcnt << " Newton steps : ";
   cout << fixed << setprecision(3) 
        << walltimesec << " seconds." << endl;
   cout << "     Time spent by all convolution kernels : "
        << totcnvlapsedms << " milliseconds." << endl;
   cout << "  Time spent by all Householder QR kernels : "
        << totqrlapsedms << " milliseconds." << endl;
   cout << "     Time spent by all Q times rhs kernels : "
        << totqtblapsedms << " milliseconds." << endl;
   cout << "Time spent by all backsubstitution kernels : "
        << totbslapsedms << " milliseconds." << endl;
   cout << "          Time spent by all update kernels : "
        << totupdlapsedms << " milliseconds." << endl;
   cout << "        Time spent by all residual kernels : "
        << totreslapsedms << " milliseconds." << endl;

   double totkerneltime = totcnvlapsedms + totqrlapsedms + totqtblapsedms
                        + totbslapsedms + totupdlapsedms + totreslapsedms;

   cout << "           Total time spent by all kernels : "
        << totkerneltime << " milliseconds." << endl;


   return 0;
}

int test_dbl8_complex_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputrehihihi_h = new double*[dim];
   double **inputrelohihi_h = new double*[dim];
   double **inputrehilohi_h = new double*[dim];
   double **inputrelolohi_h = new double*[dim];
   double **inputrehihilo_h = new double*[dim];
   double **inputrelohilo_h = new double*[dim];
   double **inputrehilolo_h = new double*[dim];
   double **inputrelololo_h = new double*[dim];
   double **inputimhihihi_h = new double*[dim];
   double **inputimlohihi_h = new double*[dim];
   double **inputimhilohi_h = new double*[dim];
   double **inputimlolohi_h = new double*[dim];
   double **inputimhihilo_h = new double*[dim];
   double **inputimlohilo_h = new double*[dim];
   double **inputimhilolo_h = new double*[dim];
   double **inputimlololo_h = new double*[dim];
   double **inputrehihihi_d = new double*[dim];
   double **inputrelohihi_d = new double*[dim];
   double **inputrehilohi_d = new double*[dim];
   double **inputrelolohi_d = new double*[dim];
   double **inputrehihilo_d = new double*[dim];
   double **inputrelohilo_d = new double*[dim];
   double **inputrehilolo_d = new double*[dim];
   double **inputrelololo_d = new double*[dim];
   double **inputimhihihi_d = new double*[dim];
   double **inputimlohihi_d = new double*[dim];
   double **inputimhilohi_d = new double*[dim];
   double **inputimlolohi_d = new double*[dim];
   double **inputimhihilo_d = new double*[dim];
   double **inputimlohilo_d = new double*[dim];
   double **inputimhilolo_d = new double*[dim];
   double **inputimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputrehihihi_h[i] = new double[degp1];
      inputrelohihi_h[i] = new double[degp1];
      inputrehilohi_h[i] = new double[degp1];
      inputrelolohi_h[i] = new double[degp1];
      inputrehihilo_h[i] = new double[degp1];
      inputrelohilo_h[i] = new double[degp1];
      inputrehilolo_h[i] = new double[degp1];
      inputrelololo_h[i] = new double[degp1];
      inputimhihihi_h[i] = new double[degp1];
      inputimlohihi_h[i] = new double[degp1];
      inputimhilohi_h[i] = new double[degp1];
      inputimlolohi_h[i] = new double[degp1];
      inputimhihilo_h[i] = new double[degp1];
      inputimlohilo_h[i] = new double[degp1];
      inputimhilolo_h[i] = new double[degp1];
      inputimlololo_h[i] = new double[degp1];
      inputrehihihi_d[i] = new double[degp1];
      inputrelohihi_d[i] = new double[degp1];
      inputrehilohi_d[i] = new double[degp1];
      inputrelolohi_d[i] = new double[degp1];
      inputrehihilo_d[i] = new double[degp1];
      inputrelohilo_d[i] = new double[degp1];
      inputrehilolo_d[i] = new double[degp1];
      inputrelololo_d[i] = new double[degp1];
      inputimhihihi_d[i] = new double[degp1];
      inputimlohihi_d[i] = new double[degp1];
      inputimhilohi_d[i] = new double[degp1];
      inputimlolohi_d[i] = new double[degp1];
      inputimhihilo_d[i] = new double[degp1];
      inputimlohilo_d[i] = new double[degp1];
      inputimhilolo_d[i] = new double[degp1];
      inputimlololo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double **accrehihihi = new double*[dim+1]; // accumulated power series
   double **accrelohihi = new double*[dim+1];
   double **accrehilohi = new double*[dim+1];
   double **accrelolohi = new double*[dim+1];
   double **accrehihilo = new double*[dim+1];
   double **accrelohilo = new double*[dim+1];
   double **accrehilolo = new double*[dim+1];
   double **accrelololo = new double*[dim+1];
   double **accimhihihi = new double*[dim+1];
   double **accimlohihi = new double*[dim+1];
   double **accimhilohi = new double*[dim+1];
   double **accimlolohi = new double*[dim+1];
   double **accimhihilo = new double*[dim+1];
   double **accimlohilo = new double*[dim+1];
   double **accimhilolo = new double*[dim+1];
   double **accimlololo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehihihi[i] = new double[degp1];
      accrelohihi[i] = new double[degp1];
      accrehilohi[i] = new double[degp1];
      accrelolohi[i] = new double[degp1];
      accrehihilo[i] = new double[degp1];
      accrelohilo[i] = new double[degp1];
      accrehilolo[i] = new double[degp1];
      accrelololo[i] = new double[degp1];
      accimhihihi[i] = new double[degp1];
      accimlohihi[i] = new double[degp1];
      accimhilohi[i] = new double[degp1];
      accimlolohi[i] = new double[degp1];
      accimhihilo[i] = new double[degp1];
      accimlohilo[i] = new double[degp1];
      accimhilolo[i] = new double[degp1];
      accimlololo[i] = new double[degp1];
   }
   double ***cffrehihihi = new double**[nbrcol]; // coefficients of monomials
   double ***cffrelohihi = new double**[nbrcol];
   double ***cffrehilohi = new double**[nbrcol];
   double ***cffrelolohi = new double**[nbrcol];
   double ***cffrehihilo = new double**[nbrcol];
   double ***cffrelohilo = new double**[nbrcol];
   double ***cffrehilolo = new double**[nbrcol];
   double ***cffrelololo = new double**[nbrcol];
   double ***cffimhihihi = new double**[nbrcol];
   double ***cffimlohihi = new double**[nbrcol];
   double ***cffimhilohi = new double**[nbrcol];
   double ***cffimlolohi = new double**[nbrcol];
   double ***cffimhihilo = new double**[nbrcol];
   double ***cffimlohilo = new double**[nbrcol];
   double ***cffimhilolo = new double**[nbrcol];
   double ***cffimlololo = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffrehihihi[i] = new double*[dim];
      cffrelohihi[i] = new double*[dim];
      cffrehilohi[i] = new double*[dim];
      cffrelolohi[i] = new double*[dim];
      cffrehihilo[i] = new double*[dim];
      cffrelohilo[i] = new double*[dim];
      cffrehilolo[i] = new double*[dim];
      cffrelololo[i] = new double*[dim];
      cffimhihihi[i] = new double*[dim];
      cffimlohihi[i] = new double*[dim];
      cffimhilohi[i] = new double*[dim];
      cffimlolohi[i] = new double*[dim];
      cffimhihilo[i] = new double*[dim];
      cffimlohilo[i] = new double*[dim];
      cffimhilolo[i] = new double*[dim];
      cffimlololo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffrehihihi[i][j] = new double[degp1];
         cffrelohihi[i][j] = new double[degp1];
         cffrehilohi[i][j] = new double[degp1];
         cffrelolohi[i][j] = new double[degp1];
         cffrehihilo[i][j] = new double[degp1];
         cffrelohilo[i][j] = new double[degp1];
         cffrehilolo[i][j] = new double[degp1];
         cffrelololo[i][j] = new double[degp1];
         cffimhihihi[i][j] = new double[degp1];
         cffimlohihi[i][j] = new double[degp1];
         cffimhilohi[i][j] = new double[degp1];
         cffimlolohi[i][j] = new double[degp1];
         cffimhihilo[i][j] = new double[degp1];
         cffimlohilo[i][j] = new double[degp1];
         cffimhilolo[i][j] = new double[degp1];
         cffimlololo[i][j] = new double[degp1];

         for(int k=0; k<deg; k++)
         {
            cffrehihihi[i][j][k] = 0.0;
            cffrelohihi[i][j][k] = 0.0;
            cffrehilohi[i][j][k] = 0.0;
            cffrelolohi[i][j][k] = 0.0;
            cffrehihilo[i][j][k] = 0.0;
            cffrelohilo[i][j][k] = 0.0;
            cffrehilolo[i][j][k] = 0.0;
            cffrelololo[i][j][k] = 0.0;
            cffimhihihi[i][j][k] = 0.0;
            cffimlohihi[i][j][k] = 0.0;
            cffimhilohi[i][j][k] = 0.0;
            cffimlolohi[i][j][k] = 0.0;
            cffimhihilo[i][j][k] = 0.0;
            cffimlohilo[i][j][k] = 0.0;
            cffimhilolo[i][j][k] = 0.0;
            cffimlololo[i][j][k] = 0.0;
         }
      }
   }
   double ***outputrehihihi_h;
   double ***outputrelohihi_h;
   double ***outputrehilohi_h;
   double ***outputrelolohi_h;
   double ***outputrehihilo_h;
   double ***outputrelohilo_h;
   double ***outputrehilolo_h;
   double ***outputrelololo_h;
   double ***outputimhihihi_h;
   double ***outputimlohihi_h;
   double ***outputimhilohi_h;
   double ***outputimlolohi_h;
   double ***outputimhihilo_h;
   double ***outputimlohilo_h;
   double ***outputimhilolo_h;
   double ***outputimlololo_h;
   double ***outputrehihihi_d;
   double ***outputrelohihi_d;
   double ***outputrehilohi_d;
   double ***outputrelolohi_d;
   double ***outputrehihilo_d;
   double ***outputrelohilo_d;
   double ***outputrehilolo_d;
   double ***outputrelololo_d;
   double ***outputimhihihi_d;
   double ***outputimlohihi_d;
   double ***outputimhilohi_d;
   double ***outputimlolohi_d;
   double ***outputimhihilo_d;
   double ***outputimlohilo_d;
   double ***outputimhilolo_d;
   double ***outputimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      outputrehihihi_h = new double**[dim];
      outputrelohihi_h = new double**[dim];
      outputrehilohi_h = new double**[dim];
      outputrelolohi_h = new double**[dim];
      outputrehihilo_h = new double**[dim];
      outputrelohilo_h = new double**[dim];
      outputrehilolo_h = new double**[dim];
      outputrelololo_h = new double**[dim];
      outputimhihihi_h = new double**[dim];
      outputimlohihi_h = new double**[dim];
      outputimhilohi_h = new double**[dim];
      outputimlolohi_h = new double**[dim];
      outputimhihilo_h = new double**[dim];
      outputimlohilo_h = new double**[dim];
      outputimhilolo_h = new double**[dim];
      outputimlololo_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputrehihihi_h[i] = new double*[dim+1];
         outputrelohihi_h[i] = new double*[dim+1];
         outputrehilohi_h[i] = new double*[dim+1];
         outputrelolohi_h[i] = new double*[dim+1];
         outputrehihilo_h[i] = new double*[dim+1];
         outputrelohilo_h[i] = new double*[dim+1];
         outputrehilolo_h[i] = new double*[dim+1];
         outputrelololo_h[i] = new double*[dim+1];
         outputimhihihi_h[i] = new double*[dim+1];
         outputimlohihi_h[i] = new double*[dim+1];
         outputimhilohi_h[i] = new double*[dim+1];
         outputimlolohi_h[i] = new double*[dim+1];
         outputimhihilo_h[i] = new double*[dim+1];
         outputimlohilo_h[i] = new double*[dim+1];
         outputimhilolo_h[i] = new double*[dim+1];
         outputimlololo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputrehihihi_h[i][j] = new double[degp1];
            outputrelohihi_h[i][j] = new double[degp1];
            outputrehilohi_h[i][j] = new double[degp1];
            outputrelolohi_h[i][j] = new double[degp1];
            outputrehihilo_h[i][j] = new double[degp1];
            outputrelohilo_h[i][j] = new double[degp1];
            outputrehilolo_h[i][j] = new double[degp1];
            outputrelololo_h[i][j] = new double[degp1];
            outputimhihihi_h[i][j] = new double[degp1];
            outputimlohihi_h[i][j] = new double[degp1];
            outputimhilohi_h[i][j] = new double[degp1];
            outputimlolohi_h[i][j] = new double[degp1];
            outputimhihilo_h[i][j] = new double[degp1];
            outputimlohilo_h[i][j] = new double[degp1];
            outputimhilolo_h[i][j] = new double[degp1];
            outputimlololo_h[i][j] = new double[degp1];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputrehihihi_d = new double**[dim];
      outputrelohihi_d = new double**[dim];
      outputrehilohi_d = new double**[dim];
      outputrelolohi_d = new double**[dim];
      outputrehihilo_d = new double**[dim];
      outputrelohilo_d = new double**[dim];
      outputrehilolo_d = new double**[dim];
      outputrelololo_d = new double**[dim];
      outputimhihihi_d = new double**[dim];
      outputimlohihi_d = new double**[dim];
      outputimhilohi_d = new double**[dim];
      outputimlolohi_d = new double**[dim];
      outputimhihilo_d = new double**[dim];
      outputimlohilo_d = new double**[dim];
      outputimhilolo_d = new double**[dim];
      outputimlololo_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputrehihihi_d[i] = new double*[dim+1];
         outputrelohihi_d[i] = new double*[dim+1];
         outputrehilohi_d[i] = new double*[dim+1];
         outputrelolohi_d[i] = new double*[dim+1];
         outputrehihilo_d[i] = new double*[dim+1];
         outputrelohilo_d[i] = new double*[dim+1];
         outputrehilolo_d[i] = new double*[dim+1];
         outputrelololo_d[i] = new double*[dim+1];
         outputimhihihi_d[i] = new double*[dim+1];
         outputimlohihi_d[i] = new double*[dim+1];
         outputimhilohi_d[i] = new double*[dim+1];
         outputimlolohi_d[i] = new double*[dim+1];
         outputimhihilo_d[i] = new double*[dim+1];
         outputimlohilo_d[i] = new double*[dim+1];
         outputimhilolo_d[i] = new double*[dim+1];
         outputimlololo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputrehihihi_d[i][j] = new double[degp1];
            outputrelohihi_d[i][j] = new double[degp1];
            outputrehilohi_d[i][j] = new double[degp1];
            outputrelolohi_d[i][j] = new double[degp1];
            outputrehihilo_d[i][j] = new double[degp1];
            outputrelohilo_d[i][j] = new double[degp1];
            outputrehilolo_d[i][j] = new double[degp1];
            outputrelololo_d[i][j] = new double[degp1];
            outputimhihihi_d[i][j] = new double[degp1];
            outputimlohihi_d[i][j] = new double[degp1];
            outputimhilohi_d[i][j] = new double[degp1];
            outputimlolohi_d[i][j] = new double[degp1];
            outputimhihilo_d[i][j] = new double[degp1];
            outputimlohilo_d[i][j] = new double[degp1];
            outputimhilolo_d[i][j] = new double[degp1];
            outputimlololo_d[i][j] = new double[degp1];
         }
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalrehihihi_h;
   double **funvalrelohihi_h;
   double **funvalrehilohi_h;
   double **funvalrelolohi_h;
   double **funvalrehihilo_h;
   double **funvalrelohilo_h;
   double **funvalrehilolo_h;
   double **funvalrelololo_h;
   double **funvalimhihihi_h;
   double **funvalimlohihi_h;
   double **funvalimhilohi_h;
   double **funvalimlolohi_h;
   double **funvalimhihilo_h;
   double **funvalimlohilo_h;
   double **funvalimhilolo_h;
   double **funvalimlololo_h;
   double **funvalrehihihi_d;
   double **funvalrelohihi_d;
   double **funvalrehilohi_d;
   double **funvalrelolohi_d;
   double **funvalrehihilo_d;
   double **funvalrelohilo_d;
   double **funvalrehilolo_d;
   double **funvalrelololo_d;
   double **funvalimhihihi_d;
   double **funvalimlohihi_d;
   double **funvalimhilohi_d;
   double **funvalimlolohi_d;
   double **funvalimhihilo_d;
   double **funvalimlohilo_d;
   double **funvalimhilolo_d;
   double **funvalimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      funvalrehihihi_h = new double*[dim];
      funvalrelohihi_h = new double*[dim];
      funvalrehilohi_h = new double*[dim];
      funvalrelolohi_h = new double*[dim];
      funvalrehihilo_h = new double*[dim];
      funvalrelohilo_h = new double*[dim];
      funvalrehilolo_h = new double*[dim];
      funvalrelololo_h = new double*[dim];
      funvalimhihihi_h = new double*[dim];
      funvalimlohihi_h = new double*[dim];
      funvalimhilohi_h = new double*[dim];
      funvalimlolohi_h = new double*[dim];
      funvalimhihilo_h = new double*[dim];
      funvalimlohilo_h = new double*[dim];
      funvalimhilolo_h = new double*[dim];
      funvalimlololo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalrehihihi_h[i] = new double[degp1];
         funvalrelohihi_h[i] = new double[degp1];
         funvalrehilohi_h[i] = new double[degp1];
         funvalrelolohi_h[i] = new double[degp1];
         funvalrehihilo_h[i] = new double[degp1];
         funvalrelohilo_h[i] = new double[degp1];
         funvalrehilolo_h[i] = new double[degp1];
         funvalrelololo_h[i] = new double[degp1];
         funvalimhihihi_h[i] = new double[degp1];
         funvalimlohihi_h[i] = new double[degp1];
         funvalimhilohi_h[i] = new double[degp1];
         funvalimlolohi_h[i] = new double[degp1];
         funvalimhihilo_h[i] = new double[degp1];
         funvalimlohilo_h[i] = new double[degp1];
         funvalimhilolo_h[i] = new double[degp1];
         funvalimlololo_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalrehihihi_d = new double*[dim];
      funvalrelohihi_d = new double*[dim];
      funvalrehilohi_d = new double*[dim];
      funvalrelolohi_d = new double*[dim];
      funvalrehihilo_d = new double*[dim];
      funvalrelohilo_d = new double*[dim];
      funvalrehilolo_d = new double*[dim];
      funvalrelololo_d = new double*[dim];
      funvalimhihihi_d = new double*[dim];
      funvalimlohihi_d = new double*[dim];
      funvalimhilohi_d = new double*[dim];
      funvalimlolohi_d = new double*[dim];
      funvalimhihilo_d = new double*[dim];
      funvalimlohilo_d = new double*[dim];
      funvalimhilolo_d = new double*[dim];
      funvalimlololo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalrehihihi_d[i] = new double[degp1];
         funvalrelohihi_d[i] = new double[degp1];
         funvalrehilohi_d[i] = new double[degp1];
         funvalrelolohi_d[i] = new double[degp1];
         funvalrehihilo_d[i] = new double[degp1];
         funvalrelohilo_d[i] = new double[degp1];
         funvalrehilolo_d[i] = new double[degp1];
         funvalrelololo_d[i] = new double[degp1];
         funvalimhihihi_d[i] = new double[degp1];
         funvalimlohihi_d[i] = new double[degp1];
         funvalimhilohi_d[i] = new double[degp1];
         funvalimlolohi_d[i] = new double[degp1];
         funvalimhihilo_d[i] = new double[degp1];
         funvalimlohilo_d[i] = new double[degp1];
         funvalimhilolo_d[i] = new double[degp1];
         funvalimlololo_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalrehihihi_h;
   double ***jacvalrelohihi_h;
   double ***jacvalrehilohi_h;
   double ***jacvalrelolohi_h;
   double ***jacvalrehihilo_h;
   double ***jacvalrelohilo_h;
   double ***jacvalrehilolo_h;
   double ***jacvalrelololo_h;
   double ***jacvalimhihihi_h;
   double ***jacvalimlohihi_h;
   double ***jacvalimhilohi_h;
   double ***jacvalimlolohi_h;
   double ***jacvalimhihilo_h;
   double ***jacvalimlohilo_h;
   double ***jacvalimhilolo_h;
   double ***jacvalimlololo_h;
   double ***jacvalrehihihi_d;
   double ***jacvalrelohihi_d;
   double ***jacvalrehilohi_d;
   double ***jacvalrelolohi_d;
   double ***jacvalrehihilo_d;
   double ***jacvalrelohilo_d;
   double ***jacvalrehilolo_d;
   double ***jacvalrelololo_d;
   double ***jacvalimhihihi_d;
   double ***jacvalimlohihi_d;
   double ***jacvalimhilohi_d;
   double ***jacvalimlolohi_d;
   double ***jacvalimhihilo_d;
   double ***jacvalimlohilo_d;
   double ***jacvalimhilolo_d;
   double ***jacvalimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      jacvalrehihihi_h = new double**[degp1];
      jacvalrelohihi_h = new double**[degp1];
      jacvalrehilohi_h = new double**[degp1];
      jacvalrelolohi_h = new double**[degp1];
      jacvalrehihilo_h = new double**[degp1];
      jacvalrelohilo_h = new double**[degp1];
      jacvalrehilolo_h = new double**[degp1];
      jacvalrelololo_h = new double**[degp1];
      jacvalimhihihi_h = new double**[degp1];
      jacvalimlohihi_h = new double**[degp1];
      jacvalimhilohi_h = new double**[degp1];
      jacvalimlolohi_h = new double**[degp1];
      jacvalimhihilo_h = new double**[degp1];
      jacvalimlohilo_h = new double**[degp1];
      jacvalimhilolo_h = new double**[degp1];
      jacvalimlololo_h = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalrehihihi_h[i] = new double*[dim];
         jacvalrelohihi_h[i] = new double*[dim];
         jacvalrehilohi_h[i] = new double*[dim];
         jacvalrelolohi_h[i] = new double*[dim];
         jacvalrehihilo_h[i] = new double*[dim];
         jacvalrelohilo_h[i] = new double*[dim];
         jacvalrehilolo_h[i] = new double*[dim];
         jacvalrelololo_h[i] = new double*[dim];
         jacvalimhihihi_h[i] = new double*[dim];
         jacvalimlohihi_h[i] = new double*[dim];
         jacvalimhilohi_h[i] = new double*[dim];
         jacvalimlolohi_h[i] = new double*[dim];
         jacvalimhihilo_h[i] = new double*[dim];
         jacvalimlohilo_h[i] = new double*[dim];
         jacvalimhilolo_h[i] = new double*[dim];
         jacvalimlololo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalrehihihi_h[i][j] = new double[dim];
            jacvalrelohihi_h[i][j] = new double[dim];
            jacvalrehilohi_h[i][j] = new double[dim];
            jacvalrelolohi_h[i][j] = new double[dim];
            jacvalrehihilo_h[i][j] = new double[dim];
            jacvalrelohilo_h[i][j] = new double[dim];
            jacvalrehilolo_h[i][j] = new double[dim];
            jacvalrelololo_h[i][j] = new double[dim];
            jacvalimhihihi_h[i][j] = new double[dim];
            jacvalimlohihi_h[i][j] = new double[dim];
            jacvalimhilohi_h[i][j] = new double[dim];
            jacvalimlolohi_h[i][j] = new double[dim];
            jacvalimhihilo_h[i][j] = new double[dim];
            jacvalimlohilo_h[i][j] = new double[dim];
            jacvalimhilolo_h[i][j] = new double[dim];
            jacvalimlololo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      jacvalrehihihi_d = new double**[degp1];
      jacvalrelohihi_d = new double**[degp1];
      jacvalrehilohi_d = new double**[degp1];
      jacvalrelolohi_d = new double**[degp1];
      jacvalrehihilo_d = new double**[degp1];
      jacvalrelohilo_d = new double**[degp1];
      jacvalrehilolo_d = new double**[degp1];
      jacvalrelololo_d = new double**[degp1];
      jacvalimhihihi_d = new double**[degp1];
      jacvalimlohihi_d = new double**[degp1];
      jacvalimhilohi_d = new double**[degp1];
      jacvalimlolohi_d = new double**[degp1];
      jacvalimhihilo_d = new double**[degp1];
      jacvalimlohilo_d = new double**[degp1];
      jacvalimhilolo_d = new double**[degp1];
      jacvalimlololo_d = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalrehihihi_d[i] = new double*[dim];
         jacvalrelohihi_d[i] = new double*[dim];
         jacvalrehilohi_d[i] = new double*[dim];
         jacvalrelolohi_d[i] = new double*[dim];
         jacvalrehihilo_d[i] = new double*[dim];
         jacvalrelohilo_d[i] = new double*[dim];
         jacvalrehilolo_d[i] = new double*[dim];
         jacvalrelololo_d[i] = new double*[dim];
         jacvalimhihihi_d[i] = new double*[dim];
         jacvalimlohihi_d[i] = new double*[dim];
         jacvalimhilohi_d[i] = new double*[dim];
         jacvalimlolohi_d[i] = new double*[dim];
         jacvalimhihilo_d[i] = new double*[dim];
         jacvalimlohilo_d[i] = new double*[dim];
         jacvalimhilolo_d[i] = new double*[dim];
         jacvalimlololo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalrehihihi_d[i][j] = new double[dim];
            jacvalrelohihi_d[i][j] = new double[dim];
            jacvalrehilohi_d[i][j] = new double[dim];
            jacvalrelolohi_d[i][j] = new double[dim];
            jacvalrehihilo_d[i][j] = new double[dim];
            jacvalrelohilo_d[i][j] = new double[dim];
            jacvalrehilolo_d[i][j] = new double[dim];
            jacvalrelololo_d[i][j] = new double[dim];
            jacvalimhihihi_d[i][j] = new double[dim];
            jacvalimlohihi_d[i][j] = new double[dim];
            jacvalimhilohi_d[i][j] = new double[dim];
            jacvalimlolohi_d[i][j] = new double[dim];
            jacvalimhihilo_d[i][j] = new double[dim];
            jacvalimlohilo_d[i][j] = new double[dim];
            jacvalimhilolo_d[i][j] = new double[dim];
            jacvalimlololo_d[i][j] = new double[dim];
         }
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solrehihihi_h;
   double **solrelohihi_h;
   double **solrehilohi_h;
   double **solrelolohi_h;
   double **solrehihilo_h;
   double **solrelohilo_h;
   double **solrehilolo_h;
   double **solrelololo_h;
   double **solimhihihi_h;
   double **solimlohihi_h;
   double **solimhilohi_h;
   double **solimlolohi_h;
   double **solimhihilo_h;
   double **solimlohilo_h;
   double **solimhilolo_h;
   double **solimlololo_h;
   double **solrehihihi_d;
   double **solrelohihi_d;
   double **solrehilohi_d;
   double **solrelolohi_d;
   double **solrehihilo_d;
   double **solrelohilo_d;
   double **solrehilolo_d;
   double **solrelololo_d;
   double **solimhihihi_d;
   double **solimlohihi_d;
   double **solimhilohi_d;
   double **solimlolohi_d;
   double **solimhihilo_d;
   double **solimlohilo_d;
   double **solimhilolo_d;
   double **solimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      solrehihihi_h = new double*[degp1];
      solrelohihi_h = new double*[degp1];
      solrehilohi_h = new double*[degp1];
      solrelolohi_h = new double*[degp1];
      solrehihilo_h = new double*[degp1];
      solrelohilo_h = new double*[degp1];
      solrehilolo_h = new double*[degp1];
      solrelololo_h = new double*[degp1];
      solimhihihi_h = new double*[degp1];
      solimlohihi_h = new double*[degp1];
      solimhilohi_h = new double*[degp1];
      solimlolohi_h = new double*[degp1];
      solimhihilo_h = new double*[degp1];
      solimlohilo_h = new double*[degp1];
      solimhilolo_h = new double*[degp1];
      solimlololo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solrehihihi_h[i] = new double[dim];
         solrelohihi_h[i] = new double[dim];
         solrehilohi_h[i] = new double[dim];
         solrelolohi_h[i] = new double[dim];
         solrehihilo_h[i] = new double[dim];
         solrelohilo_h[i] = new double[dim];
         solrehilolo_h[i] = new double[dim];
         solrelololo_h[i] = new double[dim];
         solimhihihi_h[i] = new double[dim];
         solimlohihi_h[i] = new double[dim];
         solimhilohi_h[i] = new double[dim];
         solimlolohi_h[i] = new double[dim];
         solimhihilo_h[i] = new double[dim];
         solimlohilo_h[i] = new double[dim];
         solimhilolo_h[i] = new double[dim];
         solimlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solrehihihi_d = new double*[degp1];
      solrelohihi_d = new double*[degp1];
      solrehilohi_d = new double*[degp1];
      solrelolohi_d = new double*[degp1];
      solrehihilo_d = new double*[degp1];
      solrelohilo_d = new double*[degp1];
      solrehilolo_d = new double*[degp1];
      solrelololo_d = new double*[degp1];
      solimhihihi_d = new double*[degp1];
      solimlohihi_d = new double*[degp1];
      solimhilohi_d = new double*[degp1];
      solimlolohi_d = new double*[degp1];
      solimhihilo_d = new double*[degp1];
      solimlohilo_d = new double*[degp1];
      solimhilolo_d = new double*[degp1];
      solimlololo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solrehihihi_d[i] = new double[dim];
         solrelohihi_d[i] = new double[dim];
         solrehilohi_d[i] = new double[dim];
         solrelolohi_d[i] = new double[dim];
         solrehihilo_d[i] = new double[dim];
         solrelohilo_d[i] = new double[dim];
         solrehilolo_d[i] = new double[dim];
         solrelololo_d[i] = new double[dim];
         solimhihihi_d[i] = new double[dim];
         solimlohihi_d[i] = new double[dim];
         solimhilohi_d[i] = new double[dim];
         solimlolohi_d[i] = new double[dim];
         solimhihilo_d[i] = new double[dim];
         solimlohilo_d[i] = new double[dim];
         solimhilolo_d[i] = new double[dim];
         solimlololo_d[i] = new double[dim];
      }
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhsrehihihi_h;
   double **rhsrelohihi_h;
   double **rhsrehilohi_h;
   double **rhsrelolohi_h;
   double **rhsrehihilo_h;
   double **rhsrelohilo_h;
   double **rhsrehilolo_h;
   double **rhsrelololo_h;
   double **rhsimhihihi_h;
   double **rhsimlohihi_h;
   double **rhsimhilohi_h;
   double **rhsimlolohi_h;
   double **rhsimhihilo_h;
   double **rhsimlohilo_h;
   double **rhsimhilolo_h;
   double **rhsimlololo_h;
   double **rhsrehihihi_d;
   double **rhsrelohihi_d;
   double **rhsrehilohi_d;
   double **rhsrelolohi_d;
   double **rhsrehihilo_d;
   double **rhsrelohilo_d;
   double **rhsrehilolo_d;
   double **rhsrelololo_d;
   double **rhsimhihihi_d;
   double **rhsimlohihi_d;
   double **rhsimhilohi_d;
   double **rhsimlolohi_d;
   double **rhsimhihilo_d;
   double **rhsimlohilo_d;
   double **rhsimhilolo_d;
   double **rhsimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      rhsrehihihi_h = new double*[degp1];
      rhsrelohihi_h = new double*[degp1];
      rhsrehilohi_h = new double*[degp1];
      rhsrelolohi_h = new double*[degp1];
      rhsrehihilo_h = new double*[degp1];
      rhsrelohilo_h = new double*[degp1];
      rhsrehilolo_h = new double*[degp1];
      rhsrelololo_h = new double*[degp1];
      rhsimhihihi_h = new double*[degp1];
      rhsimlohihi_h = new double*[degp1];
      rhsimhilohi_h = new double*[degp1];
      rhsimlolohi_h = new double*[degp1];
      rhsimhihilo_h = new double*[degp1];
      rhsimlohilo_h = new double*[degp1];
      rhsimhilolo_h = new double*[degp1];
      rhsimlololo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsrehihihi_h[i] = new double[dim];
         rhsrelohihi_h[i] = new double[dim];
         rhsrehilohi_h[i] = new double[dim];
         rhsrelolohi_h[i] = new double[dim];
         rhsrehihilo_h[i] = new double[dim];
         rhsrelohilo_h[i] = new double[dim];
         rhsrehilolo_h[i] = new double[dim];
         rhsrelololo_h[i] = new double[dim];
         rhsimhihihi_h[i] = new double[dim];
         rhsimlohihi_h[i] = new double[dim];
         rhsimhilohi_h[i] = new double[dim];
         rhsimlolohi_h[i] = new double[dim];
         rhsimhihilo_h[i] = new double[dim];
         rhsimlohilo_h[i] = new double[dim];
         rhsimhilolo_h[i] = new double[dim];
         rhsimlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhsrehihihi_d = new double*[degp1];
      rhsrelohihi_d = new double*[degp1];
      rhsrehilohi_d = new double*[degp1];
      rhsrelolohi_d = new double*[degp1];
      rhsrehihilo_d = new double*[degp1];
      rhsrelohilo_d = new double*[degp1];
      rhsrehilolo_d = new double*[degp1];
      rhsrelololo_d = new double*[degp1];
      rhsimhihihi_d = new double*[degp1];
      rhsimlohihi_d = new double*[degp1];
      rhsimhilohi_d = new double*[degp1];
      rhsimlolohi_d = new double*[degp1];
      rhsimhihilo_d = new double*[degp1];
      rhsimlohilo_d = new double*[degp1];
      rhsimhilolo_d = new double*[degp1];
      rhsimlololo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsrehihihi_d[i] = new double[dim];
         rhsrelohihi_d[i] = new double[dim];
         rhsrehilohi_d[i] = new double[dim];
         rhsrelolohi_d[i] = new double[dim];
         rhsrehihilo_d[i] = new double[dim];
         rhsrelohilo_d[i] = new double[dim];
         rhsrehilolo_d[i] = new double[dim];
         rhsrelololo_d[i] = new double[dim];
         rhsimhihihi_d[i] = new double[dim];
         rhsimlohihi_d[i] = new double[dim];
         rhsimhilohi_d[i] = new double[dim];
         rhsimlolohi_d[i] = new double[dim];
         rhsimhihilo_d[i] = new double[dim];
         rhsimlohilo_d[i] = new double[dim];
         rhsimhilolo_d[i] = new double[dim];
         rhsimlololo_d[i] = new double[dim];
      }
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhsrehihihi_h;
   double **urhsrelohihi_h;
   double **urhsrehilohi_h;
   double **urhsrelolohi_h;
   double **urhsrehihilo_h;
   double **urhsrelohilo_h;
   double **urhsrehilolo_h;
   double **urhsrelololo_h;
   double **urhsimhihihi_h;
   double **urhsimlohihi_h;
   double **urhsimhilohi_h;
   double **urhsimlolohi_h;
   double **urhsimhihilo_h;
   double **urhsimlohilo_h;
   double **urhsimhilolo_h;
   double **urhsimlololo_h;
   double **urhsrehihihi_d;
   double **urhsrelohihi_d;
   double **urhsrehilohi_d;
   double **urhsrelolohi_d;
   double **urhsrehihilo_d;
   double **urhsrelohilo_d;
   double **urhsrehilolo_d;
   double **urhsrelololo_d;
   double **urhsimhihihi_d;
   double **urhsimlohihi_d;
   double **urhsimhilohi_d;
   double **urhsimlolohi_d;
   double **urhsimhihilo_d;
   double **urhsimlohilo_d;
   double **urhsimhilolo_d;
   double **urhsimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      urhsrehihihi_h = new double*[degp1];
      urhsrelohihi_h = new double*[degp1];
      urhsrehilohi_h = new double*[degp1];
      urhsrelolohi_h = new double*[degp1];
      urhsrehihilo_h = new double*[degp1];
      urhsrelohilo_h = new double*[degp1];
      urhsrehilolo_h = new double*[degp1];
      urhsrelololo_h = new double*[degp1];
      urhsimhihihi_h = new double*[degp1];
      urhsimlohihi_h = new double*[degp1];
      urhsimhilohi_h = new double*[degp1];
      urhsimlolohi_h = new double*[degp1];
      urhsimhihilo_h = new double*[degp1];
      urhsimlohilo_h = new double*[degp1];
      urhsimhilolo_h = new double*[degp1];
      urhsimlololo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsrehihihi_h[i] = new double[dim];
         urhsrelohihi_h[i] = new double[dim];
         urhsrehilohi_h[i] = new double[dim];
         urhsrelolohi_h[i] = new double[dim];
         urhsrehihilo_h[i] = new double[dim];
         urhsrelohilo_h[i] = new double[dim];
         urhsrehilolo_h[i] = new double[dim];
         urhsrelololo_h[i] = new double[dim];
         urhsimhihihi_h[i] = new double[dim];
         urhsimlohihi_h[i] = new double[dim];
         urhsimhilohi_h[i] = new double[dim];
         urhsimlolohi_h[i] = new double[dim];
         urhsimhihilo_h[i] = new double[dim];
         urhsimlohilo_h[i] = new double[dim];
         urhsimhilolo_h[i] = new double[dim];
         urhsimlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhsrehihihi_d = new double*[degp1];
      urhsrelohihi_d = new double*[degp1];
      urhsrehilohi_d = new double*[degp1];
      urhsrelolohi_d = new double*[degp1];
      urhsrehihilo_d = new double*[degp1];
      urhsrelohilo_d = new double*[degp1];
      urhsrehilolo_d = new double*[degp1];
      urhsrelololo_d = new double*[degp1];
      urhsimhihihi_d = new double*[degp1];
      urhsimlohihi_d = new double*[degp1];
      urhsimhilohi_d = new double*[degp1];
      urhsimlolohi_d = new double*[degp1];
      urhsimhihilo_d = new double*[degp1];
      urhsimlohilo_d = new double*[degp1];
      urhsimhilolo_d = new double*[degp1];
      urhsimlololo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsrehihihi_d[i] = new double[dim];
         urhsrelohihi_d[i] = new double[dim];
         urhsrehilohi_d[i] = new double[dim];
         urhsrelolohi_d[i] = new double[dim];
         urhsrehihilo_d[i] = new double[dim];
         urhsrelohilo_d[i] = new double[dim];
         urhsrehilolo_d[i] = new double[dim];
         urhsrelololo_d[i] = new double[dim];
         urhsimhihihi_d[i] = new double[dim];
         urhsimlohihi_d[i] = new double[dim];
         urhsimhilohi_d[i] = new double[dim];
         urhsimlolohi_d[i] = new double[dim];
         urhsimhihilo_d[i] = new double[dim];
         urhsimlohilo_d[i] = new double[dim];
         urhsimhilolo_d[i] = new double[dim];
         urhsimlololo_d[i] = new double[dim];
      }
   }
   double *workvecrehihihi = new double[dim];
   double *workvecrelohihi = new double[dim];
   double *workvecrehilohi = new double[dim];
   double *workvecrelolohi = new double[dim];
   double *workvecrehihilo = new double[dim];
   double *workvecrelohilo = new double[dim];
   double *workvecrehilolo = new double[dim];
   double *workvecrelololo = new double[dim];
   double *workvecimhihihi = new double[dim];
   double *workvecimlohihi = new double[dim];
   double *workvecimhilohi = new double[dim];
   double *workvecimlolohi = new double[dim];
   double *workvecimhihilo = new double[dim];
   double *workvecimlohilo = new double[dim];
   double *workvecimhilolo = new double[dim];
   double *workvecimlololo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhsrehihihi = new double*[degp1];
   double **workrhsrelohihi = new double*[degp1];
   double **workrhsrehilohi = new double*[degp1];
   double **workrhsrelolohi = new double*[degp1];
   double **workrhsrehihilo = new double*[degp1];
   double **workrhsrelohilo = new double*[degp1];
   double **workrhsrehilolo = new double*[degp1];
   double **workrhsrelololo = new double*[degp1];
   double **workrhsimhihihi = new double*[degp1];
   double **workrhsimlohihi = new double*[degp1];
   double **workrhsimhilohi = new double*[degp1];
   double **workrhsimlolohi = new double*[degp1];
   double **workrhsimhihilo = new double*[degp1];
   double **workrhsimlohilo = new double*[degp1];
   double **workrhsimhilolo = new double*[degp1];
   double **workrhsimlololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhsrehihihi[i] = new double[dim];
      workrhsrelohihi[i] = new double[dim];
      workrhsrehilohi[i] = new double[dim];
      workrhsrelolohi[i] = new double[dim];
      workrhsrehihilo[i] = new double[dim];
      workrhsrelohilo[i] = new double[dim];
      workrhsrehilolo[i] = new double[dim];
      workrhsrelololo[i] = new double[dim];
      workrhsimhihihi[i] = new double[dim];
      workrhsimlohihi[i] = new double[dim];
      workrhsimhilohi[i] = new double[dim];
      workrhsimlolohi[i] = new double[dim];
      workrhsimhihilo[i] = new double[dim];
      workrhsimlohilo[i] = new double[dim];
      workrhsimhilolo[i] = new double[dim];
      workrhsimlololo[i] = new double[dim];
   }
   double **resvecrehihihi = new double*[degp1];
   double **resvecrelohihi = new double*[degp1];
   double **resvecrehilohi = new double*[degp1];
   double **resvecrelolohi = new double*[degp1];
   double **resvecrehihilo = new double*[degp1];
   double **resvecrelohilo = new double*[degp1];
   double **resvecrehilolo = new double*[degp1];
   double **resvecrelololo = new double*[degp1];
   double **resvecimhihihi = new double*[degp1];
   double **resvecimlohihi = new double*[degp1];
   double **resvecimhilohi = new double*[degp1];
   double **resvecimlolohi = new double*[degp1];
   double **resvecimhihilo = new double*[degp1];
   double **resvecimlohilo = new double*[degp1];
   double **resvecimhilolo = new double*[degp1];
   double **resvecimlololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecrehihihi[i] = new double[dim];
      resvecrelohihi[i] = new double[dim];
      resvecrehilohi[i] = new double[dim];
      resvecrelolohi[i] = new double[dim];
      resvecrehihilo[i] = new double[dim];
      resvecrelohilo[i] = new double[dim];
      resvecrehilolo[i] = new double[dim];
      resvecrelololo[i] = new double[dim];
      resvecimhihihi[i] = new double[dim];
      resvecimlohihi[i] = new double[dim];
      resvecimhilohi[i] = new double[dim];
      resvecimlolohi[i] = new double[dim];
      resvecimhihilo[i] = new double[dim];
      resvecimlohilo[i] = new double[dim];
      resvecimhilolo[i] = new double[dim];
      resvecimlololo[i] = new double[dim];
   }
   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;

   double **Qrehihihi_h;
   double **Qrelohihi_h;
   double **Qrehilohi_h;
   double **Qrelolohi_h;
   double **Qrehihilo_h;
   double **Qrelohilo_h;
   double **Qrehilolo_h;
   double **Qrelololo_h;
   double **Qimhihihi_h;
   double **Qimlohihi_h;
   double **Qimhilohi_h;
   double **Qimlolohi_h;
   double **Qimhihilo_h;
   double **Qimlohilo_h;
   double **Qimhilolo_h;
   double **Qimlololo_h;
   double **Qrehihihi_d;
   double **Qrelohihi_d;
   double **Qrehilohi_d;
   double **Qrelolohi_d;
   double **Qrehihilo_d;
   double **Qrelohilo_d;
   double **Qrehilolo_d;
   double **Qrelololo_d;
   double **Qimhihihi_d;
   double **Qimlohihi_d;
   double **Qimhilohi_d;
   double **Qimlolohi_d;
   double **Qimhihilo_d;
   double **Qimlohilo_d;
   double **Qimhilolo_d;
   double **Qimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      Qrehihihi_h = new double*[dim];
      Qrelohihi_h = new double*[dim];
      Qrehilohi_h = new double*[dim];
      Qrelolohi_h = new double*[dim];
      Qrehihilo_h = new double*[dim];
      Qrelohilo_h = new double*[dim];
      Qrehilolo_h = new double*[dim];
      Qrelololo_h = new double*[dim];
      Qimhihihi_h = new double*[dim];
      Qimlohihi_h = new double*[dim];
      Qimhilohi_h = new double*[dim];
      Qimlolohi_h = new double*[dim];
      Qimhihilo_h = new double*[dim];
      Qimlohilo_h = new double*[dim];
      Qimhilolo_h = new double*[dim];
      Qimlololo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qrehihihi_h[i] = new double[dim];
         Qrelohihi_h[i] = new double[dim];
         Qrehilohi_h[i] = new double[dim];
         Qrelolohi_h[i] = new double[dim];
         Qrehihilo_h[i] = new double[dim];
         Qrelohilo_h[i] = new double[dim];
         Qrehilolo_h[i] = new double[dim];
         Qrelololo_h[i] = new double[dim];
         Qimhihihi_h[i] = new double[dim];
         Qimlohihi_h[i] = new double[dim];
         Qimhilohi_h[i] = new double[dim];
         Qimlolohi_h[i] = new double[dim];
         Qimhihilo_h[i] = new double[dim];
         Qimlohilo_h[i] = new double[dim];
         Qimhilolo_h[i] = new double[dim];
         Qimlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qrehihihi_d = new double*[dim];
      Qrelohihi_d = new double*[dim];
      Qrehilohi_d = new double*[dim];
      Qrelolohi_d = new double*[dim];
      Qrehihilo_d = new double*[dim];
      Qrelohilo_d = new double*[dim];
      Qrehilolo_d = new double*[dim];
      Qrelololo_d = new double*[dim];
      Qimhihihi_d = new double*[dim];
      Qimlohihi_d = new double*[dim];
      Qimhilohi_d = new double*[dim];
      Qimlolohi_d = new double*[dim];
      Qimhihilo_d = new double*[dim];
      Qimlohilo_d = new double*[dim];
      Qimhilolo_d = new double*[dim];
      Qimlololo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qrehihihi_d[i] = new double[dim];
         Qrelohihi_d[i] = new double[dim];
         Qrehilohi_d[i] = new double[dim];
         Qrelolohi_d[i] = new double[dim];
         Qrehihilo_d[i] = new double[dim];
         Qrelohilo_d[i] = new double[dim];
         Qrehilolo_d[i] = new double[dim];
         Qrelololo_d[i] = new double[dim];
         Qimhihihi_d[i] = new double[dim];
         Qimlohihi_d[i] = new double[dim];
         Qimhilohi_d[i] = new double[dim];
         Qimlolohi_d[i] = new double[dim];
         Qimhihilo_d[i] = new double[dim];
         Qimlohilo_d[i] = new double[dim];
         Qimhilolo_d[i] = new double[dim];
         Qimlololo_d[i] = new double[dim];
      }
   }
   double **Rrehihihi_h;
   double **Rrelohihi_h;
   double **Rrehilohi_h;
   double **Rrelolohi_h;
   double **Rrehihilo_h;
   double **Rrelohilo_h;
   double **Rrehilolo_h;
   double **Rrelololo_h;
   double **Rimhihihi_h;
   double **Rimlohihi_h;
   double **Rimhilohi_h;
   double **Rimlolohi_h;
   double **Rimhihilo_h;
   double **Rimlohilo_h;
   double **Rimhilolo_h;
   double **Rimlololo_h;
   double **Rrehihihi_d;
   double **Rrelohihi_d;
   double **Rrehilohi_d;
   double **Rrelolohi_d;
   double **Rrehihilo_d;
   double **Rrelohilo_d;
   double **Rrehilolo_d;
   double **Rrelololo_d;
   double **Rimhihihi_d;
   double **Rimlohihi_d;
   double **Rimhilohi_d;
   double **Rimlolohi_d;
   double **Rimhihilo_d;
   double **Rimlohilo_d;
   double **Rimhilolo_d;
   double **Rimlololo_d;

   if((mode == 1) || (mode == 2))
   {
      Rrehihihi_h = new double*[dim];
      Rrelohihi_h = new double*[dim];
      Rrehilohi_h = new double*[dim];
      Rrelolohi_h = new double*[dim];
      Rrehihilo_h = new double*[dim];
      Rrelohilo_h = new double*[dim];
      Rrehilolo_h = new double*[dim];
      Rrelololo_h = new double*[dim];
      Rimhihihi_h = new double*[dim];
      Rimlohihi_h = new double*[dim];
      Rimhilohi_h = new double*[dim];
      Rimlolohi_h = new double*[dim];
      Rimhihilo_h = new double*[dim];
      Rimlohilo_h = new double*[dim];
      Rimhilolo_h = new double*[dim];
      Rimlololo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rrehihihi_h[i] = new double[dim];
         Rrelohihi_h[i] = new double[dim];
         Rrehilohi_h[i] = new double[dim];
         Rrelolohi_h[i] = new double[dim];
         Rrehihilo_h[i] = new double[dim];
         Rrelohilo_h[i] = new double[dim];
         Rrehilolo_h[i] = new double[dim];
         Rrelololo_h[i] = new double[dim];
         Rimhihihi_h[i] = new double[dim];
         Rimlohihi_h[i] = new double[dim];
         Rimhilohi_h[i] = new double[dim];
         Rimlolohi_h[i] = new double[dim];
         Rimhihilo_h[i] = new double[dim];
         Rimlohilo_h[i] = new double[dim];
         Rimhilolo_h[i] = new double[dim];
         Rimlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rrehihihi_d = new double*[dim];
      Rrelohihi_d = new double*[dim];
      Rrehilohi_d = new double*[dim];
      Rrelolohi_d = new double*[dim];
      Rrehihilo_d = new double*[dim];
      Rrelohilo_d = new double*[dim];
      Rrehilolo_d = new double*[dim];
      Rrelololo_d = new double*[dim];
      Rimhihihi_d = new double*[dim];
      Rimlohihi_d = new double*[dim];
      Rimhilohi_d = new double*[dim];
      Rimlolohi_d = new double*[dim];
      Rimhihilo_d = new double*[dim];
      Rimlohilo_d = new double*[dim];
      Rimhilolo_d = new double*[dim];
      Rimlololo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rrehihihi_d[i] = new double[dim];
         Rrelohihi_d[i] = new double[dim];
         Rrehilohi_d[i] = new double[dim];
         Rrelolohi_d[i] = new double[dim];
         Rrehihilo_d[i] = new double[dim];
         Rrelohilo_d[i] = new double[dim];
         Rrehilolo_d[i] = new double[dim];
         Rrelololo_d[i] = new double[dim];
         Rimhihihi_d[i] = new double[dim];
         Rimlohihi_d[i] = new double[dim];
         Rimhilohi_d[i] = new double[dim];
         Rimlolohi_d[i] = new double[dim];
         Rimhihilo_d[i] = new double[dim];
         Rimlohilo_d[i] = new double[dim];
         Rimhilolo_d[i] = new double[dim];
         Rimlololo_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **startsolrehihihi = new double*[dim];
   double **startsolrelohihi = new double*[dim];
   double **startsolrehilohi = new double*[dim];
   double **startsolrelolohi = new double*[dim];
   double **startsolrehihilo = new double*[dim];
   double **startsolrelohilo = new double*[dim];
   double **startsolrehilolo = new double*[dim];
   double **startsolrelololo = new double*[dim];
   double **startsolimhihihi = new double*[dim];
   double **startsolimlohihi = new double*[dim];
   double **startsolimhilohi = new double*[dim];
   double **startsolimlolohi = new double*[dim];
   double **startsolimhihilo = new double*[dim];
   double **startsolimlohilo = new double*[dim];
   double **startsolimhilolo = new double*[dim];
   double **startsolimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      startsolrehihihi[i] = new double[degp1];
      startsolrelohihi[i] = new double[degp1];
      startsolrehilohi[i] = new double[degp1];
      startsolrelolohi[i] = new double[degp1];
      startsolrehihilo[i] = new double[degp1];
      startsolrelohilo[i] = new double[degp1];
      startsolrehilolo[i] = new double[degp1];
      startsolrelololo[i] = new double[degp1];
      startsolimhihihi[i] = new double[degp1];
      startsolimlohihi[i] = new double[degp1];
      startsolimhilohi[i] = new double[degp1];
      startsolimlolohi[i] = new double[degp1];
      startsolimhihilo[i] = new double[degp1];
      startsolimlohilo[i] = new double[degp1];
      startsolimhilolo[i] = new double[degp1];
      startsolimlololo[i] = new double[degp1];
   }
   make_complex8_exponentials
      (dim,deg,
       startsolrehihihi,startsolrelohihi,startsolrehilohi,startsolrelolohi,
       startsolrehihilo,startsolrelohilo,startsolrehilolo,startsolrelololo,
       startsolimhihihi,startsolimlohihi,startsolimhilohi,startsolimlolohi,
       startsolimhihilo,startsolimlohilo,startsolimhilolo,startsolimlololo);

   if(nbrcol != 1) // generate coefficients for the columns
      make_complex8_coefficients
         (nbrcol,dim,
          cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
          cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
          cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
          cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);

   // compute the right hand sides via evaluation

   double **mbrhsrehihihi = new double*[dim];
   double **mbrhsrelohihi = new double*[dim];
   double **mbrhsrehilohi = new double*[dim];
   double **mbrhsrelolohi = new double*[dim];
   double **mbrhsrehihilo = new double*[dim];
   double **mbrhsrelohilo = new double*[dim];
   double **mbrhsrehilolo = new double*[dim];
   double **mbrhsrelololo = new double*[dim];
   double **mbrhsimhihihi = new double*[dim];
   double **mbrhsimlohihi = new double*[dim];
   double **mbrhsimhilohi = new double*[dim];
   double **mbrhsimlolohi = new double*[dim];
   double **mbrhsimhihilo = new double*[dim];
   double **mbrhsimlohilo = new double*[dim];
   double **mbrhsimhilolo = new double*[dim];
   double **mbrhsimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhsrehihihi[i] = new double[degp1];
      mbrhsrelohihi[i] = new double[degp1];
      mbrhsrehilohi[i] = new double[degp1];
      mbrhsrelolohi[i] = new double[degp1];
      mbrhsrehihilo[i] = new double[degp1];
      mbrhsrelohilo[i] = new double[degp1];
      mbrhsrehilolo[i] = new double[degp1];
      mbrhsrelololo[i] = new double[degp1];
      mbrhsimhihihi[i] = new double[degp1];
      mbrhsimlohihi[i] = new double[degp1];
      mbrhsimhilohi[i] = new double[degp1];
      mbrhsimlolohi[i] = new double[degp1];
      mbrhsimhihilo[i] = new double[degp1];
      mbrhsimlohilo[i] = new double[degp1];
      mbrhsimhilolo[i] = new double[degp1];
      mbrhsimlololo[i] = new double[degp1];

      mbrhsrehihihi[i][0] = 1.0;     // initialize product to one
      mbrhsrelohihi[i][0] = 0.0;
      mbrhsrehilohi[i][0] = 0.0;
      mbrhsrelolohi[i][0] = 0.0;
      mbrhsrehihilo[i][0] = 0.0; 
      mbrhsrelohilo[i][0] = 0.0;
      mbrhsrehilolo[i][0] = 0.0;
      mbrhsrelololo[i][0] = 0.0;
      mbrhsimhihihi[i][0] = 0.0;
      mbrhsimlohihi[i][0] = 0.0;
      mbrhsimhilohi[i][0] = 0.0;
      mbrhsimlolohi[i][0] = 0.0;
      mbrhsimhihilo[i][0] = 0.0; 
      mbrhsimlohilo[i][0] = 0.0;
      mbrhsimhilolo[i][0] = 0.0;
      mbrhsimlololo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhsrehihihi[i][k] = 0.0; mbrhsrelohihi[i][k] = 0.0;
         mbrhsrehilohi[i][k] = 0.0; mbrhsrelolohi[i][k] = 0.0;
         mbrhsrehihilo[i][k] = 0.0; mbrhsrelohilo[i][k] = 0.0;
         mbrhsrehilolo[i][k] = 0.0; mbrhsrelololo[i][k] = 0.0;
         mbrhsimhihihi[i][k] = 0.0; mbrhsimlohihi[i][k] = 0.0;
         mbrhsimhilohi[i][k] = 0.0; mbrhsimlolohi[i][k] = 0.0;
         mbrhsimhihilo[i][k] = 0.0; mbrhsimlohilo[i][k] = 0.0;
         mbrhsimhilolo[i][k] = 0.0; mbrhsimlololo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_complex8_monomials
         (dim,deg,rowsA,
          startsolrehihihi,startsolrelohihi,startsolrehilohi,startsolrelolohi,
          startsolrehihilo,startsolrelohilo,startsolrehilolo,startsolrelololo,
          startsolimhihihi,startsolimlohihi,startsolimhilohi,startsolimlolohi,
          startsolimhihilo,startsolimlohilo,startsolimhilolo,startsolimlololo,
          mbrhsrehihihi,mbrhsrelohihi,mbrhsrehilohi,mbrhsrelolohi,
          mbrhsrehihilo,mbrhsrelohilo,mbrhsrehilolo,mbrhsrelololo,
          mbrhsimhihihi,mbrhsimlohihi,mbrhsimhilohi,mbrhsimlolohi,
          mbrhsimhihilo,mbrhsimlohilo,mbrhsimhilolo,mbrhsimlololo);
   else
      evaluate_complex8_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
          cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
          cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
          cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
          startsolrehihihi,startsolrelohihi,startsolrehilohi,startsolrelolohi,
          startsolrehihilo,startsolrelohilo,startsolrehilolo,startsolrelololo,
          startsolimhihihi,startsolimlohihi,startsolimhilohi,startsolimlolohi,
          startsolimhihilo,startsolimlohilo,startsolimhilolo,startsolimlololo,
          mbrhsrehihihi,mbrhsrelohihi,mbrhsrehilohi,mbrhsrelolohi,
          mbrhsrehihilo,mbrhsrelohilo,mbrhsrehilolo,mbrhsrelololo,
          mbrhsimhihihi,mbrhsimlohihi,mbrhsimhilohi,mbrhsimlolohi,
          mbrhsimhihilo,mbrhsimlohilo,mbrhsimhilolo,mbrhsimlololo,vrblvl);

   // rhs coefficients are c(t) = (1-t)*c(t) = c(t) - t*c(t)
   for(int i=0; i<dim; i++)
      for(int j=1; j<degp1; j++) // mbrhs[i][j] = mbrhs[i][j] - mbrhs[i][j-1];
      {
         double acchihihi,acclohihi,acchilohi,acclolohi;
         double acchihilo,acclohilo,acchilolo,acclololo;

         odf_sub(mbrhsrehihihi[i][j],  mbrhsrelohihi[i][j],
                 mbrhsrehilohi[i][j],  mbrhsrelolohi[i][j],
                 mbrhsrehihilo[i][j],  mbrhsrelohilo[i][j],
                 mbrhsrehilolo[i][j],  mbrhsrelololo[i][j],
                 mbrhsrehihihi[i][j-1],mbrhsrelohihi[i][j-1],
                 mbrhsrehilohi[i][j-1],mbrhsrelolohi[i][j-1],
                 mbrhsrehihilo[i][j-1],mbrhsrelohilo[i][j-1],
                 mbrhsrehilolo[i][j-1],mbrhsrelololo[i][j-1],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         mbrhsrehihihi[i][j] = acchihihi;
         mbrhsrelohihi[i][j] = acclohihi;
         mbrhsrehilohi[i][j] = acchilohi;
         mbrhsrelolohi[i][j] = acclolohi;
         mbrhsrehihilo[i][j] = acchihilo;
         mbrhsrelohilo[i][j] = acclohilo;
         mbrhsrehilolo[i][j] = acchilolo;
         mbrhsrelololo[i][j] = acclololo;

         odf_sub(mbrhsimhihihi[i][j],  mbrhsimlohihi[i][j],
                 mbrhsimhilohi[i][j],  mbrhsimlolohi[i][j],
                 mbrhsimhihilo[i][j],  mbrhsimlohilo[i][j],
                 mbrhsimhilolo[i][j],  mbrhsimlololo[i][j],
                 mbrhsimhihihi[i][j-1],mbrhsimlohihi[i][j-1],
                 mbrhsimhilohi[i][j-1],mbrhsimlolohi[i][j-1],
                 mbrhsimhihilo[i][j-1],mbrhsimlohilo[i][j-1],
                 mbrhsimhilolo[i][j-1],mbrhsimlololo[i][j-1],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         mbrhsimhihihi[i][j] = acchihihi;
         mbrhsimlohihi[i][j] = acclohihi;
         mbrhsimhilohi[i][j] = acchilohi;
         mbrhsimlolohi[i][j] = acclolohi;
         mbrhsimhihilo[i][j] = acchihilo;
         mbrhsimlohilo[i][j] = acclohilo;
         mbrhsimhilolo[i][j] = acchilolo;
         mbrhsimlololo[i][j] = acclololo;
      }

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhsrehihihi[i][j] << "  "
                 << mbrhsrelohihi[i][j] << endl
                 << "  "
                 << mbrhsrehilohi[i][j] << "  "
                 << mbrhsrelolohi[i][j] << endl
                 << "  "
                 << mbrhsrehihilo[i][j] << "  "
                 << mbrhsrelohilo[i][j] << endl
                 << "  "
                 << mbrhsrehilolo[i][j] << "  "
                 << mbrhsrelololo[i][j] << endl
                 << mbrhsimhihihi[i][j] << "  "
                 << mbrhsimlohihi[i][j] << endl
                 << "  "
                 << mbrhsimhilohi[i][j] << "  "
                 << mbrhsimlolohi[i][j] << endl
                 << "  "
                 << mbrhsimhihilo[i][j] << "  "
                 << mbrhsimlohilo[i][j] << endl
                 << "  "
                 << mbrhsimhilolo[i][j] << "  "
                 << mbrhsimlololo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehihihi_h[i][j] = startsolrehihihi[i][j];
         inputrelohihi_h[i][j] = startsolrelohihi[i][j];
         inputrehilohi_h[i][j] = startsolrehilohi[i][j];
         inputrelolohi_h[i][j] = startsolrelolohi[i][j];
         inputrehihilo_h[i][j] = startsolrehihilo[i][j];
         inputrelohilo_h[i][j] = startsolrelohilo[i][j];
         inputrehilolo_h[i][j] = startsolrehilolo[i][j];
         inputrelololo_h[i][j] = startsolrelololo[i][j];
         inputimhihihi_h[i][j] = startsolimhihihi[i][j];
         inputimlohihi_h[i][j] = startsolimlohihi[i][j];
         inputimhilohi_h[i][j] = startsolimhilohi[i][j];
         inputimlolohi_h[i][j] = startsolimlolohi[i][j];
         inputimhihilo_h[i][j] = startsolimhihilo[i][j];
         inputimlohilo_h[i][j] = startsolimlohilo[i][j];
         inputimhilolo_h[i][j] = startsolimhilolo[i][j];
         inputimlololo_h[i][j] = startsolimlololo[i][j];
         inputrehihihi_d[i][j] = startsolrehihihi[i][j];
         inputrelohihi_d[i][j] = startsolrelohihi[i][j];
         inputrehilohi_d[i][j] = startsolrehilohi[i][j];
         inputrelolohi_d[i][j] = startsolrelolohi[i][j];
         inputrehihilo_d[i][j] = startsolrehihilo[i][j];
         inputrelohilo_d[i][j] = startsolrelohilo[i][j];
         inputrehilolo_d[i][j] = startsolrehilolo[i][j];
         inputrelololo_d[i][j] = startsolrelololo[i][j];
         inputimhihihi_d[i][j] = startsolimhihihi[i][j];
         inputimlohihi_d[i][j] = startsolimlohihi[i][j];
         inputimhilohi_d[i][j] = startsolimhilohi[i][j];
         inputimlolohi_d[i][j] = startsolimlolohi[i][j];
         inputimhihilo_d[i][j] = startsolimhihilo[i][j];
         inputimlohilo_d[i][j] = startsolimlohilo[i][j];
         inputimhilolo_d[i][j] = startsolimhilolo[i][j];
         inputimlololo_d[i][j] = startsolimlololo[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputrehihihi_h[i][0] << "  "
                            << inputrelohihi_h[i][0] << endl;
         cout << "     " << inputrehilohi_h[i][0] << "  "
                         << inputrelolohi_h[i][0] << endl;
         cout << "     " << inputrehihilo_h[i][0] << "  "
                         << inputrelohilo_h[i][0] << endl;
         cout << "     " << inputrehilolo_h[i][0] << "  "
                         << inputrelololo_h[i][0] << endl;
         cout << "     " << inputimhihihi_h[i][0] << "  "
                         << inputimlohihi_h[i][0] << endl;
         cout << "     " << inputimhilohi_h[i][0] << "  "
                         << inputimlolohi_h[i][0] << endl;
         cout << "     " << inputimhihilo_h[i][0] << "  "
                         << inputimlohilo_h[i][0] << endl;
         cout << "     " << inputimhilolo_h[i][0] << "  "
                         << inputimlololo_h[i][0] << endl;
      }
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   cmplx8_run_newton
      (szt,nbt,dim,deg,nbrcol,nbsteps,nvr,idx,exp,nbrfac,expfac,
       mbrhsrehihihi,mbrhsrelohihi,mbrhsrehilohi,mbrhsrelolohi,
       mbrhsrehihilo,mbrhsrelohilo,mbrhsrehilolo,mbrhsrelololo,
       mbrhsimhihihi,mbrhsimlohihi,mbrhsimhilohi,mbrhsimlolohi,
       mbrhsimhihilo,mbrhsimlohilo,mbrhsimhilolo,mbrhsimlololo,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       accrehihihi,accrelohihi,accrehilohi,accrelolohi,
       accrehihilo,accrelohilo,accrehilolo,accrelololo,
       accimhihihi,accimlohihi,accimhilohi,accimlolohi,
       accimhihilo,accimlohilo,accimhilolo,accimlololo,
       inputrehihihi_h,inputrelohihi_h,inputrehilohi_h,inputrelolohi_h,
       inputrehihilo_h,inputrelohilo_h,inputrehilolo_h,inputrelololo_h,
       inputimhihihi_h,inputimlohihi_h,inputimhilohi_h,inputimlolohi_h,
       inputimhihilo_h,inputimlohilo_h,inputimhilolo_h,inputimlololo_h,
       inputrehihihi_d,inputrelohihi_d,inputrehilohi_d,inputrelolohi_d,
       inputrehihilo_d,inputrelohilo_d,inputrehilolo_d,inputrelololo_d,
       inputimhihihi_d,inputimlohihi_d,inputimhilohi_d,inputimlolohi_d,
       inputimhihilo_d,inputimlohilo_d,inputimhilolo_d,inputimlololo_d,
       outputrehihihi_h,outputrelohihi_h,outputrehilohi_h,outputrelolohi_h,
       outputrehihilo_h,outputrelohilo_h,outputrehilolo_h,outputrelololo_h,
       outputimhihihi_h,outputimlohihi_h,outputimhilohi_h,outputimlolohi_h,
       outputimhihilo_h,outputimlohilo_h,outputimhilolo_h,outputimlololo_h,
       outputrehihihi_d,outputrelohihi_d,outputrehilohi_d,outputrelolohi_d,
       outputrehihilo_d,outputrelohilo_d,outputrehilolo_d,outputrelololo_d,
       outputimhihihi_d,outputimlohihi_d,outputimhilohi_d,outputimlolohi_d,
       outputimhihilo_d,outputimlohilo_d,outputimhilolo_d,outputimlololo_d,
       funvalrehihihi_h,funvalrelohihi_h,funvalrehilohi_h,funvalrelolohi_h,
       funvalrehihilo_h,funvalrelohilo_h,funvalrehilolo_h,funvalrelololo_h,
       funvalimhihihi_h,funvalimlohihi_h,funvalimhilohi_h,funvalimlolohi_h,
       funvalimhihilo_h,funvalimlohilo_h,funvalimhilolo_h,funvalimlololo_h,
       funvalrehihihi_d,funvalrelohihi_d,funvalrehilohi_d,funvalrelolohi_d,
       funvalrehihilo_d,funvalrelohilo_d,funvalrehilolo_d,funvalrelololo_d,
       funvalimhihihi_d,funvalimlohihi_d,funvalimhilohi_d,funvalimlolohi_d,
       funvalimhihilo_d,funvalimlohilo_d,funvalimhilolo_d,funvalimlololo_d,
       jacvalrehihihi_h,jacvalrelohihi_h,jacvalrehilohi_h,jacvalrelolohi_h,
       jacvalrehihilo_h,jacvalrelohilo_h,jacvalrehilolo_h,jacvalrelololo_h,
       jacvalimhihihi_h,jacvalimlohihi_h,jacvalimhilohi_h,jacvalimlolohi_h,
       jacvalimhihilo_h,jacvalimlohilo_h,jacvalimhilolo_h,jacvalimlololo_h,
       jacvalrehihihi_d,jacvalrelohihi_d,jacvalrehilohi_d,jacvalrelolohi_d,
       jacvalrehihilo_d,jacvalrelohilo_d,jacvalrehilolo_d,jacvalrelololo_d,
       jacvalimhihihi_d,jacvalimlohihi_d,jacvalimhilohi_d,jacvalimlolohi_d,
       jacvalimhihilo_d,jacvalimlohilo_d,jacvalimhilolo_d,jacvalimlololo_d,
       rhsrehihihi_h,rhsrelohihi_h,rhsrehilohi_h,rhsrelolohi_h,
       rhsrehihilo_h,rhsrelohilo_h,rhsrehilolo_h,rhsrelololo_h,
       rhsimhihihi_h,rhsimlohihi_h,rhsimhilohi_h,rhsimlolohi_h,
       rhsimhihilo_h,rhsimlohilo_h,rhsimhilolo_h,rhsimlololo_h,
       rhsrehihihi_d,rhsrelohihi_d,rhsrehilohi_d,rhsrelolohi_d,
       rhsrehihilo_d,rhsrelohilo_d,rhsrehilolo_d,rhsrelololo_d,
       rhsimhihihi_d,rhsimlohihi_d,rhsimhilohi_d,rhsimlolohi_d,
       rhsimhihilo_d,rhsimlohilo_d,rhsimhilolo_d,rhsimlololo_d,
       urhsrehihihi_h,urhsrelohihi_h,urhsrehilohi_h,urhsrelolohi_h,
       urhsrehihilo_h,urhsrelohilo_h,urhsrehilolo_h,urhsrelololo_h,
       urhsimhihihi_h,urhsimlohihi_h,urhsimhilohi_h,urhsimlolohi_h,
       urhsimhihilo_h,urhsimlohilo_h,urhsimhilolo_h,urhsimlololo_h,
       urhsrehihihi_d,urhsrelohihi_d,urhsrehilohi_d,urhsrelolohi_d,
       urhsrehihilo_d,urhsrelohilo_d,urhsrehilolo_d,urhsrelololo_d,
       urhsimhihihi_d,urhsimlohihi_d,urhsimhilohi_d,urhsimlolohi_d,
       urhsimhihilo_d,urhsimlohilo_d,urhsimhilolo_d,urhsimlololo_d,
       solrehihihi_h,solrelohihi_h,solrehilohi_h,solrelolohi_h,
       solrehihilo_h,solrelohilo_h,solrehilolo_h,solrelololo_h,
       solimhihihi_h,solimlohihi_h,solimhilohi_h,solimlolohi_h,
       solimhihilo_h,solimlohilo_h,solimhilolo_h,solimlololo_h,
       solrehihihi_d,solrelohihi_d,solrehilohi_d,solrelolohi_d,
       solrehihilo_d,solrelohilo_d,solrehilolo_d,solrelololo_d,
       solimhihihi_d,solimlohihi_d,solimhilohi_d,solimlolohi_d,
       solimhihilo_d,solimlohilo_d,solimhilolo_d,solimlololo_d,
       Qrehihihi_h,Qrelohihi_h,Qrehilohi_h,Qrelolohi_h,
       Qrehihilo_h,Qrelohilo_h,Qrehilolo_h,Qrelololo_h,
       Qimhihihi_h,Qimlohihi_h,Qimhilohi_h,Qimlolohi_h,
       Qimhihilo_h,Qimlohilo_h,Qimhilolo_h,Qimlololo_h,
       Qrehihihi_d,Qrelohihi_d,Qrehilohi_d,Qrelolohi_d,
       Qrehihilo_d,Qrelohilo_d,Qrehilolo_d,Qrelololo_d,
       Qimhihihi_d,Qimlohihi_d,Qimhilohi_d,Qimlolohi_d,
       Qimhihilo_d,Qimlohilo_d,Qimhilolo_d,Qimlololo_d,
       Rrehihihi_h,Rrelohihi_h,Rrehilohi_h,Rrelolohi_h,
       Rrehihilo_h,Rrelohilo_h,Rrehilolo_h,Rrelololo_h,
       Rimhihihi_h,Rimlohihi_h,Rimhilohi_h,Rimlolohi_h,
       Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,
       Rrehihihi_d,Rrelohihi_d,Rrehilohi_d,Rrelolohi_d,
       Rrehihilo_d,Rrelohilo_d,Rrehilolo_d,Rrelololo_d,
       Rimhihihi_d,Rimlohihi_d,Rimhilohi_d,Rimlolohi_d,
       Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,
       workvecrehihihi,workvecrelohihi,workvecrehilohi,workvecrelolohi,
       workvecrehihilo,workvecrelohilo,workvecrehilolo,workvecrelololo,
       workvecimhihihi,workvecimlohihi,workvecimhilohi,workvecimlolohi,
       workvecimhihilo,workvecimlohilo,workvecimhilolo,workvecimlololo,
       resvecrehihihi,resvecrelohihi,resvecrehilohi,resvecrelolohi,
       resvecrehihilo,resvecrelohilo,resvecrehilolo,resvecrelololo,
       resvecimhihihi,resvecimlohihi,resvecimhilohi,resvecimlolohi,
       resvecimhihilo,resvecimlohilo,resvecimhilolo,resvecimlololo,
       &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
       &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
       vrblvl,mode);

   return 0;
}
