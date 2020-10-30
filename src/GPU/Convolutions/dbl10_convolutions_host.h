/* The file dbl10_convolutions_host.h specifies functions for the product of 
 * two series, truncated to the same degree, in deca double precision. */

#ifndef __dbl10_convolutions_host_h__
#define __dbl10_convolutions_host_h__

void CPU_dbl10_product
 ( int deg,
   double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for real coefficients in deca double precision.
 *
 * REQUIRED :
 *   All arrays have space allocated for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xrtb     highest parts of the coefficients of x;
 *   xrix     second highest parts of the coefficients of x;
 *   xrmi     third highest parts of the coefficients of x;
 *   xrrg     fourth highest parts of the coefficients of x;
 *   xrpk     fifth highest parts of the coefficients of x;
 *   xltb     fifth lowest parts of the coefficients of x;
 *   xlix     fourth lowest parts of the coefficients of x;
 *   xlmi     third lowest parts of the coefficients of x;
 *   xlrg     second lowest parts of the coefficients of x;
 *   xlpk     lowest parts of the coefficients of x;
 *   yrtb     highest parts of the coefficients of y;
 *   yrix     second highest parts of the coefficients of y;
 *   yrmi     middle parts of the coefficients of y;
 *   yrrg     second lowest parts of the coefficients of y;
 *   yrpk     lowest parts of the coefficients of y;
 *   yltb     highest parts of the coefficients of y;
 *   ylix     second highest parts of the coefficients of y;
 *   ylmi     middle parts of the coefficients of y;
 *   ylrg     second lowest parts of the coefficients of y;
 *   ylpk     lowest parts of the coefficients of y.
 *
 * ON RETURN :
 *   zrtb     highest parts of the coefficients of the product x*y;
 *   zrix     second highest parts of the coefficients of the product x*y;
 *   zrmi     middle parts of the coefficients of the product x*y;
 *   zrrg     second lowest parts of the coefficients of the product x*y;
 *   zrpk     lowest parts of the coefficients of the product x*y;
 *   zltb     highest parts of the coefficients of the product x*y;
 *   zlix     second highest parts of the coefficients of the product x*y;
 *   zlmi     middle parts of the coefficients of the product x*y;
 *   zlrg     second lowest parts of the coefficients of the product x*y;
 *   zlpk     lowest parts of the coefficients of the product x*y. */

void CPU_cmplx10_product
 ( int deg,
   double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk, double *ximrtb, double *ximrix,
   double *ximrmi, double *ximrrg, double *ximrpk, double *ximltb,
   double *ximlix, double *ximlmi, double *ximlrg, double *ximlpk,
   double *yrertb, double *yrerix, double *yrermi, double *yrerrg,
   double *yrerpk, double *yreltb, double *yrelix, double *yrelmi,
   double *yrelrg, double *yrelpk, double *yimrtb, double *yimrix,
   double *yimrmi, double *yimrrg, double *yimrpk, double *yimltb,
   double *yimlix, double *yimlmi, double *yimlrg, double *yimlpk,
   double *zrertb, double *zrerix, double *zrermi, double *zrerrg,
   double *zrerpk, double *zreltb, double *zrelix, double *zrelmi,
   double *zrelrg, double *zrelpk, double *zimrtb, double *zimrix,
   double *zimrmi, double *zimrrg, double *zimrpk, double *zimltb,
   double *zimlix, double *zimlmi, double *zimlrg, double *zimlpk );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for complex coefficients in deca double precision.
 *
 * REQUIRED :
 *   All arrays have space allocated for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xrertb   highest real parts of the coefficients of x;
 *   xrerix   second highest real parts of the coefficients of x;
 *   xrermi   third highest real parts of the coefficients of x;
 *   xrerrg   fourth highest real parts of the coefficients of x;
 *   xrerpk   fifth highest real parts of the coefficients of x;
 *   xreltb   fifth lowest real parts of the coefficients of x;
 *   xrelix   fourth lowest real parts of the coefficients of x;
 *   xrelmi   third lowest real parts of the coefficients of x;
 *   xrelrg   second lowest real parts of the coefficients of x;
 *   xrelpk   lowest real parts of the coefficients of x;
 *   ximrtb   highest imaginary parts of the coefficients of x;
 *   ximrix   second highest imaginary parts of the coefficients of x;
 *   ximrmi   third highest imaginary parts of the coefficients of x;
 *   ximrrg   fourth highest imaginary parts of the coefficients of x;
 *   ximrpk   fifth highest imaginary parts of the coefficients of x;
 *   ximltb   fifth lowest imaginary parts of the coefficients of x;
 *   ximlix   fourth lowest imaginary parts of the coefficients of x;
 *   ximlmi   third lowest imaginary parts of the coefficients of x;
 *   ximlrg   second lowest imaginary parts of the coefficients of x;
 *   ximlpk   lowest imaginary parts of the coefficients of x;
 *   yrertb   highest real parts of the coefficients of y;
 *   yrerix   second highest real parts of the coefficients of y;
 *   yrermi   third highest real parts of the coefficients of y;
 *   yrerrg   fourth highest real parts of the coefficients of y;
 *   yrerpk   fifth highest real parts of the coefficients of y;
 *   yreltb   fifth lowest real parts of the coefficients of y;
 *   yrelix   fourth lowest real parts of the coefficients of y;
 *   yrelmi   third lowest real parts of the coefficients of y;
 *   yrelrg   second lowest real parts of the coefficients of y;
 *   yrelpk   lowest real parts of the coefficients of y;
 *   yimrtb   highest imaginary parts of the coefficients of y;
 *   yimrix   second highest imaginary parts of the coefficients of y;
 *   yimrmi   third highest imaginary parts of the coefficients of y;
 *   yimrrg   fourth highest imaginary parts of the coefficients of y;
 *   yimrpk   fifth highest imaginary parts of the coefficients of y;
 *   yimltb   fifth lowest imaginary parts of the coefficients of y;
 *   yimlix   fourth lowest imaginary parts of the coefficients of y;
 *   yimlmi   third lowest imaginary parts of the coefficients of y;
 *   yimlrg   second lowest imaginary parts of the coefficients of y;
 *   yimlpk   lowest imaginary parts of the coefficients of y.
 *
 * ON RETURN :
 *   zrertb   highest real parts of the coefficients of x*y;
 *   zrerix   second highest real parts of the coefficients of x*y;
 *   zrermi   third highest real parts of the coefficients of x*y;
 *   zrerrg   fourth highest real parts of the coefficients of x*y;
 *   zrerpk   fifth highest real parts of the coefficients of x*y;
 *   zreltb   fifth lowest real parts of the coefficients of x*y;
 *   zrelix   fourth lowest real parts of the coefficients of x*y;
 *   zrelmi   third lowest real parts of the coefficients of x*y;
 *   zrelrg   second lowest real parts of the coefficients of x*y;
 *   zrelpk   lowest real parts of the coefficients of x*y;
 *   zimrtb   highest imaginary parts of the coefficients of x*y;
 *   zimrix   second highest imaginary parts of the coefficients of x*y;
 *   zimrmi   third highest imaginary parts of the coefficients of x*y;
 *   zimrrg   fourth highest imaginary parts of the coefficients of x*y;
 *   zimrpk   fifth highest imaginary parts of the coefficients of x*y;
 *   zimltb   fifth lowest imaginary parts of the coefficients of x*y;
 *   zimlix   fourth lowest imaginary parts of the coefficients of x*y;
 *   zimlmi   third lowest imaginary parts of the coefficients of x*y;
 *   zimlrg   second lowest imaginary parts of the coefficients of x*y;
 *   zimlpk   lowest imaginary parts of the coefficients of x*y. */

#endif
