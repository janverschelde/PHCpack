/* The file dbl5_convolutions_host.h specifies functions for the product of 
 * two series, truncated to the same degree, in penta double precision. */

#ifndef __dbl5_convolutions_host_h__
#define __dbl5_convolutions_host_h__

void CPU_dbl5_product
 ( int deg, double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
            double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
            double *ztb, double *zix, double *zmi, double *zrg, double *zpk );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for real coefficients in penta double precision.
 *
 * REQUIRED :
 *   The arrays xtb, xix, xmi, xrg, xpk, ytb, yix, ymi, yrg, ypk, ztb, zix,
 *   zmi, zrg, and zpk have space allocated for deg+1 doubles, 
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xtb      highest parts of the coefficients of x;
 *   xix      second highest parts of the coefficients of x;
 *   xmi      middle parts of the coefficients of x;
 *   xrg      second lowest parts of the coefficients of x;
 *   xpk      lowest parts of the coefficients of x;
 *   ytb      highest parts of the coefficients of y;
 *   yix      second highest parts of the coefficients of y;
 *   ymi      middle parts of the coefficients of y;
 *   yrg      second lowest parts of the coefficients of y;
 *   ypk      lowest parts of the coefficients of y.
 *
 * ON RETURN :
 *   ztb      highest parts of the coefficients of the product x*y;
 *   zix      second highest parts of the coefficients of the product x*y;
 *   zmi      middle parts of the coefficients of the product x*y;
 *   zrg      second lowest parts of the coefficients of the product x*y;
 *   zpk      lowest parts of the coefficients of the product x*y. */

void CPU_cmplx5_product
 ( int deg,
   double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *yretb, double *yreix, double *yremi, double *yrerg, double *yrepk,
   double *yimtb, double *yimix, double *yimmi, double *yimrg, double *yimpk,
   double *zretb, double *zreix, double *zremi, double *zrerg, double *zrepk,
   double *zimtb, double *zimix, double *zimmi, double *zimrg, double *zimpk );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for complex coefficients in penta double precision.
 *
 * REQUIRED :
 *   All arrays xretb, xreix, xremi, xrerg, xrepk, ximtb, ximix, ximmi, ximrg,
 *   ximpk, yretb, yreix, yremi, yrerg, yrepk, yimtb, yimix, yimmi, yimrg,
 *   yimpk, zretb, zreix, zremi, zrerg, zrepk, zimtb, zimix, zimmi, zimrg,
 *   and zimpk have space allocated for deg+1 doubles, for range 0 up to deg,
 *   deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xretb    highest real parts of the coefficients of x;
 *   xreix    second highest real parts of the coefficients of x;
 *   xremi    middle real parts of the coefficients of x;
 *   xrerg    second lowest real parts of the coefficients of x;
 *   xrepk    lowest real parts of the coefficients of x;
 *   ximtb    highest imaginary parts of the coefficients of x;
 *   ximix    second highest imaginary parts of the coefficients of x;
 *   ximmi    middle imaginary parts of the coefficients of x;
 *   ximrg    second lowest imaginary parts of the coefficients of x;
 *   ximpk    lowest imaginary parts of the coefficients of x;
 *   yretb    highest real parts of the coefficients of y;
 *   yreix    second highest real parts of the coefficients of y;
 *   yremi    middle real parts of the coefficients of y;
 *   yrerg    second lowest real parts of the coefficients of y;
 *   yrepk    lowest real parts of the coefficients of y;
 *   yimtb    highest imaginary parts of the coefficients of y;
 *   yimix    second highest imaginary parts of the coefficients of y;
 *   yimmi    middle imaginary parts of the coefficients of y;
 *   yimrg    second lowest imaginary parts of the coefficients of y;
 *   yimpk    lowest imaginary parts of the coefficients of y.
 *
 * ON RETURN :
 *   zretb    highest real parts of the coefficients of x*y;
 *   zreix    second highest real parts of the coefficients of x*y;
 *   zremi    middle real parts of the coefficients of x*y;
 *   zrerg    second lowest real parts of the coefficients of x*y;
 *   zrepk    lowest real parts of the coefficients of x*y;
 *   zimtb    highest imaginary parts of the coefficients of x*y;
 *   zimix    second highest imaginary parts of the coefficients of x*y;
 *   zimmi    middle imaginary parts of the coefficients of x*y;
 *   zimrg    second lowest imaginary parts of the coefficients of x*y;
 *   zimpk    lowest imaginary parts of the coefficients of x*y. */

#endif
