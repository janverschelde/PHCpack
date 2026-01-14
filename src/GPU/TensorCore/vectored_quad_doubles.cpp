/* Collection of functions for vectored quad double arithmetic. */

#include <stdio.h>
#include "quad_double_functions.h"
#include "splitting_doubles.h"

void quarter_quad_double
 ( double xhihi, double xlohi, double xhilo, double xlolo,
   double *xhihi0, double *xhihi1, double *xhihi2, double *xhihi3,
   double *xlohi0, double *xlohi1, double *xlohi2, double *xlohi3,
   double *xhilo0, double *xhilo1, double *xhilo2, double *xhilo3,
   double *xlolo0, double *xlolo1, double *xlolo2, double *xlolo3 )
{
   quarter_split(xhihi, xhihi0, xhihi1, xhihi2, xhihi3);
   quarter_split(xlohi, xlohi0, xlohi1, xlohi2, xlohi3);
   quarter_split(xhilo, xhilo0, xhilo1, xhilo2, xhilo3);
   quarter_split(xlolo, xlolo0, xlolo1, xlolo2, xlolo3);
}

void to_quad_double
 ( double xhihi0, double xhihi1, double xhihi2, double xhihi3,
   double xlohi0, double xlohi1, double xlohi2, double xlohi3,
   double xhilo0, double xhilo1, double xhilo2, double xhilo3,
   double xlolo0, double xlolo1, double xlolo2, double xlolo3,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo )
{
   *xhihi = xhihi0;
   *xlohi = 0.0;
   *xhilo = 0.0;
   *xlolo = 0.0;

   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhihi1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhihi2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhihi3);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi0);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlohi3);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo0);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xhilo3);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo0);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo1);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo2);
   qdf_inc_d(xhihi, xlohi, xhilo, xlolo, xlolo3);
}
