/* Collection of functions for vectored hexa double arithmetic. */

#ifndef __VECTORED_HEXA_DOUBLES_H__
#define __VECTORED_HEXA_DOUBLES_H__

void quarter_hexa_double
 ( double xhihihihi, double xlohihihi, double xhilohihi, double xlolohihi,
   double xhihilohi, double xlohilohi, double xhilolohi, double xlololohi,
   double xhihihilo, double xlohihilo, double xhilohilo, double xlolohilo,
   double xhihilolo, double xlohilolo, double xhilololo, double xlolololo,
   double *xhihihihi0, double *xhihihihi1,
   double *xhihihihi2, double *xhihihihi3,
   double *xlohihihi0, double *xlohihihi1,
   double *xlohihihi2, double *xlohihihi3,
   double *xhilohihi0, double *xhilohihi1,
   double *xhilohihi2, double *xhilohihi3,
   double *xlolohihi0, double *xlolohihi1,
   double *xlolohihi2, double *xlolohihi3,
   double *xhihilohi0, double *xhihilohi1,
   double *xhihilohi2, double *xhihilohi3,
   double *xlohilohi0, double *xlohilohi1,
   double *xlohilohi2, double *xlohilohi3,
   double *xhilolohi0, double *xhilolohi1,
   double *xhilolohi2, double *xhilolohi3,
   double *xlololohi0, double *xlololohi1,
   double *xlololohi2, double *xlololohi3,
   double *xhihihilo0, double *xhihihilo1,
   double *xhihihilo2, double *xhihihilo3,
   double *xlohihilo0, double *xlohihilo1,
   double *xlohihilo2, double *xlohihilo3,
   double *xhilohilo0, double *xhilohilo1,
   double *xhilohilo2, double *xhilohilo3,
   double *xlolohilo0, double *xlolohilo1,
   double *xlolohilo2, double *xlolohilo3,
   double *xhihilolo0, double *xhihilolo1,
   double *xhihilolo2, double *xhihilolo3,
   double *xlohilolo0, double *xlohilolo1,
   double *xlohilolo2, double *xlohilolo3,
   double *xhilololo0, double *xhilololo1,
   double *xhilololo2, double *xhilololo3,
   double *xlolololo0, double *xlolololo1,
   double *xlolololo2, double *xlolololo3 );
/*
 * Quarters the parts of a hexa double
 * (xhihihihi, xlohihihi, xhilohihi, xlolohihi, xhihilohi, xlohilohi,
 *  xhilolohi, xlololohi, xhihihilo, xlohihilo, xhilohilo, xlolohilo,
 *  xhihilolo, xlohilolo, xhilololo, xlolololo), resulting in 
 * (xhihihi0hi, xhihihihi1, xhihihihi2, xhihihihi3),
 * (xlohihi0hi, xlohihihi1, xlohihihi2, xlohihihi3),
 * (xhilohi0hi, xhilohihi1, xhilohihi2, xhilohihi3),
 * (xlolohi0hi, xlolohihi1, xlolohihi2, xlolohihi3),
 * (xhihilo0hi, xhihilohi1, xhihilohi2, xhihilohi3),
 * (xlohilo0hi, xlohilohi1, xlohilohi2, xlohilohi3),
 * (xhilolo0hi, xhilolohi1, xhilolohi2, xhilolohi3),
 * (xlololo0hi, xlololohi1, xlololohi2, xlololohi3),
 * (xhihihi0lo, xhihihi1lo, xhihihi2lo, xhihihi3lo),
 * (xlohihi0lo, xlohihi1lo, xlohihi2lo, xlohihi3lo),
 * (xhilohi0lo, xhilohi1lo, xhilohi2lo, xhilohi3lo),
 * (xlolohi0lo, xlolohi1lo, xlolohi2lo, xlolohi3lo),
 * (xhihilo0lo, xhihilo1lo, xhihilo2lo, xhihilo3lo),
 * (xlohilo0lo, xlohilo1lo, xlohilo2lo, xlohilo3lo),
 * (xhilolo0lo, xhilolo1lo, xhilolo2lo, xhilolo3lo), and
 * (xlololo0lo, xlololo1lo, xlololo2lo, xlololo3lo). */

void to_hexa_double
 ( double xhihihihi0, double xhihihihi1, double xhihihihi2, double xhihihihi3,
   double xlohihihi0, double xlohihihi1, double xlohihihi2, double xlohihihi3,
   double xhilohihi0, double xhilohihi1, double xhilohihi2, double xhilohihi3,
   double xlolohihi0, double xlolohihi1, double xlolohihi2, double xlolohihi3,
   double xhihilohi0, double xhihilohi1, double xhihilohi2, double xhihilohi3,
   double xlohilohi0, double xlohilohi1, double xlohilohi2, double xlohilohi3,
   double xhilolohi0, double xhilolohi1, double xhilolohi2, double xhilolohi3,
   double xlololohi0, double xlololohi1, double xlololohi2, double xlololohi3,
   double xhihihilo0, double xhihihilo1, double xhihihilo2, double xhihihilo3,
   double xlohihilo0, double xlohihilo1, double xlohihilo2, double xlohihilo3,
   double xhilohilo0, double xhilohilo1, double xhilohilo2, double xhilohilo3,
   double xlolohilo0, double xlolohilo1, double xlolohilo2, double xlolohilo3,
   double xhihilolo0, double xhihilolo1, double xhihilolo2, double xhihilolo3,
   double xlohilolo0, double xlohilolo1, double xlohilolo2, double xlohilolo3,
   double xhilololo0, double xhilololo1, double xhilololo2, double xhilololo3,
   double xlolololo0, double xlolololo1, double xlolololo2, double xlolololo3,
   double *xhihihihi, double *xlohihihi, double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi, double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo, double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo, double *xhilololo, double *xlolololo );
/*
 * Given the quarters of the parts of a hexa double, returns in
 * (xhihihihi, xlohihihi, xhilohihi, xlolohihi, xhihilohi, xlohilohi,
 *  xhilolohi, xlololohi, xhihihilo, xlohihilo, xhilohilo, xlolohilo,
 *  xhihilolo, xlohilolo, xhilololo, xlolololo) 
 * the parts of a hexa double, using hexa double arithmetic. */

#endif
