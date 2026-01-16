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

void quarter_hd_vector
 ( int dim,
   double *xhihihihi, double *xlohihihi,
   double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi,
   double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo,
   double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo,
   double *xhilololo, double *xlolololo,
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
 * Given a vector of size dim in 16 arrays,
 * quarters the parts of the hexa doubles,
 * resulting in 64 arrays of size dim. */

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

void hd_write_vector
 ( int dim,
   double *xhihihihi, double *xlohihihi, double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi, double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo, double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo, double *xhilololo, double *xlolololo );
/*
 * Writes the hexa double vector of size dim with parts in
 * xhihihihi, xlohihihi, xhilohihi, xlolohihi, xhihilohi, xlohilohi,
 * xhilolohi, xlololohi, xhihihilo, xlohihilo, xhilohilo, xlolohilo,
 * xhihilolo, xlohilolo, xhilololo, xlolololo. */

void hexa_double_product
 ( int dim,
   double *xhihihihi, double *xlohihihi, double *xhilohihi, double *xlolohihi,
   double *xhihilohi, double *xlohilohi, double *xhilolohi, double *xlololohi,
   double *xhihihilo, double *xlohihilo, double *xhilohilo, double *xlolohilo,
   double *xhihilolo, double *xlohilolo, double *xhilololo, double *xlolololo,
   double *yhihihihi, double *ylohihihi, double *yhilohihi, double *ylolohihi,
   double *yhihilohi, double *ylohilohi, double *yhilolohi, double *ylololohi,
   double *yhihihilo, double *ylohihilo, double *yhilohilo, double *ylolohilo,
   double *yhihilolo, double *ylohilolo, double *yhilololo, double *ylolololo,
   double *phihihihi, double *plohihihi, double *philohihi, double *plolohihi,
   double *phihilohi, double *plohilohi, double *philolohi, double *plololohi,
   double *phihihilo, double *plohihilo, double *philohilo, double *plolohilo,
   double *phihilolo, double *plohilolo, double *philololo, double *plolololo );
/*
 * Makes the product of two hexa double vectors of size dim, given in
 * (xhihihihi, xlohihihi, xhilohihi, xlolohihi, xhihilohi, xlohilohi,
 *  xhilolohi, xlololohi, xhihihilo, xlohihilo, xhilohilo, xlolohilo,
 *  xhihilolo, xlohilolo, xhilololo, xlolololo) 
 * and
 * (yhihihihi, ylohihihi, yhilohihi, ylolohihi, yhihilohi, ylohilohi,
 *  yhilolohi, ylololohi, yhihihilo, ylohihilo, yhilohilo, ylolohilo,
 *  yhihilolo, ylohilolo, yhilololo, ylolololo),
 * resulting in
 * (phihihihi, plohihihi, philohihi, plolohihi, phihilohi, plohilohi,
 *  philolohi, plololohi, phihihilo, plohihilo, philohilo, plolohilo,
 *  phihilolo, plohilolo, philololo, plolololo). */

void vectored_hd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *x8, double *x9, double *x10, double *x11,
   double *x12, double *x13, double *x14, double *x15,
   double *x16, double *x17, double *x18, double *x19,
   double *x20, double *x21, double *x22, double *x23,
   double *x24, double *x25, double *x26, double *x27,
   double *x28, double *x29, double *x30, double *x31,
   double *x32, double *x33, double *x34, double *x35,
   double *x36, double *x37, double *x38, double *x39,
   double *x40, double *x41, double *x42, double *x43,
   double *x44, double *x45, double *x46, double *x47,
   double *x48, double *x49, double *x50, double *x51,
   double *x52, double *x53, double *x54, double *x55,
   double *x56, double *x57, double *x58, double *x59,
   double *x60, double *x61, double *x62, double *x63,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *y8, double *y9, double *y10, double *y11,
   double *y12, double *y13, double *y14, double *y15,
   double *y16, double *y17, double *y18, double *y19,
   double *y20, double *y21, double *y22, double *y23,
   double *y24, double *y25, double *y26, double *y27,
   double *y28, double *y29, double *y30, double *y31,
   double *y32, double *y33, double *y34, double *y35,
   double *y36, double *y37, double *y38, double *y39,
   double *y40, double *y41, double *y42, double *y43,
   double *y44, double *y45, double *y46, double *y47,
   double *y48, double *y49, double *y50, double *y51,
   double *y52, double *y53, double *y54, double *y55,
   double *y56, double *y57, double *y58, double *y59,
   double *y60, double *y61, double *y62, double *y63,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7,
   double *s8, double *s9, double *s10, double *s11,
   double *s12, double *s13, double *s14, double *s15,
   double *s16, double *s17, double *s18, double *s19,
   double *s20, double *s21, double *s22, double *s23,
   double *s24, double *s25, double *s26, double *s27,
   double *s28, double *s29, double *s30, double *s31,
   double *s32, double *s33, double *s34, double *s35,
   double *s36, double *s37, double *s38, double *s39,
   double *s40, double *s41, double *s42, double *s43,
   double *s44, double *s45, double *s46, double *s47,
   double *s48, double *s49, double *s50, double *s51,
   double *s52, double *s53, double *s54, double *s55,
   double *s56, double *s57, double *s58, double *s59,
   double *s60, double *s61, double *s62, double *s63 );
/*
 * Makes the vectored product of x and y, with the sums of the product
 * in s0, s1, etc ... */

#endif
