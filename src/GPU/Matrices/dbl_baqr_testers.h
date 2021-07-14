// The file dbl_baqr_testers.h specifies test functions on
// blocked accelerated QR decomposition in double precision.

#ifndef __dbl_baqr_testers_h__
#define __dbl_baqr_testers_h__

void test_real_blocked_qr ( void );
/*
 * DESCRIPTION :
 *   Generates a random real upper triangular matrix
 *   to test the computation of the QR decomposition. */

void test_cmplx_blocked_qr ( void );
/*
 * DESCRIPTION :
 *   Generates a random complex upper triangular matrix
 *   to test the computation of the QR decomposition. */

#endif
