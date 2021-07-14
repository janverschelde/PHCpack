// The file dbl2_baqr_testers.h specifies test functions on
// blocked accelerated QR decomposition in double double precision.

#ifndef __dbl2_baqr_testers_h__
#define __dbl2_baqr_testers_h__

void test_real2_blocked_qr ( void );
/*
 * DESCRIPTION :
 *   Generates a random real matrix to test the blocked
 *   computation of the QR decomposition. */

void test_cmplx2_blocked_qr ( void );
/*
 * DESCRIPTION :
 *   Generates a random complex matrix to test the blocked
 *   computation of the QR decomposition. */

#endif
