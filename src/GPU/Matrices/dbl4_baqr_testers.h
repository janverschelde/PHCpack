// The file dbl4_baqr_testers.h specifies test functions on
// blocked accelerated QR decomposition in quad double precision.

#ifndef __dbl4_baqr_testers_h__
#define __dbl4_baqr_testers_h__

void test_real4_blocked_qr
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode );
/*
 * DESCRIPTION :
 *   Generates a random real matrix to test the blocked
 *   computation of the QR decomposition.
 *
 * ON ENTRY :
 *   seed     the seed for the random number generator (0 for time);
 *   szt      size of one tile;
 *   nbt      number of tiles, szt*nbt equals the number of columns;
 *   nrows    number of rows >= number of columns;
 *   vrb      the verbose level;
 *   mode     execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU). */

void test_cmplx4_blocked_qr
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode );
/*
 * DESCRIPTION :
 *   Generates a random complex matrix to test the blocked
 *   computation of the QR decomposition.
 *
 * ON ENTRY :
 *   seed     the seed for the random number generator (0 for time);
 *   szt      size of one tile;
 *   nbt      number of tiles, szt*nbt equals the number of columns;
 *   nrows    number of rows >= number of columns;
 *   vrb      the verbose level;
 *   mode     execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU). */

#endif
