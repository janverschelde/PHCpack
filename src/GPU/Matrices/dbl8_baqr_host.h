/* The file dbl8_baqr_host.h specifies functions on the host for the
 * blocked accelerated qr decomposition in octo double precision. */

#ifndef __dbl8_baqr_host_h__
#define __dbl8_baqr_host_h__

void CPU_dbl8_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double **Vhihihi, double **Vlohihi, double **Vhilohi, double **Vlolohi,
   double **Vhihilo, double **Vlohilo, double **Vhilolo, double **Vlololo,
   double **Whihihi, double **Wlohihi, double **Whilohi, double **Wlolohi,
   double **Whihilo, double **Wlohilo, double **Whilolo, double **Wlololo );
/*
 * DESCRIPTION :
 *   Computes the W in the WY representation of the Householder
 *   transformations defined by V and B, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   Bhihihi  are the highest doubles of the beta computed by house;
 *   Blohihi  are the second highest doubles of the beta computed by house;
 *   Bhilohi  are the third highest doubles of the beta computed by house;
 *   Blolohi  are the fourth highest doubles of the beta computed by house;
 *   Bhihilo  are the fourth lowest doubles of the beta computed by house;
 *   Blohilo  are the third lowest doubles of the beta computed by house;
 *   Bhilolo  are the second lowest doubles of the beta computed by house;
 *   Blololo  are the lowest doubles of the beta computed by house;
 *   Vhihihi  are the highest doubles of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vlohihi  are the second highest doubles of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vhilohi  are the third highest doubles of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vlolohi  are the fourth highest doubles of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vhihilo  are the fourth lowest double of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vlohilo  are the third lowest double of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vhilolo  are the second lowest double of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vlololo  are the lowest double of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Whihihi  has space for ncols columns with rows from 0 to nrows-1;
 *   Wlohihi  has space for ncols columns with rows from 0 to nrows-1;
 *   Whilohi  has space for ncols columns with rows from 0 to nrows-1;
 *   Wlolohi  has space for ncols columns with rows from 0 to nrows-1;
 *   Whihilo  has space for ncols columns with rows from 0 to nrows-1;
 *   Wlohilo  has space for ncols columns with rows from 0 to nrows-1;
 *   Whilolo  has space for ncols columns with rows from 0 to nrows-1;
 *   Wlololo  has space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Whihihi  are the highest doubles of W in the WY representation;
 *   Wlohihi  are the second highest doubles of W;
 *   Whilohi  are the third highest doubles of W;
 *   Wlolohi  are the fourth highest doubles of W;
 *   Whihilo  are the fourth lowest doubles of W;
 *   Wlohilo  are the third lowest doubles of W;
 *   Whilolo  are the second lowest doubles of W;
 *   Wlololo  are the lowest doubles of W. */

void CPU_cmplx8_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double **Vrehihihi, double **Vrelohihi,
   double **Vrehilohi, double **Vrelolohi,
   double **Vrehihilo, double **Vrelohilo,
   double **Vrehilolo, double **Vrelololo,
   double **Vimhihihi, double **Vimlohihi,
   double **Vimhilohi, double **Vimlolohi,
   double **Vimhihilo, double **Vimlohilo,
   double **Vimhilolo, double **Vimlololo,
   double **Wrehihihi, double **Wrelohihi,
   double **Wrehilohi, double **Wrelolohi,
   double **Wrehihilo, double **Wrelohilo,
   double **Wrehilolo, double **Wrelololo,
   double **Wimhihihi, double **Wimlohihi,
   double **Wimhilohi, double **Wimlolohi,
   double **Wimhihilo, double **Wimlohilo,
   double **Wimhilolo, double **Wimlololo );
/*
 * DESCRIPTION :
 *   Computes the W in the WY representation of the Householder
 *   transformations defined by V and B, on complex data.
 *
 * ON ENTRY :
 *   nrows     number of rows in the matrices V, Y, and W;
 *   ncols     equals the size of one tile, or equivalently,
 *             is the number of elements in B,
 *             and the number of columns in V, Y, and W;
 *   Bhihihi   are the highest doubles of the beta computed by house;
 *   Blohihi   are the second highest doubles of the beta computed by house;
 *   Bhilohi   are the third highest doubles of the beta computed by house;
 *   Blolohi   are the fourth highest doubles of the beta computed by house;
 *   Bhihilo   are the fourth lowest doubles of the beta computed by house;
 *   Blohilo   are the third lowest doubles of the beta computed by house;
 *   Bhilolo   are the second lowest doubles of the beta computed by house;
 *   Blololo   are the lowest doubles of the beta computed by house;
 *   Vrehihihi are the highest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrelohihi are the second highest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrehilohi are the third highest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrelolohi are the fourth highest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrehihilo are the fourth lowest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrelohilo are the third lowest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrehilolo are the second lowest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrelololo are the lowest doubles of the real parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimhihihi are the highest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimlohihi are the second highest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimhilohi are the third highest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimlolohi are the fourth highest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimhihilo are the fourth lowest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimlohilo are the third lowest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimhilolo are the second lowest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimlololo are the lowest doubles of the imaginary parts of the
 *             Householder vectors, with zeros inserted so V is trapezoidal;
 *   Wrehihihi has space for ncols columns with rows from 0 to nrows-1;
 *   Wrelohihi has space for ncols columns with rows from 0 to nrows-1;
 *   Wrehilohi has space for ncols columns with rows from 0 to nrows-1;
 *   Wrelolohi has space for ncols columns with rows from 0 to nrows-1;
 *   Wrehihilo has space for ncols columns with rows from 0 to nrows-1;
 *   Wrelohilo has space for ncols columns with rows from 0 to nrows-1;
 *   Wrehilolo has space for ncols columns with rows from 0 to nrows-1;
 *   Wrelololo has space for ncols columns with rows from 0 to nrows-1;
 *   Wimhihihi has space for ncols columns with rows from 0 to nrows-1;
 *   Wimlohihi has space for ncols columns with rows from 0 to nrows-1;
 *   Wimhilohi has space for ncols columns with rows from 0 to nrows-1;
 *   Wimlolohi has space for ncols columns with rows from 0 to nrows-1;
 *   Wimhihilo has space for ncols columns with rows from 0 to nrows-1;
 *   Wimlohilo has space for ncols columns with rows from 0 to nrows-1;
 *   Wimhilolo has space for ncols columns with rows from 0 to nrows-1;
 *   Wimlololo has space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Wrehihihi are the highest doubles of the real parts of W;
 *   Wrelohihi are the second highest doubles of the real parts of W;
 *   Wrehilohi are the third highest doubles of the real parts of W;
 *   Wrelolohi are the fourth highest doubles of the real parts of W;
 *   Wrehihilo are the fourth lowest doubles of the real parts of W;
 *   Wrelohilo are the third lowest doubles of the real parts of W;
 *   Wrehilolo are the second lowest doubles of the real parts of W;
 *   Wrelololo are the lowest doubles of the real parts of W;
 *   Wimhihihi are the highest doubles of the imaginary parts of W;
 *   Wimlohihi are the second highest doubles of the imaginary parts of W;
 *   Wimhilohi are the third highest doubles of the imaginary parts of W;
 *   Wimlolohi are the fourth highest doubles of the imaginary parts of W;
 *   Wimhihilo are the fourth lowest doubles of the imaginary parts of W;
 *   Wimlohilo are the third lowest doubles of the imaginary parts of W;
 *   Wimhilolo are the second lowest doubles of the imaginary parts of W;
 *   Wimlololo are the lowest doubles of the imaginary parts of W. */

void CPU_dbl8_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Chihihi, double **Clohihi, double **Chilohi, double **Clolohi,
   double **Chihilo, double **Clohilo, double **Chilolo, double **Clololo,
   double **Yhihihi, double **Ylohihi, double **Yhilohi, double **Ylolohi,
   double **Yhihilo, double **Ylohilo, double **Yhilolo, double **Ylololo,
   double **Whihihi, double **Wlohihi, double **Whilohi, double **Wlolohi,
   double **Whihilo, double **Wlohilo, double **Whilolo, double **Wlololo,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the blocked Householder update to C,
 *   as C = C + (Y*W^T)*C, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix C;
 *   ncols    number of columns in the matrix C;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Chihihi  are the highest doubles of C;
 *   Clohihi  are the second highest doubles of C;
 *   Chilohi  are the third highest doubles of C;
 *   Clolohi  are the fourth highest doubles of C;
 *   Chihilo  are the fourth lowest doubles of C;
 *   Clohilo  are the third lowest doubles of C;
 *   Chilolo  are the second lowest doubles of C;
 *   Clololo  are the lowest doubles of C;
 *   Yhihihi  are the highest doubles of the Householder vectors;
 *   Ylohihi  are the second highest doubles of the Householder vectors;
 *   Yhilohi  are the third highest doubles of the Householder vectors;
 *   Ylohihi  are the fourth highest doubles of the Householder vectors;
 *   Yhihilo  are the fourth lowest doubles of the Householder vectors;
 *   Ylohilo  are the third lowest doubles of the Householder vectors;
 *   Yhilolo  are the second lowest doubles of the Householder vectors;
 *   Ylololo  are the lowest doubles of the Householder vectors;
 *   Whihihi  are the highest doubles of W in the WY representation;
 *   Wlohihi  are the second highest doubles of W in the WY representation;
 *   Whilohi  are the third highest doubles of W in the WY representation;
 *   Wlolohi  are the fourth highest doubles of W in the WY representation;
 *   Whihilo  are the fourth lowest doubles of W in the WY representation;
 *   Wlohilo  are the third lowest doubles of W in the WY representation;
 *   Whilolo  are the second lowest doubles of W in the WY representation;
 *   Wlololo  are the lowest doubles of W in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Chihihi  are the highest doubles of the update
 *            with the Householder matrix;
 *   Clohihi  are the second highest doubles of the update;
 *   Chilohi  are the third highest doubles of the update;
 *   Clolohi  are the fourth highest doubles of the update;
 *   Chihilo  are the fourth lowest doubles of the update;
 *   Clohilo  are the third lowest doubles of the update;
 *   Chilolo  are the second lowest doubles of the update;
 *   Clololo  are the lowest doubles of the update. */

void CPU_cmplx8_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Crehihihi, double **Crelohihi,
   double **Crehilohi, double **Crelolohi,
   double **Crehihilo, double **Crelohilo,
   double **Crehilolo, double **Crelololo,
   double **Cimhihihi, double **Cimlohihi,
   double **Cimhilohi, double **Cimlolohi,
   double **Cimhihilo, double **Cimlohilo,
   double **Cimhilolo, double **Cimlololo,
   double **Yrehihihi, double **Yrelohihi,
   double **Yrehilohi, double **Yrelolohi,
   double **Yrehihilo, double **Yrelohilo,
   double **Yrehilolo, double **Yrelololo,
   double **Yimhihihi, double **Yimlohihi,
   double **Yimhilohi, double **Yimlolohi,
   double **Yimhihilo, double **Yimlohilo,
   double **Yimhilolo, double **Yimlololo,
   double **Wrehihihi, double **Wrelohihi,
   double **Wrehilohi, double **Wrelolohi,
   double **Wrehihilo, double **Wrelohilo,
   double **Wrehilolo, double **Wrelololo,
   double **Wimhihihi, double **Wimlohihi,
   double **Wimhilohi, double **Wimlolohi,
   double **Wimhihilo, double **Wimlohilo,
   double **Wimhilolo, double **Wimlololo,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the blocked Householder update to C,
 *   as C = C + (Y*W^T)*C, on complex data.
 *
 * ON ENTRY :
 *   nrows     number of rows in the matrix C;
 *   ncols     number of columns in the matrix C;
 *   szt       size of each block;
 *   idx       index of the current block;
 *   Crehihihi are the highest doubles of the real parts of C;
 *   Crelohihi are the second highest doubles of the real parts of C;
 *   Crehilohi are the third highest doubles of the real parts of C;
 *   Crelolohi are the fourth highest doubles of the real parts of C;
 *   Crehihilo are the fourth lowest doubles of the real parts of C;
 *   Crelohilo are the third lowest doubles of the real parts of C;
 *   Crehilolo are the second lowest doubles of the real parts of C;
 *   Crelololo are the lowest doubles of the real parts of C;
 *   Cimhihihi are the highest doubles of the imaginary parts of C;
 *   Cimlohihi are the second highest doubles of the imaginary parts of C;
 *   Cimhilohi are the third highest doubles of the imaginary parts of C;
 *   Cimlolohi are the fourth highest doubles of the imaginary parts of C;
 *   Cimhihilo are the fourth lowest doubles of the imaginary parts of C;
 *   Cimlohilo are the third lowest doubles of the imaginary parts of C;
 *   Cimhilolo are the second lowest doubles of the imaginary parts of C;
 *   Cimlololo are the lowest doubles of the imaginary parts of C;
 *   Yrehihihi are the highest doubles of the real parts
 *             of the Householder vectors Y;
 *   Yrelohihi are the second highest doubles of the real parts of Y;
 *   Yrehilohi are the third highest doubles of the real parts of Y;
 *   Yrelolohi are the fourth highest doubles of the real parts of Y;
 *   Yrehihilo are the fourth lowest doubles of the real parts of Y;
 *   Yrelohilo are the third lowest doubles of the real parts of Y;
 *   Yrehilolo are the second lowest doubles of the real parts of Y;
 *   Yrelololo are the lowest doubles of the real parts of Y;
 *   Yimhihihi are the highest doubles of the imaginary parts of Y;
 *   Yimlohihi are the second highest doubles of the imaginary parts of Y;
 *   Yimhilohi are the third highest doubles of the imaginary parts of Y;
 *   Yimlolohi are the fourth highest doubles of the imaginary parts of Y;
 *   Yimhihilo are the fourth lowest doubles of the imaginary parts of Y;
 *   Yimlohilo are the third lowest doubles of the imaginary parts of Y;
 *   Yimhilolo are the second lowest doubles of the imaginary parts of Y;
 *   Yimlololo are the lowest doubles of the imaginary parts of Y;
 *   Wrehihihi are the highest doubles of the real parts of W;
 *   Wrelohihi are the second highest doubles of the real parts of W;
 *   Wrehilohi are the third highest doubles of the real parts of W;
 *   Wrelolohi are the fourth highest doubles of the real parts of W;
 *   Wrehihilo are the fourth lowest doubles of the real parts of W;
 *   Wrelohilo are the third lowest doubles of the real parts of W;
 *   Wrehilolo are the second lowest doubles of the real parts of W;
 *   Wrelololo are the lowest doubles of the real parts of W;
 *   Wimhihihi are the highest doubles of the imaginary parts of W;
 *   Wimlohihi are the second highest doubles of the imaginary parts of W;
 *   Wimhilohi are the third highest doubles of the imaginary parts of W;
 *   Wimlolohi are the fourth highest doubles of the imaginary parts of W;
 *   Wimhihilo are the fourth lowest doubles of the imaginary parts of W;
 *   Wimlohilo are the third lowest doubles of the imaginary parts of W;
 *   Wimhilolo are the second lowest doubles of the imaginary parts of W;
 *   Wimlololo are the lowest doubles of the imaginary parts of W;
 *   verbose   is the verbose flag.
 *
 * ON RETURN :
 *   Crehihihi are the highest doubles of the real parts of the update;
 *   Crelohihi are the second highest doubles of the real parts of the update;
 *   Crehilohi are the third highest doubles of the real parts of the update;
 *   Crelolohi are the fourth highest doubles of the real parts of the update;
 *   Crehihilo are the fourth lowest doubles of the real parts of the update;
 *   Crelohilo are the third lowest doubles of the real parts of the update;
 *   Crehilolo are the second lowest doubles of the real parts of the update;
 *   Crelololo are the lowest doubles of the real parts of the update;
 *   Cimhihihi are the highest doubles of the imaginary parts of the update;
 *   Cimlohihi are the second highest doubles of the imag parts of the update;
 *   Cimhilohi are the third highest doubles of the imag parts of the update;
 *   Cimlolohi are the fourth highest doubles of the imag parts of the update;
 *   Cimhihilo are the fourth lowest doubles of the imag parts of the update;
 *   Cimlohilo are the third lowest doubles of the imag parts of the update;
 *   Cimhilolo are the second lowest doubles of the imag parts of the update;
 *   Cimlololo are the lowest doubles of the imaginary parts of the update. */

void CPU_dbl8_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Yhihihi, double **Ylohihi, double **Yhilohi, double **Ylolohi,
   double **Yhihilo, double **Ylohilo, double **Yhilolo, double **Ylololo,
   double **Whihihi, double **Wlohihi, double **Whilohi, double **Wlolohi,
   double **Whihilo, double **Wlohilo, double **Whilolo, double **Wlololo,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix Q;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Qhihihi  are the highest doubles of Q;
 *   Qlohihi  are the second highest doubles of Q;
 *   Qhilohi  are the third highest doubles of Q;
 *   Qlolohi  are the fourth highest doubles of Q;
 *   Qhihilo  are the fourth lowest doubles of Q;
 *   Qlohilo  are the third lowest doubles of Q;
 *   Qhilolo  are the second lowest doubles of Q;
 *   Qlololo  are the lowest doubles of Q;
 *   Yhihihi  are the highest doubles of the Householder vectors Y;
 *   Ylohihi  are the second highest doubles of Y;
 *   Yhilohi  are the third highest doubles of Y;
 *   Ylolohi  are the fourth highest doubles of Y;
 *   Yhihilo  are the fourth lowest doubles of Y;
 *   Ylohilo  are the third lowest doubles of Y;
 *   Yhilolo  are the second lowest doubles of Y;
 *   Ylololo  are the lowest doubles of Y;
 *   Whihihi  are the highest doubles of the W in the WY representation;
 *   Wlohihi  are the second highest doubles of W;
 *   Whilohi  are the third highest doubles of W;
 *   Wlolohi  are the fourth highest doubles of W;
 *   Whihilo  are the fourth lowest doubles of W;
 *   Wlohilo  are the third lowest doubles of W;
 *   Whilolo  are the second lowest doubles of W;
 *   Wlololo  are the lowest doubles of W;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhihihi  are the highest doubles of the updated Q;
 *   Qlohihi  are the second highest doubles of the updated Q;
 *   Qhilohi  are the third highest doubles of the updated Q;
 *   Qlolohi  are the fourth highest doubles of the updated Q;
 *   Qhihilo  are the fourth lowest doubles of the updated Q;
 *   Qlohilo  are the third lowest doubles of the updated Q;
 *   Qhihilo  are the second lowest doubles of the updated Q;
 *   Qlololo  are the lowest doubles of the updated Q. */

void CPU_cmplx8_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Yrehihihi, double **Yrelohihi,
   double **Yrehilohi, double **Yrelolohi,
   double **Yrehihilo, double **Yrelohilo,
   double **Yrehilolo, double **Yrelololo,
   double **Yimhihihi, double **Yimlohihi,
   double **Yimhilohi, double **Yimlolohi,
   double **Yimhihilo, double **Yimlohilo,
   double **Yimhilolo, double **Yimlololo,
   double **Wrehihihi, double **Wrelohihi,
   double **Wrehilohi, double **Wrelolohi,
   double **Wrehihilo, double **Wrelohilo,
   double **Wrehilolo, double **Wrelololo,
   double **Wimhihihi, double **Wimlohihi,
   double **Wimhilohi, double **Wimlolohi,
   double **Wimhihilo, double **Wimlohilo,
   double **Wimhilolo, double **Wimlololo, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on complex data.
 *
 * ON ENTRY :
 *   dim       number of rows and columns in the matrix Q;
 *   szt       size of each block;
 *   idx       index of the current block;
 *   Qrehihihi are the highest doubles of the real parts of Q;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   Yrehihihi are the highest doubles of the real parts
 *             of the Householder vectors Y;
 *   Yrelohihi are the second highest doubles of the real parts of Y;
 *   Yrehilohi are the third highest doubles of the real parts of Y;
 *   Yrelolohi are the fourth highest doubles of the real parts of Y;
 *   Yrehihilo are the fourth lowest doubles of the real parts of Y;
 *   Yrelohilo are the third lowest doubles of the real parts of Y;
 *   Yrehilolo are the second lowest doubles of the real parts of Y;
 *   Yrelololo are the lowest doubles of the real parts of Y;
 *   Yimhihihi are the highest doubles of the imaginary parts of Y;
 *   Yimlohihi are the second highest doubles of the imaginary parts of Y;
 *   Yimhilohi are the third highest doubles of the imaginary parts of Y;
 *   Yimlolohi are the fourth highest doubles of the imaginary parts of Y;
 *   Yimhihilo are the fourth lowest doubles of the imaginary parts of Ys;
 *   Yimlohilo are the third lowest doubles of the imaginary parts of Ys;
 *   Yimhilolo are the second lowest doubles of the imaginary parts of Ys;
 *   Yimlololo are the lowest doubles of the imaginary parts of Ys;
 *   Wrehihihi are the highest doubles of the real parts of W;
 *   Wrelohihi are the second highest doubles of the real parts of W;
 *   Wrehilohi are the third highest doubles of the real parts of W;
 *   Wrelolohi are the fourth highest doubles of the real parts of W;
 *   Wrehihilo are the fourth lowest doubles of the real parts of W;
 *   Wrelohilo are the thirdd lowest doubles of the real parts of W;
 *   Wrehilolo are the second lowest doubles of the real parts of W;
 *   Wrelololo are the lowest doubles of the real parts of W;
 *   Wimhihihi are the highest doubles of the imaginary parts of W;
 *   Wimlohihi are the second highest doubles of the imaginary parts of W;
 *   Wimhilohi are the third highest doubles of the imaginary parts of W;
 *   Wimlolohi are the fourth highest doubles of the imaginary parts of W;
 *   Wimhihilo are the fourth lowest doubles of the imaginary parts of W;
 *   Wimlohilo are the third lowest doubles of the imaginary parts of W;
 *   Wimhilolo are the second lowest doubles of the imaginary parts of W;
 *   Wimlololo are the lowest doubles of the imaginary parts of W;
 *   verbose   is the verbose flag.
 *
 * ON RETURN :
 *   Qrehihihi are the highest doubles of the real parts of the updated Q;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of the updated Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q. */

void CPU_dbl8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *lapsed, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies Householder transformations in a blocked manner
 *   to compute a QR decomposition of A, on real data.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   szt      size of each block;
 *   nbt      number of tiles, ncols = szt*nbt;
 *   Ahihihi  are the highest doubles of A,
 *            stored as nrows arrays of ncols numbers;
 *   Alohihi  are the second highest doubles of A;
 *   Ahilohi  are the third highest doubles of A;
 *   Alolohi  are the fourth highest doubles of A;
 *   Ahihilo  are the fourth lowest doubles of A;
 *   Alohilo  are the third lowest doubles of A;
 *   Ahilolo  are the second lowest doubles of A;
 *   Alololo  are the lowest doubles of A;
 *   Qhihihi  has space for an nrows-by-nrows matrix;
 *   Qlohihi  has space for an nrows-by-nrows matrix;
 *   Qhilohi  has space for an nrows-by-nrows matrix;
 *   Qlolohi  has space for an nrows-by-nrows matrix;
 *   Qhihilo  has space for an nrows-by-nrows matrix;
 *   Qlohilo  has space for an nrows-by-nrows matrix;
 *   Qhilolo  has space for an nrows-by-nrows matrix;
 *   Qlololo  has space for an nrows-by-nrows matrix;
 *   Rhihihi  has space for an nrows-by-ncols matrix;
 *   Rlohihi  has space for an nrows-by-ncols matrix;
 *   Rhilohi  has space for an nrows-by-ncols matrix;
 *   Rlolohi  has space for an nrows-by-ncols matrix;
 *   Rhihilo  has space for an nrows-by-ncols matrix;
 *   Rlohilo  has space for an nrows-by-ncols matrix;
 *   Rhilolo  has space for an nrows-by-ncols matrix;
 *   Rlololo  has space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhihihi  are the highest doubles of Q,
 *            Q is orthogonal and transpose(Q)*A = R;
 *   Qlohihi  are the second highest doubles of Q;
 *   Qhilohi  are the third highest doubles of Q;
 *   Qlolohi  are the fourth highest doubles of Q;
 *   Qhihilo  are the fourth lowest doubles of Q;
 *   Qlohilo  are the third lowest doubles of Q;
 *   Qhilolo  are the second lowest doubles of Q;
 *   Qlololo  are the lowest doubles of Q;
 *   Rhihihi  are the highest doubles of R,
 *            the reduced upper triangular form of A;
 *   Rlohihi  are the second highest doubles of R;
 *   Rhilohi  are the third highest doubles of R;
 *   Rlolohi  are the fourth highest doubles of R;
 *   Rhihilo  are the fourth lowest doubles of R;
 *   Rlohilo  are the third lowest doubles of R;
 *   Rhilolo  are the second lowest doubles of R;
 *   Rlololo  are the lowest doubles of R;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *lapsec, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies Householder transformations in a blocked manner
 *   to compute a QR decomposition of A, on complex data.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows     number of rows of A;
 *   ncols     number of columns of A;
 *   szt       size of each block;
 *   nbt       number of tiles, ncols = szt*nbt;
 *   Arehihihi are the highest doubles of the real parts of A,
 *             stored as nrows arrays of ncols numbers;
 *   Arelohihi are the second highest doubles of the real parts of A;
 *   Arehilohi are the third highest doubles of the real parts of A;
 *   Arelolohi are the fourth highest doubles of the real parts of A;
 *   Arehihilo are the fourth lowest doubles of the real parts of A;
 *   Arelohilo are the third lowest doubles of the real parts of A;
 *   Arehilolo are the second lowest doubles of the real parts of A;
 *   Arelololo are the lowest doubles of the real parts of A;
 *   Aimhihihi are the highest doubles of the imaginary parts of A;
 *   Aimlohihi are the second highest doubles of the imaginary parts of A;
 *   Aimhilohi are the third highest doubles of the imaginary parts of A;
 *   Aimlolohi are the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo are the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo are the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo are the second lowest doubles of the imaginary parts of A;
 *   Aimlololo are the lowest doubles of the imaginary parts of A;
 *   Qrehihihi has space for an nrows-by-nrows matrix;
 *   Qrelohihi has space for an nrows-by-nrows matrix;
 *   Qrehilohi has space for an nrows-by-nrows matrix;
 *   Qrelolohi has space for an nrows-by-nrows matrix;
 *   Qrehihilo has space for an nrows-by-nrows matrix;
 *   Qrelohilo has space for an nrows-by-nrows matrix;
 *   Qrehilolo has space for an nrows-by-nrows matrix;
 *   Qrelololo has space for an nrows-by-nrows matrix;
 *   Qimhihihi has space for an nrows-by-nrows matrix;
 *   Qimlohihi has space for an nrows-by-nrows matrix;
 *   Qimhilohi has space for an nrows-by-nrows matrix;
 *   Qimlolohi has space for an nrows-by-nrows matrix;
 *   Qimhihilo has space for an nrows-by-nrows matrix;
 *   Qimlohilo has space for an nrows-by-nrows matrix;
 *   Qimhilolo has space for an nrows-by-nrows matrix;
 *   Qimlololo has space for an nrows-by-nrows matrix;
 *   Rrehihihi has space for an nrows-by-ncols matrix;
 *   Rrelohihi has space for an nrows-by-ncols matrix;
 *   Rrehilohi has space for an nrows-by-ncols matrix;
 *   Rrelolohi has space for an nrows-by-ncols matrix;
 *   Rrehihilo has space for an nrows-by-ncols matrix;
 *   Rrelohilo has space for an nrows-by-ncols matrix;
 *   Rrehilolo has space for an nrows-by-ncols matrix;
 *   Rrelololo has space for an nrows-by-ncols matrix;
 *   Rimhihihi has space for an nrows-by-ncols matrix;
 *   Rimlohihi has space for an nrows-by-ncols matrix;
 *   Rimhilohi has space for an nrows-by-ncols matrix;
 *   Rimlolohi has space for an nrows-by-ncols matrix;
 *   Rimhihilo has space for an nrows-by-ncols matrix;
 *   Rimlohilo has space for an nrows-by-ncols matrix;
 *   Rimhilolo has space for an nrows-by-ncols matrix;
 *   Rimlololo has space for an nrows-by-ncols matrix;
 *   verbose   is the verbose flag.
 *
 * ON RETURN :
 *   Qrehihihi are the highest doubles of the real parts of Q,
 *             Q is orthogonal and transpose(Q)*A = R;
 *   Qrelohihi are the second highest doubles of the real parts of Q;
 *   Qrehilohi are the third highest doubles of the real parts of Q;
 *   Qrelolohi are the fourth highest doubles of the real parts of Q;
 *   Qrehihilo are the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo are the third lowest doubles of the real parts of Q;
 *   Qrehilolo are the second lowest doubles of the real parts of Q;
 *   Qrelololo are the lowest doubles of the real parts of Q;
 *   Qimhihihi are the highest doubles of the imaginary parts of Q;
 *   Qimlohihi are the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi are the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi are the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo are the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo are the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo are the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo are the lowest doubles of the imaginary parts of Q;
 *   Rrehihihi are the highest doubles of the real parts of R,
 *             the reduced upper triangular form of A;
 *   Rrelohihi are the second highest doubles of the real parts of R;
 *   Rrehilohi are the third highest doubles of the real parts of R;
 *   Rrelolohi are the fourth highest doubles of the real parts of R;
 *   Rrehihilo are the fourth lowest doubles of the real parts of R;
 *   Rrelohilo are the third lowest doubles of the real parts of R;
 *   Rrehilolo are the second lowest doubles of the real parts of R;
 *   Rrelololo are the lowest doubles of the real parts of R;
 *   Rimhihihi are the highest doubles of the imaginary parts of R;
 *   Rimlohihi are the second highest doubles of the imaginary parts of R;
 *   Rimhilohi are the third highest doubles of the imaginary parts of R;
 *   Rimlolohi are the fourth highest doubles of the imaginary parts of R;
 *   Rimhihilo are the fourth lowest doubles of the imaginary parts of R;
 *   Rimlohilo are the third lowest doubles of the imaginary parts of R;
 *   Rimhilolo are the second lowest doubles of the imaginary parts of R;
 *   Rimlololo are the lowest doubles of the imaginary parts of R;
 *   lapsec    elapsed time in seconds. */

#endif
