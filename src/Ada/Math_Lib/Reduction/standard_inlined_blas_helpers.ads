with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;

package Standard_Inlined_BLAS_Helpers is

-- DESCRIPTION :
--   Functions and procedures to help the inlined SVD computation,
--   with complex arithmetic inlined with the vector arithmetic.

  function cdabs ( zr,zi : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the modulus of the complex number z,
  --   with real part in zr and imaginary part in zi.

  procedure csign ( pr,pi,qr,qi : in double_float; zr,zi : out double_float );

  -- DESCRIPTION :
  --   Given are two complex numbers p and q,
  --   as pairs of real and imaginary parts (pr, pi), (qr, qi).
  --   Returns in zr and zi respectively the real and imaginary
  --   parts of(cdabs(pr,pi)/cdabs(qr,qi)) multiplied with q.

  function dznrm2 ( n : integer32;
                    xre : Standard_Floating_Vectors.Link_to_Vector;
                    xim : Standard_Floating_Vectors.Link_to_Vector;
                    ind,incx : integer32 ) return double_float;

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector given as a pair of vectors,
  --   as a vector of its real parts and a vector of its imaginary parts,
  --   starting at xre(ind), xim(ind), with n steps of increment incx.

  -- REQUIRED :
  --   xre'range = xim'range and xre'last >= ind + (n-1)*incx.

  function dznrm2 ( n : integer32;
                    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
                    row,col,incx : integer32 ) return double_float;

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a column of a matrix given
  --   as columns of real and imaginary parts.

  -- REQUIRED : col is in xre'range and col is in xim'range,
  --   xre(col)'range = xim(col)'range 
  --   and xre(col)'last >= row + (n-1)*incx.

  -- ON ENTRY :
  --   n        number of steps done in the matrix;
  --   rvv      real parts of the columns of the matrix;
  --   ivv      imaginary parts of the columns of the matrix;
  --   row      start row in the column col of the matrix;
  --   col      index to the column in the matrix;
  --   incx     increment in the matrix.

  procedure zscal ( n : in integer32; zr,zi : in double_float;
                    xre : in Standard_Floating_Vectors.Link_to_Vector;
                    xim : in Standard_Floating_Vectors.Link_to_Vector;
                    ind,incx : in integer32 );

  -- DESCRIPTION :
  --   Scales the vector given in splitted form by multiplication with 
  --   a complex number, given by its real and imaginary part.

  -- REQUIRED : xre'range = xim'range and xre'last >= ind + (n-1)*incx.

  -- ON ENTRY :
  --   n       number of steps;
  --   zr      real part of the multiplier;
  --   zi      imaginary part of the multiplier;
  --   xre     real parts of a complex vector;
  --   xim     imaginary parts of a complex vector;
  --   ind     start index in xre and xim;
  --   incx    increment in the vector.

  -- ON RETURN :
  --   xre     scaled real parts of the vector;
  --   xim     scaled imaginary parts of the vector.

  procedure zscal ( n : in integer32; zr,zi : in double_float;
                    rvv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    ivv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    row,col,incx : in integer32 );

  -- DESCRIPTION :
  --   Scales the column of a matrix given in splitted form by 
  --   multiplication with  a complex number,
  --   given by its real and imaginary part.

  -- REQUIRED : col is in xre'range and col is in xim'range,
  --   xre(col)'range = xim(col)'range 
  --   and xre(col)'last >= row + (n-1)*incx.

  -- ON ENTRY :
  --   n       number of steps;
  --   zr      real part of the multiplier;
  --   zi      imaginary part of the multiplier;
  --   rvv     real parts of a complex column;
  --   ivv     imaginary parts of a complex colum;
  --   row     start index in rvv(col) and ivv(col);
  --   incx    increment in the vector.

  -- ON RETURN :
  --   rvv     scaled real parts of the column rvv(col);
  --   ivv     scaled imaginary parts of the column ivv(col).

-- UPDATE COLUMN OF A MATRIX WITH A MULTIPLE OF A VECTOR :

  procedure zaxpy ( n : in integer32; zr,zi : in double_float;
                    xre : in Standard_Floating_Vectors.Link_to_Vector;
                    xim : in Standard_Floating_Vectors.Link_to_Vector;
                    ind,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32 );

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at ind and in y
  --   at (rwy,cly), using increments incx and incy to advance.

  -- ON ENTRY :
  --   n        number of steps;
  --   zr       real part of the multiplier;
  --   zi       imaginary part of the multiplier;
  --   xre      real parts of a complex vector;
  --   xim      imaginary parts of a complex vector;
  --   ind      start index in the vectors xre and xim;
  --   incx     increment of xre and xim;
  --   yrv      real parts of the columns of a complex matrix;
  --   yiv      imaginary parts of the columns of a complex matrix;
  --   rwy      start row in the columns on yrv and yiv;
  --   cly      column index in yrv and yiv;
  --   incy     increment for the rows of yrv and yiv.

  -- ON RETURN :
  --   yrv      updated real parts of the column yrv(cly);
  --   yiv      updated imaginary parts of the column yrv(cly).

-- UPDATE VECTOR WITH A MULTIPLE OF A COLUMN OF A MATRIX :

  procedure zaxpy ( n : in integer32; zr,zi : in double_float;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yre : in  Standard_Floating_Vectors.Link_to_Vector;
                    yim : in  Standard_Floating_Vectors.Link_to_Vector;
                    ind,incy : in integer32 );

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at (rwx,clx)
  --   and in y at ind, using increments incx and incy to advance.

  -- ON ENTRY :
  --   n        number of steps;
  --   zr       real part of the multiplier;
  --   zi       imaginary part of the multiplier;
  --   xrv      real parts of the columns of a complex matrix;
  --   xiv      imaginary parts of the columns of a complex matrix;
  --   rwx      start row in the columns on xrv and xiv;
  --   clx      column index in xrv and xiv;
  --   incx     increment for the rows of xrv and xiv;
  --   yre      real parts of a complex vector;
  --   yim      imaginary parts of a complex vector;
  --   ind      start index in the vectors yre and yim;
  --   incy     increment of yre and yim.

  -- ON RETURN :
  --   yre      updated real parts of a complex vector;
  --   yie      updated imaginary parts of a complex vector.

-- UPDATE COLUMN WITH A MULTIPLE OF A COLUMN OF A MATRIX :

  procedure zaxpy ( n : in integer32; zr,zi : in double_float;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32 );

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x and y at the
  --   respective (rows,columns): (rwx,clx) and (rwy,cly), using
  --   increments incx and incy to advance in the rows.

  -- ON ENTRY :
  --   n        number of steps;
  --   zr       real part of the multiplier;
  --   zi       imaginary part of the multiplier;
  --   xrv      real parts of the columns of a complex matrix;
  --   xiv      imaginary parts of the columns of a complex matrix;
  --   rwx      start row in the columns on xrv and xiv;
  --   clx      column index in xrv and xiv;
  --   incx     increment for the rows of xrv and xiv;
  --   yrv      real parts of the columns of a complex matrix;
  --   yiv      imaginary parts of the columns of a complex matrix;
  --   rwy      start row in the columns on yrv and yiv;
  --   cly      column index in yrv and yiv;
  --   incy     increment for the rows of yrv and yiv.

  -- ON RETURN :
  --   yrv      updated real parts of a the column yrv(cly);
  --   yiv      updated imaginary parts of a the column yrv(cly).

-- DOT PRODUCT OF TWO COLUMNS OF A MATRIX :

  procedure zdotc ( n : in integer32;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32;
                    zr,zi : out double_float );

  -- DESCRIPTION :
  --   Returns the dot product of two vectors in two columns of matrices
  --   x and y, starting at rows rwx and rwy respectively, with respective
  --   increments in incx and incy.

  -- ON ENTRY :
  --   n        number of steps;
  --   xrv      real parts of the columns of a complex matrix;
  --   xiv      imaginary parts of the columns of a complex matrix;
  --   rwx      start row in the columns on xrv and xiv;
  --   clx      column index in xrv and xiv;
  --   incx     increment for the rows of xrv and xiv;
  --   yrv      real parts of the columns of a complex matrix;
  --   yiv      imaginary parts of the columns of a complex matrix;
  --   rwy      start row in the columns on yrv and yiv;
  --   cly      column index in yrv and yiv;
  --   incy     increment for the rows of yrv and yiv.

  -- ON RETURN :
  --   zr       real part of the dot product;
  --   zi       imaginary part of the dot product.

-- PLANE ROTATION ON TWO COLUMNS OF A MATRIX :

  procedure zdrot ( n : in integer32;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32;
                    c,s : in double_float );

  -- DESCRIPTION :
  --   Applies a plane rotation where the cos and sin are c and s
  --   and the vectors are in the columns of the matrices x and y,
  --   starting at (rwx,clx) and (rwy,cly) advancing in the rows with
  --   increments incx and incy respectively.

  -- ON ENTRY :
  --   n        number of steps;
  --   xrv      real parts of the columns of a complex matrix;
  --   xiv      imaginary parts of the columns of a complex matrix;
  --   rwx      start row in the columns on xrv and xiv;
  --   clx      column index in xrv and xiv;
  --   incx     increment for the rows of xrv and xiv;
  --   yrv      real parts of the columns of a complex matrix;
  --   yiv      imaginary parts of the columns of a complex matrix;
  --   rwy      start row in the columns on yrv and yiv;
  --   cly      column index in yrv and yiv;
  --   incy     increment for the rows of yrv and yiv;
  --   c        cosine of the rotation angle;
  --   s        sine of the rotation angle.

  -- ON RETURN :
  --   xrv      real parts of the updated column clx;
  --   xiv      imaginary parts of the updated column clx;
  --   yrv      real parts of the updated column cly;
  --   yiv      imaginary parts of the updated column cly.

  procedure zswap ( n : in integer32;
                    xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwx,clx,incx : in integer32;
                    yrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    yiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwy,cly,incy : in integer32 );

  -- DESCRIPTION :
  --   Interchanges two vectors in the columns of the matrices x and y,
  --   respectively starting at (rwx,clx) and (rwy,cly), and advancing
  --   in the rows with the respective increments incx and incy.

  -- ON ENTRY :
  --   n        number of steps;
  --   xrv      real parts of the columns of a complex matrix;
  --   xiv      imaginary parts of the columns of a complex matrix;
  --   rwx      start row in the columns on xrv and xiv;
  --   clx      column index in xrv and xiv;
  --   incx     increment for the rows of xrv and xiv;
  --   yrv      real parts of the columns of a complex matrix;
  --   yiv      imaginary parts of the columns of a complex matrix;
  --   rwy      start row in the columns on yrv and yiv;
  --   cly      column index in yrv and yiv;
  --   incy     increment for the rows of yrv and yiv.

  -- ON RETURN :
  --   xrv      real parts of the column clx swapped with yrv;
  --   xiv      imaginary parts of the column clx swapped with yiv;
  --   yrv      real parts of the column cly swapped with xrv;
  --   yiv      imaginary parts of the column cly swapped with xiv.

end Standard_Inlined_BLAS_Helpers;
