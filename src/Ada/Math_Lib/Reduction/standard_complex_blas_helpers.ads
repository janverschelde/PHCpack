with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Complex_Vectors;         use Standard_Complex_Vectors;
with Standard_Complex_Matrices;        use Standard_Complex_Matrices;

package Standard_Complex_BLAS_Helpers is

-- DESCRIPTION :
--   Functions and procedures are helping the SVD implementation,
--   translated from the linpack code.

  function Max0 ( a,b : integer32 ) return integer32;

  -- DESCRIPTION : returns the maximum of a and b.
  --   Note that this should correspond to the fortran max0 function.

  function dmax1 ( x1,x2 : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the maximum of x1 and x2.

  function dmax1 ( x1,x2,x3 : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the maximum of the three floating-point numbers.

  function dmax1 ( x1,x2,x3,x4 : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the maximum of the four floating-point numbers.

  function dmax1 ( x1,x2,x3,x4,x5 : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the maximum of the five floating-point numbers.

  function cabs1 ( z : Complex_Number ) return double_float;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of real and imaginary
  --   part of the complex number z.  Translation of
  --     complex*16 zdum
  --     double precision cabs1
  --     cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))

  function cdabs ( z : Complex_Number ) return double_float;

  -- DESCRIPTION :
  --   Computes the modulus of the complex number, let us hope this
  --   corresponds to the `cdabs' fortran function.

  function csign ( z1,z2 : Complex_Number ) return Complex_Number;

  -- DESCRIPTION : translated from
  --       csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2)) 

  function dsign ( a,b : double_float ) return double_float;

  -- DESCRIPTION :
  --   The implementation of this routine is written from web page
  --   documentation of sign...

  function dznrm2 ( n : integer32; x : Vector; ind,incx : integer32 )
                  return double_float;

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector x, starting at x(ind)
  --   and continueing n steps with increment incx.

  function dznrm2 ( n : integer32; x : Matrix; row,col,incx : integer32 )
                  return double_float;

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector x, starting at x(row,col)
  --   and continueing n steps with increment incx in the same column.

  procedure zscal ( n : in integer32; za : in Complex_Number;
                    zx : in out Vector; ind,incx : in integer32 );

  -- DESCRIPTION :
  --   Scales the vector starting at zx(ind) with the constant za.

  procedure zscal ( n : in integer32; za : in Complex_Number;
                    zx : in out Matrix; row,col,incx : in integer32 );

  -- DESCRIPTION :
  --   Scales the vector starting at zx(row,col) with the constant za.

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Vector; ind,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 );

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at ind and in y
  --   at (rwy,cly), using increments incx and incy to advance.

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Vector; ind,incy : in integer32 );

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x at (rwx,clx)
  --   and in y at ind, using increments incx and incy to advance.

  procedure zaxpy ( n : in integer32; z : in Complex_Number;
                    x : in Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 );

  -- DESCRIPTION :
  --   Add to y the vector x times z, starting in x and y at the
  --   respective (rows,columns): (rwx,clx) and (rwy,cly), using
  --   increments incx and incy to advance in the rows.

  function zdotc ( n : in integer32; x : in Matrix;
                   rwx,clx,incx : in integer32;
                   y : in Matrix; rwy,cly,incy : in integer32 )
                 return Complex_Number;

  -- DESCRIPTION :
  --   Returns the dot product of two vectors in two columns of matrices
  --   x and y, starting at rows rwx and rwy respectively, with respective
  --   increments in incx and incy.

  procedure drotg ( da,db,c,s : in out double_float );

  -- DESCRIPTION :
  --   Constructs Givens plane rotation.

  procedure zdrot ( n : in integer32;
                    x : in out Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32;
                    c,s : in double_float );

  -- DESCRIPTION :
  --   Applies a plane rotation where the cos and sin are c and s
  --   and the vectors are in the columns of the matrices x and y,
  --   starting at (rwx,clx) and (rwy,cly) advancing in the rows with
  --   increments incx and incy respectively.

  procedure zswap ( n : in integer32;
                    x : in out Matrix; rwx,clx,incx : in integer32;
                    y : in out Matrix; rwy,cly,incy : in integer32 );

  -- DESCRIPTION :
  --   Interchanges two vectors in the columns of the matrices x and y,
  --   respectively starting at (rwx,clx) and (rwy,cly), and advancing
  --   in the rows with the respective increments incx and incy.

end Standard_Complex_BLAS_Helpers;
