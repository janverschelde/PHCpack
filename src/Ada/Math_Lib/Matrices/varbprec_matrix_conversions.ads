with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Double_Double_Matrices;
with Quad_Double_Matrices;
with Multprec_Floating_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;

package VarbPrec_Matrix_Conversions is

-- DESCRIPTION :
--   Often we want to convert matrices of various precisions.
--   This package collects routines to convert between matrices of
--   different types of precision for use in variable precision solvers,
--   for real and complex numbers.

  function d2dd ( mtx : Standard_Floating_Matrices.Matrix )
                return Double_Double_Matrices.Matrix;
  function d2dd ( mtx : Standard_Complex_Matrices.Matrix )
                return DoblDobl_Complex_Matrices.Matrix;
  function d2qd ( mtx : Standard_Floating_Matrices.Matrix )
                return Quad_Double_Matrices.Matrix;
  function d2qd ( mtx : Standard_Complex_Matrices.Matrix )
                return QuadDobl_Complex_Matrices.Matrix;
  function d2mp ( mtx : Standard_Floating_Matrices.Matrix )
                return Multprec_Floating_Matrices.Matrix;
  function d2mp ( mtx : Standard_Complex_Matrices.Matrix )
                return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Converts a floating-point matrix in standard double precision
  --   to a matrix in double double (dd), quad double (qd) precision,
  --   or arbitrary multiprecision (mp).

  function dd2d ( mtx : Double_Double_Matrices.Matrix )
                return Standard_Floating_Matrices.Matrix;
  function dd2d ( mtx : DoblDobl_Complex_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix;
  function dd2qd ( mtx : Double_Double_Matrices.Matrix )
                 return Quad_Double_Matrices.Matrix;
  function dd2qd ( mtx : DoblDobl_Complex_Matrices.Matrix )
                 return QuadDobl_Complex_Matrices.Matrix;
  function dd2mp ( mtx : Double_Double_Matrices.Matrix )
                 return Multprec_Floating_Matrices.Matrix;
  function dd2mp ( mtx : DoblDobl_Complex_Matrices.Matrix )
                 return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Converts a matrix in double double precision to a matrix
  --   in standard double (d), quad double (qd) precision,
  --   or arbitrary multiprecision (mp).

  function qd2d ( mtx : Quad_Double_Matrices.Matrix )
                return Standard_Floating_Matrices.Matrix;
  function qd2d ( mtx : QuadDobl_Complex_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix;
  function qd2dd ( mtx : Quad_Double_Matrices.Matrix )
                 return Double_Double_Matrices.Matrix;
  function qd2dd ( mtx : QuadDobl_Complex_Matrices.Matrix )
                 return DoblDobl_Complex_Matrices.Matrix;
  function qd2mp ( mtx : Quad_Double_Matrices.Matrix )
                 return Multprec_Floating_Matrices.Matrix;
  function qd2mp ( mtx : QuadDobl_Complex_Matrices.Matrix )
                 return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Converts a matrix in quad double precision to a matrix
  --   in standard double (d) or double double (dd) precision,
  --   or arbitrary multiprecision (mp).

  procedure Set_Size ( mtx : in out Multprec_Floating_Matrices.Matrix;
                       size : in natural32 );
  procedure Set_Size ( mtx : in out Multprec_Complex_Matrices.Matrix;
                       size : in natural32 );

  -- DESCRIPTION :
  --   Sets the size of the matrix mtx to the given value of size.

end VarbPrec_Matrix_Conversions;
