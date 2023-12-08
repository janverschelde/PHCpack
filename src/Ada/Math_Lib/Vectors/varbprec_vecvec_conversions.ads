with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Double_Double_VecVecs;
with Triple_Double_VecVecs;
with Quad_Double_VecVecs;
with Multprec_Floating_VecVecs;
with DoblDobl_Complex_VecVecs;
with TripDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with PentDobl_Complex_VecVecs;
with OctoDobl_Complex_VecVecs;
with DecaDobl_Complex_VecVecs;
with HexaDobl_Complex_VecVecs;
with Multprec_Complex_VecVecs;

package Varbprec_VecVec_Conversions is

-- DESCRIPTION :
--   Often we want to convert vectors of vectors of various precisions.
--   This package collects routines to convert between vectors of vectors of
--   different types of precision for use in variable precision solvers,
--   for real and complex numbers.

  function d2dd ( mtx : Standard_Floating_VecVecs.VecVec )
                return Double_Double_VecVecs.VecVec;
  function d2dd ( mtx : Standard_Complex_VecVecs.VecVec )
                return DoblDobl_Complex_VecVecs.VecVec;
  function d2qd ( mtx : Standard_Floating_VecVecs.VecVec )
                return Quad_Double_VecVecs.VecVec;
  function d2qd ( mtx : Standard_Complex_VecVecs.VecVec )
                return QuadDobl_Complex_VecVecs.VecVec;
  function d2mp ( mtx : Standard_Floating_VecVecs.VecVec )
                return Multprec_Floating_VecVecs.VecVec;
  function d2mp ( mtx : Standard_Complex_VecVecs.VecVec )
                return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Converts a floating-point matrix in standard double precision
  --   to a matrix in double double (dd), quad double (qd) precision,
  --   or arbitrary multiprecision (mp).

  function dd2d ( mtx : Double_Double_VecVecs.VecVec )
                return Standard_Floating_VecVecs.VecVec;
  function dd2d ( mtx : DoblDobl_Complex_VecVecs.VecVec )
                return Standard_Complex_VecVecs.VecVec;
  function dd2qd ( mtx : Double_Double_VecVecs.VecVec )
                 return Quad_Double_VecVecs.VecVec;
  function dd2qd ( mtx : DoblDobl_Complex_VecVecs.VecVec )
                 return QuadDobl_Complex_VecVecs.VecVec;
  function dd2mp ( mtx : Double_Double_VecVecs.VecVec )
                 return Multprec_Floating_VecVecs.VecVec;
  function dd2mp ( mtx : DoblDobl_Complex_VecVecs.VecVec )
                 return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Converts a matrix in double double precision to a matrix
  --   in standard double (d), quad double (qd) precision,
  --   or arbitrary multiprecision (mp).

  function qd2d ( mtx : Quad_Double_VecVecs.VecVec )
                return Standard_Floating_VecVecs.VecVec;
  function qd2d ( mtx : QuadDobl_Complex_VecVecs.VecVec )
                return Standard_Complex_VecVecs.VecVec;
  function qd2dd ( mtx : Quad_Double_VecVecs.VecVec )
                 return Double_Double_VecVecs.VecVec;
  function qd2dd ( mtx : QuadDobl_Complex_VecVecs.VecVec )
                 return DoblDobl_Complex_VecVecs.VecVec;
  function qd2td ( mtx : Quad_Double_VecVecs.VecVec )
                 return Triple_Double_VecVecs.VecVec;
  function qd2td ( mtx : QuadDobl_Complex_VecVecs.VecVec )
                 return TripDobl_Complex_VecVecs.VecVec;
  function qd2mp ( mtx : Quad_Double_VecVecs.VecVec )
                 return Multprec_Floating_VecVecs.VecVec;
  function qd2mp ( mtx : QuadDobl_Complex_VecVecs.VecVec )
                 return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Converts a matrix in quad double precision to a matrix
  --   in standard double (d) precision, double double (dd) precision,
  --   triple double (td) precision, or arbitrary multiprecision (mp).

  function da2d ( v : DecaDobl_Complex_VecVecs.VecVec )
                return Standard_Complex_VecVecs.VecVec;
  function da2dd ( v : DecaDobl_Complex_VecVecs.VecVec )
                 return DoblDobl_Complex_VecVecs.VecVec;
  function da2td ( v : DecaDobl_Complex_VecVecs.VecVec )
                 return TripDobl_Complex_VecVecs.VecVec;
  function da2qd ( v : DecaDobl_Complex_VecVecs.VecVec )
                 return QuadDobl_Complex_VecVecs.VecVec;
  function da2pd ( v : DecaDobl_Complex_VecVecs.VecVec )
                 return PentDobl_Complex_VecVecs.VecVec;
  function da2od ( v : DecaDobl_Complex_VecVecs.VecVec )
                 return OctoDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Converts the coefficients in v from deca double
  --   to a lower precision.

  function hd2d ( v : HexaDobl_Complex_VecVecs.VecVec )
                return Standard_Complex_VecVecs.VecVec;
  function hd2dd ( v : HexaDobl_Complex_VecVecs.VecVec )
                 return DoblDobl_Complex_VecVecs.VecVec;
  function hd2td ( v : HexaDobl_Complex_VecVecs.VecVec )
                 return TripDobl_Complex_VecVecs.VecVec;
  function hd2qd ( v : HexaDobl_Complex_VecVecs.VecVec )
                 return QuadDobl_Complex_VecVecs.VecVec;
  function hd2pd ( v : HexaDobl_Complex_VecVecs.VecVec )
                 return PentDobl_Complex_VecVecs.VecVec;
  function hd2od ( v : HexaDobl_Complex_VecVecs.VecVec )
                 return OctoDobl_Complex_VecVecs.VecVec;
  function hd2da ( v : HexaDobl_Complex_VecVecs.VecVec )
                 return DecaDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Converts the coefficients in v from hexa double
  --   to a lower precision.

  procedure Set_Size ( mtx : in out Multprec_Floating_VecVecs.VecVec;
                       size : in natural32 );
  procedure Set_Size ( mtx : in out Multprec_Complex_VecVecs.VecVec;
                       size : in natural32 );

  -- DESCRIPTION :
  --   Sets the size of the matrix mtx to the given value of size.

end Varbprec_VecVec_Conversions;
