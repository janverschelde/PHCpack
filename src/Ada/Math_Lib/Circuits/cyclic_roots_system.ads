with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;

package Cyclic_Roots_System is

-- DESCRIPTION :
--   The cyclic n-roots system is an important benchmark system.
--   This package defines the polynomials in the system,
--   via their support sets.

  function Number_of_Monomials ( n,i : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of monomials in the i-th polynomial of
  --   the cyclic n-roots problem.

  function Support_of_Cyclic
             ( n,i : integer32 )
             return Standard_Natural_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the support of the i-th polynomial in the
  --   cyclic n-roots system.

  -- REQUIRED : 0 < i < n

  function Supports_of_Cyclic
             ( n : integer32 )
             return Standard_Natural_VecVecs.Array_of_VecVecs;

  -- DESCRIPTION :
  --   Returns the supports of the cyclic n-roots problem.

  function Standard_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the cyclic polynomial defined by the support in s.
  --   as a polynomial with coefficients in standard double precision.

  function DoblDobl_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return DoblDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the cyclic polynomial defined by the support in s.
  --   as a polynomial with coefficients in double double precision.

  function QuadDobl_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the cyclic polynomial defined by the support in s.
  --   as a polynomial with coefficients in quad double precision.

  function Multprec_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the cyclic polynomial defined by the support in s.
  --   as a polynomial with coefficients of arbitrary precision type.

  function Standard_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system defined by the supports in s,
  --   with coefficients in standard double precision.

  function DoblDobl_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system defined by the supports in s,
  --   with coefficients in double double precision.

  function QuadDobl_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system defined by the supports in s,
  --   with coefficients in quad double precision.

  function Multprec_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system defined by the supports in s,
  --   with coefficients of arbitrary precision type.

end Cyclic_Roots_System;
