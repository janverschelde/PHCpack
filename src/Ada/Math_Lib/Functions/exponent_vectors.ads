with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Exponent_Vectors is

-- DESCRIPTION :
--   This package facilitates the management of exponent vectors.

-- DATA STRUCTURE : array of exponent vectors

  type Exponent_Vectors_Array is
    array ( integer32 range <> ) of Standard_Integer_VecVecs.Link_to_VecVec;
  type Link_to_Exponent_Vectors_Array is access Exponent_Vectors_Array;

-- CREATORS :

  function Create ( p : Standard_Complex_Laurentials.Poly )
                  return Standard_Integer_VecVecs.VecVec;
  function Create ( p : Standard_Complex_Polynomials.Poly )
                  return Standard_Integer_VecVecs.VecVec;
  function Create ( p : DoblDobl_Complex_Laurentials.Poly )
                  return Standard_Integer_VecVecs.VecVec;
  function Create ( p : DoblDobl_Complex_Polynomials.Poly )
                  return Standard_Integer_VecVecs.VecVec;
  function Create ( p : QuadDobl_Complex_Laurentials.Poly )
                  return Standard_Integer_VecVecs.VecVec;
  function Create ( p : QuadDobl_Complex_Polynomials.Poly )
                  return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   The range of the vector on return is 1..Number_of_Terms(p).
  --   This vector contains copies of all exponents of p, as ordered in p.

  function Create ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Exponent_Vectors_Array;
  function Create ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Exponent_Vectors_Array;
  function Create ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return Exponent_Vectors_Array;
  function Create ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return Exponent_Vectors_Array;
  function Create ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return Exponent_Vectors_Array;
  function Create ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return Exponent_Vectors_Array;

-- SELECTOR :

  function Position ( ev : Standard_Integer_VecVecs.VecVec;
                      v : Standard_Integer_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the position of v in the vector ev.
  --   If v does not occur in ev, then ev'last+1 will be returned.

-- EVALUATORS :

  function Eval ( e : Standard_Integer_Vectors.Vector;
                  c : Standard_Complex_Numbers.Complex_Number;
                  x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the term c*x^e.

  function Eval ( ev : Standard_Integer_VecVecs.VecVec;
                  c,x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the polynomial with coefficients in c and exponents in ev.

  function Eval ( ev : Exponent_Vectors_Array;
                  c : Standard_Complex_VecVecs.VecVec;
                  x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the system with coefficients in c and exponents in ev.

-- DESTRUCTORS :

  procedure Clear ( v : in out Exponent_Vectors_Array );
  procedure Clear ( v : in out Link_to_Exponent_Vectors_Array );

  -- DESCRIPTION :
  --   Clears the allocated memory.

end Exponent_Vectors;
