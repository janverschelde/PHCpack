with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
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

package Random_Coefficient_Systems is

-- DESCRIPTION :
--   This package offers functions to create polynomials and systems
--   with given support but with random coefficients.

  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return Standard_Complex_Polynomials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return DoblDobl_Complex_Polynomials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return QuadDobl_Complex_Polynomials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return Standard_Complex_Polynomials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return DoblDobl_Complex_Polynomials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return QuadDobl_Complex_Polynomials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return Standard_Complex_Laurentials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return DoblDobl_Complex_Laurentials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return QuadDobl_Complex_Laurentials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return Standard_Complex_Laurentials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return DoblDobl_Complex_Laurentials.Poly;
  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns a polynomial in n variables with random coefficients
  --   whose support equals the first n-vectors of L.

  -- REQUIRED :
  --   The range of the vectors in L must allow a range of 1..n,
  --   but may be longer than n.

  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns a system of L'range with random polynomials using
  --   the first n elements in the lists L.

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns a system of n random polynomials using the first n elements
  --   in the lists L, according to the type of mixture in mix.

  procedure Add_Constant ( p : in out Standard_Complex_Polynomials.Poly );
  procedure Add_Constant ( p : in out DoblDobl_Complex_Polynomials.Poly );
  procedure Add_Constant ( p : in out QuadDobl_Complex_Polynomials.Poly );
  procedure Add_Constant ( p : in out Standard_Complex_Laurentials.Poly );
  procedure Add_Constant ( p : in out DoblDobl_Complex_Laurentials.Poly );
  procedure Add_Constant ( p : in out QuadDobl_Complex_Laurentials.Poly );
  procedure Add_Constant ( p : in out Standard_Complex_Polynomials.Poly;
                           c : in Standard_Complex_Numbers.Complex_Number );
  procedure Add_Constant ( p : in out DoblDobl_Complex_Polynomials.Poly;
                           c : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Add_Constant ( p : in out QuadDobl_Complex_Polynomials.Poly;
                           c : in QuadDobl_Complex_Numbers.Complex_Number );
  procedure Add_Constant ( p : in out Standard_Complex_Laurentials.Poly;
                           c : in Standard_Complex_Numbers.Complex_Number );
  procedure Add_Constant ( p : in out DoblDobl_Complex_Laurentials.Poly;
                           c : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Add_Constant ( p : in out QuadDobl_Complex_Laurentials.Poly;
                           c : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   If the polynomial p does not yet already contain a constant,
  --   then the number c will be added.
  --   When c is omitted, a random complex number will be added.

  procedure Add_Constants ( p : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Add_Constants ( p : in out Standard_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Adds random constants to the polynomials in p in case the
  --   polynomials have not yet a constant coefficient in them.

end Random_Coefficient_Systems;
