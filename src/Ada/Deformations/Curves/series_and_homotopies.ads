with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Systems;

package Series_and_Homotopies is

-- DESCRIPTION :
--   A homotopy in one parameter is naturally encoded as a polynomial
--   system with truncated power series as coefficients.

  function Create ( h : in Standard_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return Standard_Series_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Given in h the output of Standard_Homotopy.Homotopy_System
  --   and in idx the index to the homotopy continuation parameter,
  --   on return is polynomial system with coefficients power series
  --   in the homotopy continuation parameter.

  function Eval ( p : Standard_Series_Polynomials.Poly;
                  t : double_float )
                return Standard_Complex_Polynomials.Poly;
  function Eval ( p : Standard_Series_Polynomials.Poly;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Every series coefficient in p is evaluated at t.
  --   On return is a polynomial with complex coefficients
  --   in standard double precision.

  function Eval ( h : Standard_Series_Poly_Systems.Poly_Sys;
                  t : double_float )
                return Standard_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : Standard_Series_Poly_Systems.Poly_Sys;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Every series coefficient in h is evaluated at t.
  --   On return is the evaluated polynomial system,
  --   with complex coefficients in standard double precision.

end Series_and_Homotopies;
