with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
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

end Series_and_Homotopies;
