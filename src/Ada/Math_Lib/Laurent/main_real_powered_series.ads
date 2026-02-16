with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;

package Main_Real_Powered_Series is

-- DESCRIPTION :
--   Defines the access to computing real powered series
--   from the main procedure phc.

  function Is_Linear ( p : Standard_Complex_Laurentials.Poly )
                     return boolean;
  function Is_Linear ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                     return boolean;

  -- DESCRIPTION :
  --   Returns true if the monomials in p define a linear system.
  --   Returns false otherwise.

  procedure Main ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Main procedure to compute real powered series.
  --   The verbose level is given by vrblvl.

end Main_Real_Powered_Series;
