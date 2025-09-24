with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package DEMiCs_Translated is

-- DESCRIPTION :
--   This package provides an interface to the translated DEMiCs,
--   to compute all mixed cells by dynamic enumeration.
--   DEMiCs was developed by Tomohiko Mizutani, Akiko Takeda, and
--   Masakazu Kojima and licensed under GNU GPL Version 2 or higher.
--   The translated DEMiCs is a literal translation of the C++ into Ada.

  function Mixed_Volume
             ( p : Poly_Sys; vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the mixed volume of the polynomials in p.

end DEMiCs_Translated;
