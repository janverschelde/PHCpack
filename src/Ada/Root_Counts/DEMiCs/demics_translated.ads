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

  function Mixed_Labels
             ( p : Poly_Sys; monitor : boolean := true;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the mixed volume and the labels to the mixed cell indices.
  --   If monitor is true, then the cell is written to screen each time
  --   the cell is added to the data stored in demics_output_cells.

end DEMiCs_Translated;
