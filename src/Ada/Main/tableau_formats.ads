with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Tableau_Formats is

-- DESCRIPTION :
--   Conversion between tableau and symbolic formats of polynomial systems.
--
-- SYMBOLIC FORMAT :
--
--   2
--    x*y**2 - x**2 + 3;
--    x + 2;
--
-- TABLEAU FORMAT :
--
--   2
--    x y
--   3
--    1 2
--    2 0
--    0 0
--   3
--    1 0
--    0 0
--    1.0 0.0
--   -1.0 0.0
--    3.0 0.0
--    1.0 0.0
--    2.0 0.0

  procedure get ( file : in file_type; realcoeff : in boolean;
                  p : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Reads a polynomial system in tableau format from file.
  --   If realcoeff, then the imaginary part must be omitted.

  procedure put ( file : in file_type; realcoeff : in boolean;
                  p : in Poly_Sys );

  -- DESCRIPTION :
  --   Writes a polynomial system in tableau format on file.
  --   If realcoeff, then the imaginary part is not written.

end Tableau_Formats;
