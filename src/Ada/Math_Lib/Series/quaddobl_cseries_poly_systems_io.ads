with text_io;                           use text_io;
with QuadDobl_CSeries_Poly_Systems;     use QuadDobl_CSeries_Poly_Systems;

package QuadDobl_CSeries_Poly_Systems_io is

-- DESCRIPTION :
--   Provides basic output procedures for polynomial systems,
--   with series coefficients in double precision.

  procedure put ( p : in Poly_Sys );
  procedure put ( file : in file_type; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Write the polynomial p to standard output or to file.

end QuadDobl_CSeries_Poly_Systems_io;
