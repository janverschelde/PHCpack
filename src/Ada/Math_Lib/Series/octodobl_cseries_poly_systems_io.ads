with text_io;                           use text_io;
with OctoDobl_CSeries_Poly_Systems;     use OctoDobl_CSeries_Poly_Systems;

package OctoDobl_CSeries_Poly_Systems_io is

-- DESCRIPTION :
--   Provides basic output procedures for polynomial systems,
--   with series coefficients in octo double precision.

  procedure put ( p : in Poly_Sys );
  procedure put ( file : in file_type; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Write the polynomial p to standard output or to file.

end OctoDobl_CSeries_Poly_Systems_io;
