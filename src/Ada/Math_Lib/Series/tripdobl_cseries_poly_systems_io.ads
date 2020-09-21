with text_io;                           use text_io;
with TripDobl_CSeries_Poly_Systems;     use TripDobl_CSeries_Poly_Systems;

package TripDobl_CSeries_Poly_Systems_io is

-- DESCRIPTION :
--   Provides basic output procedures for polynomial systems,
--   with series coefficients in triple double precision.

  procedure put ( p : in Poly_Sys );
  procedure put ( file : in file_type; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Write the polynomial p to standard output or to file.

end TripDobl_CSeries_Poly_Systems_io;
