with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

procedure bablsolve ( p : in Poly_Sys; outname : in string;
                      verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This routine is called by the blackbox solver of PHCpack,
--   to apply the progressive equation-by-equation solver,
--   in case the system has more than one equation and a
--   number of variables different from the number of equations.

-- ON ENTRY :
--   p              system of polynomial equations;
--   outname        the name of the output file (may be empty);
--   verbose        the verbose level.
