with text_io;                            use text_io;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;      use QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_Tracked_Solutions_io is

-- DESCRIPTION :
--   The output file of a path tracker contains the target system,
--   start system, the start solutions, and the solutions of the target.
--   The procedures in this package scan the output file for those data.
--   Coefficients are read in quad double precision.

  procedure get ( file : in file_type; 
                  lp,lq : out Link_to_Poly_Sys;
                  psols,qsols : out Solution_List;
                  verbose : in boolean := false );
  procedure get ( file : in file_type; 
                  lp,lq : out Link_to_Laur_Sys;
                  psols,qsols : out Solution_List;
                  verbose : in boolean := false );

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.

  -- ON ENTRY :
  --   file       file must be opened for reading;
  --   verbose    flag for optional output during reading.

  -- ON RETURN :
  --   lp         target system;
  --   lq         start system;
  --   psols      the solutions on file of the target system;
  --   qsols      start solutions.

end QuadDobl_Tracked_Solutions_io;
