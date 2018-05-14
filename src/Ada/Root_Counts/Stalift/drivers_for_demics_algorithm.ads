with text_io;                            use text_io;
with Standard_Complex_Laur_Systems;

package Drivers_for_DEMiCs_Algorithm is

-- DESCRIPTION :
--   This package offers an interface to the DEMiCs Algorithm,
--   for dynamic enumeration of all mixed cells.

  procedure DEMiCs_Algorithm_Info;

  -- DESCRIPTION :
  --   Displays information on the DEMiCs Algorithm to screen.

  procedure Driver_for_DEMiCs_Algorithm
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Interactive driver to the MixedVol Algorithm,
  --   as called by phc -m.
  --   All output is written to file.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   p        a polynomial system.

end Drivers_for_DEMiCs_Algorithm;
