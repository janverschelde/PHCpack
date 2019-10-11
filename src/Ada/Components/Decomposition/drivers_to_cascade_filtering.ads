with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Drivers_to_Cascade_Filtering is

-- DESCRIPTION :
--   This package contains drivers to the cascade homotopies.

  procedure Standard_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );
  procedure DoblDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );
  procedure QuadDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );

  -- DESCRIPTION :
  --   Interactive driver to call the Witness_Generate procedure,
  --   in standard double, double double, or quad double precision.
  --   Names for the input and output file are passed at the command line.

  -- ON ENTRY :
  --   nt       number of tasks, if zero, then no tasking;
  --   inpname  name of the input file;
  --   outname  name of the output file.

  procedure Driver_to_Witness_Generate
              ( nt : in natural32; inpname,outname : in string );

  -- DESCRIPTION :
  --   Prompts the user for the level of the working precision and
  --   then calls the Standard, DoblDobl, or QuadDobl_Witness_Generate.
  --   Names for the input and output file are passed at the command line.

  -- ON ENTRY :
  --   nt       number of tasks, if zero, then no tasking;
  --   inpname  name of the input file;
  --   outname  name of the output file.

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 );
  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 );
  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does the embedding of the top dimension and runs the cascade,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   nt       the number of tasks, if 0 then no multitasking,
  --            otherwise nt tasks will be used to track the paths;
  --   inpname  name of the input file;
  --   outname  name of the output file;
  --   verbose  the verbose level.

  procedure Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for the level of working precision,
  --   does the embedding of the top dimension and runs the cascade.

  -- ON ENTRY :
  --   nt       the number of tasks, if 0 then no multitasking,
  --            otherwise nt tasks will be used to track the paths;
  --   inpname  name of the input file;
  --   outname  name of the output file;
  --   verbose  is the verbose level.

end Drivers_to_Cascade_Filtering;
