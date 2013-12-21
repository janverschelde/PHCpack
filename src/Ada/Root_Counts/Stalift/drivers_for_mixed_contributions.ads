with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;

package Drivers_for_Mixed_Contributions is

-- DESCRIPTION :
--   This package provides some drivers to compute essential sets of
--   a tuple of support sets.

-- SWEEPING (once/full) WITH (simple/exhaustive) CRITERION :

  procedure Once_Simple_Sweep 
                ( file : in file_type; L : in out Array_of_Lists;
                  nred : out natural32 );

  procedure Once_Exhaustive_Sweep 
                ( file : in file_type; L : in out Array_of_Lists;
                  nred : out natural32 );

  procedure Full_Simple_Sweep
                ( file : in file_type; L : in out Array_of_Lists;
                  nred : out natural32 );

  procedure Full_Exhaustive_Sweep
                ( file : in file_type; L : in out Array_of_Lists;
                  nred : out natural32 );

  -- DESCRIPTION :
  --   Applies simple/exhaustive criterion once/full to the point lists.
  --   Full means that the lists are scanned until the criterion fails
  --   for all points, while with the prefix Once, only one sweep is performed.

  -- ON RETURN :
  --   l          reduced set of supports;
  --   nred       number of eliminated vectors with zero contribution.

end Drivers_for_Mixed_Contributions;
