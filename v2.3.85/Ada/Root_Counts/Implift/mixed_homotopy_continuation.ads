with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Trees_of_Vectors;                   use Trees_of_Vectors;

package Mixed_Homotopy_Continuation is

-- DESCRIPTION :
--   This package provides solvers based on the principle
--   of mixed continuation.

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    bkk : out natural32; sols : in out Solution_List );

  -- ON ENTRY :
  --   file         a file for writing intermediate results;
  --   p            a Laurent polynomial system.

  -- ON RETURN :
  --   bkk          the BKK bound of p;
  --   sols         the computed solutions.

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    tv : in Tree_of_Vectors; bkk : out natural32;
                    sols : in out Solution_List );

  -- ON ENTRY :
  --   file         a file for writing intermediate results;
  --   p            a Laurent polynomial system;
  --   tv           the tree of vectors containing the useful directions.


  -- ON RETURN :
  --   bkk          the BKK bound of p;
  --   sols         the computed solutions.

end Mixed_Homotopy_Continuation;
