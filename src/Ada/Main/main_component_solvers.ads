with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Main_Component_Solvers is

-- DESCRIPTION :
--   A numerical irreducible decomposition in blackbox mode (phc -B)
--   is available in double, double double, and quad double precision.

  procedure Standard_Main 
              ( nt : in natural32; infilename,outfilename : in string;
                verbose : in integer32 := 0 );
  procedure DoblDobl_Main 
              ( nt : in natural32; infilename,outfilename : in string;
                verbose : in integer32 := 0 );
  procedure QuadDobl_Main 
              ( nt : in natural32; infilename,outfilename : in string;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines the phc -B, to compute a numerical irreducible decomposition
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   verbose        the verbose level.

end Main_Component_Solvers;
