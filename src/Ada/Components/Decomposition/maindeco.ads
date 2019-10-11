with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

procedure maindeco ( nt : in natural32; infilename,outfilename : in string;
                     verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This routine is called by the central dispatcher of phc,
--   to perform a numerical irreducible decomposition.

-- ON ENTRY :
--   nt             the number of tasks, if 0 then no multitasking,
--                  otherwise nt tasks will be used to track the paths;
--   infilename     file name for input (a polynomial system),
--                  if empty, then the user will be prompted to supply
--                  the name of an input file;
--   outfilename    name of file for output, if empty, then the user will
--                  be asked to supply the name of an output file;
--   verbose        the verbose level.
