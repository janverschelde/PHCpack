with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

procedure compsolve2 ( nt : in natural32; infilename,outfilename : in string;
                       verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the main interactive driver to compute a numerical irreducible
--   decomposition for polynomial systems in the blackbox option,
--   with computations done with double double arithmetic.
--   This routine is executed with the option -B2 of phc.

-- ON ENTRY :
--   nt             the number of tasks, if 0 then no multitasking,
--                  otherwise nt tasks will be used to track the paths;
--   infilename     the name of the input file;
--   outfilename    the name of the output file;
--   verbose        the verbose level.
