with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

procedure bablroco ( nt : in natural32; infilename,outfilename : in string;
                     verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the routine for counting the roots of a polynomial system,
--   as called by the top level procedure, in its batch or black box version.


-- ON ENTRY :
--   nt             the number of tasks, if 0 then no multitasking,
--                  otherwise nt tasks will be used to track the paths;
--   infilename     the name of the input file;
--   outfilename    the name of the output file;
--   verbose        the verbose level.
