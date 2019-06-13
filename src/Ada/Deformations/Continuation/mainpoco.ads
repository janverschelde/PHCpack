with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

procedure mainpoco ( nt : in natural32; infilename,outfilename : in string;
                     prclvl : in natural32; verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the polynomial continuation routine.

-- ON ENTRY :
--   nt             the number of tasks, if 0 then no multitasking,
--                  otherwise nt tasks will be used to track the paths;
--   infilename     the name of the input file;
--   outfilename    the name of the output file.
--   prclvl         level of precision, passed as the number following
--                  phc -p at the command line, either 1, 2, or 4
--                  for double, double double, or quad double precision;
--   verbose        the verbose level.
