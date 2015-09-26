with Standard_Natural_Numbers;          use Standard_Natural_Numbers;

procedure bablphc2 ( nt : in natural32; infilename,outfilename : in string );

-- DESCRIPTION :
--   This is the main interactive driver for the homotopy continuation
--   package for the blackbox solution of polynomial systems,
--   with computations done with double double arithmetic.
--   This routine is executed with the option -b2 of phc.

-- ON ENTRY :
--   nt             the number of tasks, if 0 then no multitasking,
--                  otherwise nt tasks will be used to track the paths;
--   infilename     the name of the input file;
--   outfilename    the name of the output file.
