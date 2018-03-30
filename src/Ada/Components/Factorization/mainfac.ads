with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

procedure mainfac ( nt : in natural32; infilename,outfilename : in string );

-- DESCRIPTION :
--   This routine is called by the central dispatcher of phc,
--   to perform a factorization of a pure dimensional solution
--   set into irreducible components. 

-- ON ENTRY :
--   nt             the number of tasks, if 0 then no multitasking,
--                  otherwise nt tasks will be used wherever defined;
--   infilename     file name for input (embedded system+generic points),
--                  if empty, then the user will be prompted to supply
--                  the name of an input file;
--   outfilename    name of file for output, if empty, then the user will
--                  be asked to supply the name of an output file.
