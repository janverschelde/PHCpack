procedure mainfac ( infilename,outfilename : in string );

-- DESCRIPTION :
--   This routine is called by the central dispatcher of phc,
--   to perform a factorization of a pure dimensional solution
--   set into irreducible components. 

-- ON ENTRY :
--   infilename     file name for input (embedded system+generic points),
--                  if empty, then the user will be prompted to supply
--                  the name of an input file;
--   outfilename    name of file for output, if empty, then the user will
--                  be asked to supply the name of an output file.
