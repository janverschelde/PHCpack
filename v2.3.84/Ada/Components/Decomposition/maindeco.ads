procedure maindeco ( infilename,outfilename : in string );

-- DESCRIPTION :
--   This routine is called by the central dispatcher of phc,
--   to perform a numerical irreducible decomposition.

-- ON ENTRY :
--   infilename     file name for input (a polynomial system),
--                  if empty, then the user will be prompted to supply
--                  the name of an input file;
--   outfilename    name of file for output, if empty, then the user will
--                  be asked to supply the name of an output file.
