procedure mainsolve ( infilename,outfilename : in string );

-- DESCRIPTION :
--   This routine is called by the central dispatcher of phc,
--   to apply the equation-by-equation solver.

-- ON ENTRY :
--   infilename     file name for input for a polynomial system,
--                  if empty, then the user will be prompted to supply
--                  the name of an input file;
--   outfilename    name of file for output, if empty, then the user will
--                  be asked to supply the name of an output file.
