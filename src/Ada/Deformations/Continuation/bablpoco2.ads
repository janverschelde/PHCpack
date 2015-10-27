procedure bablpoco2 ( targetname,startname,outfilename : in string );

-- DESCRIPTION :
--   This is the polynomial continuation routine, to run in batch
--   processing mode or as a black box routine, in double double precision.
--   If the names of the file do not lead to the proper data,
--   then the user is prompted to provide file names.

-- ON INPUT :
--   targetname      name of the file where the target system is;
--   startname       name of the file where the start system is;
--   outfilename     name of the output file.
