with Standard_Natural_Numbers;          use Standard_Natural_Numbers;

procedure babldmvc ( nt : in natural32; infilename,outfilename : in string );

-- DESCRIPTION :
--   This is the mixed volume computation procedure, to run in batch
--   processing mode or as a black box routine.

-- ON ENTRY :
--   nt            the number of tasks for the polyhedral path trackers,
--                 set nt to zero for the sequential code;
--   infilename    the name of the input file;
--   outfilename   the name of the output file.
