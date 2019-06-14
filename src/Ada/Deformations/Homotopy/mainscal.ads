with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

procedure mainscal ( infilename,outfilename : in string;
                     verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the routine for scaling a polynomial system or a solution list,
--   as called by the option handlers of the phc executable.
--   The first two arguments are the names of the input and output files.
--   The last argument is the verbose level.
