with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

procedure mainred ( infilename,outfilename : in string;
                    verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the routine for reducing a polynomial system,
--   as called by the central dispatcher.
--   The first two arguments are the names of the input and output files.
--   The last argument is the verbose level.
