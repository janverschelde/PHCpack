with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

procedure bablvali ( infilename,outfilename : in string;
                     verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the routine for validating the solutions of a polynomial system,
--   as called by the central dispatcher, in its batch or black box version.
--   The first two arguments are the names of the input and output file.
--   The last argument is the verbose level.
