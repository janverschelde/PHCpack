with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

procedure mainfeed ( infilename,outfilename : in string;
                     verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the main driver to the algorithms written by Yusong Wang,
--   to realize the feedback laws computed by the Pieri homotopies.
--   The first two arguments are the names of the input and output files.
--   The last argument is the verbose level.
