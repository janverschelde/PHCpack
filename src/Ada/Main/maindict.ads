with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

procedure maindict ( infilename,outfilename : in string;
                     verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This routine converts output solutions into a Python dictionary
--   as called by the central dispatcher.
--   The first two arguments are the names of the input and output file.
--   The last argument is the verbose level.

-- FUNCTIONALITY :
--   The program strips the output file of phc (with name infilename)
--   of all intermediate output, scans that file till "THE SOLUTIONS"
--   is reached and then converts the solutions into a Python dictionary.
--   If the outfilename is the empty string, then the output will
--   be written to the screen, otherwise to the file name.
