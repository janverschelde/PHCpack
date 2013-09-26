procedure mainzip ( infilename,outfilename : in string );

-- DESCRIPTION :
--   This routine converts output solutions into Maple format,
--   as called by the central dispatcher.
--   The arguments are the names of the input and output file respectively.

-- FUNCTIONALITY :
--   The program strips the output file of phc (with name infilename)
--   of all intermediate output, scans that file till "THE SOLUTIONS"
--   is reached and then converts the solutions into a Maple format.
--   If the outfilename is the empty string, then the output will
--   be written to the screen, otherwise to the file name.
