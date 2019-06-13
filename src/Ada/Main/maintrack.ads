with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

procedure maintrack ( targetfilename,startfilename,outfilename : in string;
                      verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This is the main interface to the path trackers.
--   The arguments are the respective names of the files
--   for the target system, start system, output file;
--   and at last the verbose level.
