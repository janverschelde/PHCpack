with text_io;                            use text_io;

procedure valipoco ( pocofile,resultfile : in file_type );

-- DESCRIPTION :
--   This is an intermediate stage in the polyhedral end game.
--   Scans the output file of poco for the computed directions,
--   residuals and estimated multiplicities of the solutions.
--   Computes the frequency table of path directions.
