with text_io;                            use text_io;
with Standard_Complex_Laur_Systems;

procedure Driver_to_Rank_Supports
            ( outfile : in file_type;
              p : in Standard_Complex_Laur_Systems.Laur_Sys );

-- DESCRIPTION :
--   Interactive driver to determine the rank of the supports of 
--   a Laurent system and to transform the system for positive corank.

-- ON ENTRY :
--   outfile  file that must be opened for output;
--   p        a Laurent polynomial system.
