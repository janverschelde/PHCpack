with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;

package File_Management is

-- DESCRIPTION :
--   This package provides produces and functions to manage files,
--   for use in incremental read and write operations.
--   The benefit of this package is that one function may open or create
--   a file, while another function may then use the file, without having
--   to explicitly use a variable that refers to that file.

  type link_to_file_type is access file_type;

  procedure Silent_Open_Input_File;
  procedure Silent_Open_Input_File ( k : in natural32 );
  procedure Silent_Open_Input_File ( filename : in string );
  procedure Silent_Open_Input_File ( k : in natural32; filename : in string );
  procedure Open_Input_File;
  procedure Open_Input_File ( k : in natural32 );

  -- DESCRIPTION :
  --   Reads the name of a file for input of solutions.
  --   The "Silent_" one just asks for a string of characters,
  --   while the other one prints first a message before prompting.
  --   When the filename is provided, no messages are written.
  --   If k is provided, then it is for reading witness set k.

  procedure Create_Output_File;

  -- DESCRIPTION :
  --   Reads the name of a file for output of solutions.

  function Link_to_Input return link_to_file_type;
  function Link_to_Input ( k : natural32 ) return link_to_file_type;

  -- DESCRIPTION :
  --   Returns the reference to the input file.

  function Link_to_Output return link_to_file_type;

  -- DESCRIPTION :
  --   Returns the reference to the output file.

  procedure Reset_Input_File ( k : in natural32 );

  -- DESCRIPTION :
  --   Resets the file for the k-th witness set for reading.
  --   The number k must be 1 or 2 to have an effect.

  procedure Close_Input_File;
  procedure Close_Input_File ( k : in natural32 );
  procedure Close_Output_File;

  -- DESCRIPTION :
  --   Closes file used for input or output.

end File_Management;
