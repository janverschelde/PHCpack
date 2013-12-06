with text_io;                            use text_io;

package File_Scanning is

-- DESCRIPTION :
--   This package collects some scanning operations on text files.

  function Scan_Line_for_Number ( file : in file_type ) return natural;

  -- DESCRIPTION :
  --   Until the end of line is reached, a symbol is read and
  --   the corresponding number is returned.
  --   If no number on the line is present, zero is returned.

  procedure Scan ( file : in file_type; ch : in character;
                   found : out boolean );

  -- DESCRIPTION :
  --   Scans the file on the search of a character.
  --   When the character has been found, then the procedure stops
  --   and sets found to true, otherwise, found will be false and
  --   End_of_File(file) will be true.

  procedure Scan ( file : in file_type; banner : in string;
                   found : out boolean );
  procedure Scan_Line ( file : in file_type; banner : in string;
                        found : out boolean );
  -- DESCRIPTION :
  --   Scans the file on the search of a text banner.
  --   The Scan_Line limits the search to the current line.

  -- ON INPUT :
  --   file        a file opened for input;
  --   banner      a string.

  -- ON RETURN :
  --   found       true if the banner has been found, false otherwise.

  procedure Scan_and_Skip ( file : in file_type; banner : in string;
                            found : out boolean );

  -- DESCRIPTION :
  --   The line where the banner has been found will be skipped.

end File_Scanning;
