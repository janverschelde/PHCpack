with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;

package Communications_with_User is

-- DESCRIPTION :
--   The routines of this package are used to
--    1) determine the names of input/output files:
--       + create a file: prevent over-writing of important files;
--       + open a file: check first its existence.
--    2) make selections from menu's: choose one of the alternatives.

  generic
    with function Valid_Alternative ( alt : character ) return boolean;
  procedure Ask ( ans : out character );

  -- DESCRIPTION :
  --   This procedure keeps reading a character from standard input,
  --   until a valid one has been given.
  --   The function Valid_Alternative decides whether a certain
  --   character corresponds to a valid answer to the question.

  procedure Ask_Yes_or_No ( ans : out character );

  -- DESCRIPTION :
  --   Keeps reading a character from standard output,
  --   until the user gives a 'y' or a 'n'.

  procedure Ask_Alternative ( ans : out character; alternatives : in string );

  -- DESCRIPTION :
  --   This procedure keeps reading a character from standard input,
  --   until a character that belongs to the string s has been given.

  procedure Ask_Alternative
                ( ans : in out string; alternatives : string;
                  prefix : in character );

  -- DESCRIPTION :
  --   Ask the user to give a character that occurs in the string of
  --   alternatives, eventually preceded by the given prefix character.
  --   This procedure keeps reading till a valid choice has been made.

  -- REQUIRED : ans'range = 1..2.

  procedure Read_Name_and_Open_File ( file : in out file_type );

  -- DESCRIPTION :
  --   This procedure prompts the user the enter a file name
  --   and tries to open the file in the input mode.
  --   If case of file does not exist, another file name will be asked.

  procedure Read_Name_and_Open_File
              ( file : in out file_type; name : out Link_to_String );

  -- DESCRIPTION :
  --   Does the same as the other Read_Name_and_Open_File, and returns
  --   in addition the name of the file that was finally opened.

  procedure Read_Name_and_Create_File ( file : in out file_type );

  -- DESCRIPTION :
  --   This procedure prompts the user the enter a file name
  --   and tries to create the file.  In case file already exists,
  --   permission will be asked to destroy the existing file.

  procedure Read_Name_and_Create_File
              ( file : in out file_type; name : out Link_to_String );

  -- DESCRIPTION :
  --   Does the same as the other Read_Name_and_Create_File, and
  --   returns in addition the name of the file that was created.

  procedure Read_Name_and_Append_File ( file : in out file_type );

  -- DESCRIPTION :
  --   This procedure prompts the user the enter a file name
  --   and tries to open the file in the append mode.
  --   If case of failure, another file name will be asked.

  procedure Open_Input_File
               ( file : in out file_type; filename : in string );
  procedure Open_Input_File
               ( file : in out file_type; filename : in string;
                 name : out Link_to_String );

  -- DESCRIPTION :
  --   Tries to open a file for input, starting with the given file name.
  --   If the opening of the file with the given name fails,
  --   then the procedure `Read_Name_and_Open_File' will be invoked.
  --   When name is supplied as parameter, then the name the user gave
  --   is returned as well.

  procedure Create_Output_File
               ( file : in out file_type; filename : in string );
  procedure Create_Output_File
               ( file : in out file_type; filename : in string;
                 name : out Link_to_String );

  -- DESCRIPTION :
  --   This procedure creates an output file, starting with the given
  --   file name.  If the creation of a file with this name fails,
  --   then the procedure `Read_Name_and_Create_File' will be invoked.
  --   If desired, the name the user gives is returned as well.

  procedure Open_Append_File
               ( file : in out file_type; filename : in string );

  -- DESCRIPTON :
  --   This procedure opens a file in the append mode, starting with
  --   the given name.  If case of failure, the interactive procedure
  --   `Read_Name_and_Append_File' will be invoked.

  function Prompt_for_Precision return character;

  -- DESCRIPTION :
  --   Prompts the user for the working precision and returns '0'
  --   for standard double, '1' for double double, or '2' for quad double.

  procedure End_of_Input_Message;

  -- DESCRIPTON :
  --   After a new_line, writes
  --   "No more input expected.  See output file for results."
  --   and then executes another new_line.

end Communications_with_User;
