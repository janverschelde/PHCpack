with text_io;                            use text_io;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Parse_Strings_to_Polynomials is

-- DESCRIPTION :
--   Below is a collection of procedures to check whether a collection
--   of strings (often read from file) represents a valid input to phc.
--   The procedures have been developed for use in phc -g.
--   The parsing of strings into polynomials happens in two stages:
--   (1) Reading of the first line with the number of polynomials
--       and the number of variables, followed by the reading of
--       as many strings as the first number on the first line.
--       Each polynomial must end with a semicolon.
--   (2) Parsing from the strings into valid polynomials.

  procedure First_Line_Format ( file : in file_type );

  -- DESCRIPTION :
  --   Writes the format for the first line to file.
  --   Called when something already goes wrong on the first line.

  procedure Parse_Dimensions
              ( outfile : in file_type; to_file : in boolean;
                s : in string; nq,nv : out natural32;
                fail : out boolean );

  -- DESCRIPTION :
  --   Parses the number of equations and variable from the first
  --   line given on input in the string s.

  -- ON ENTRY :
  --   outfile  must be opened for output if to_file is true;
  --   to_file  if true, then all messages are written to outfile,
  --            if false, then all error message are written to screen
  --            and if to_file is off more diagnostics are written,
  --            also when all parsing goes well;
  --   s        first line on the input.

  -- ON RETURN :
  --   nq       the first number on the line is the number of equations;
  --   nv       the second optional number on the line is the number
  --            of variables, if nv is omitted, then nv = nq.
  --   fail     the first line is expected to contain either one
  --            or two positive natural numbers, for any other input,
  --            fail will be printed.

  procedure Read_First_Line
              ( infile,outfile : in file_type; to_file : in boolean;
                nq,nv : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Reads the first line from file and tries to extract
  --   the number of equations and variables.

  -- ON ENTRY :
  --   infile   input file which must be opened for input;
  --   outfile  output file which must be opened for output,
  --            if the to_file variable is set to true;
  --   to_file  if true, then error messages are written to file,
  --            if false, then error messages will appear on screen.
  
  -- ON RETURN :
  --   nq       the number of equations;
  --   nv       the number of variables;
  --   fail     true if the format of the first line was wrong.

  procedure Read_Polynomials
             ( infile,outfile : in file_type; to_file : in boolean;
               n : in natural32;
               la : out String_Splitters.Link_to_Array_of_Strings;
               fail : out boolean );

  -- DESCRIPTION :
  --   Reads an array of n strings from file.

  -- ON ENTRY :
  --   infile  input file must be opened for input;
  --   outfile output file must be opened for output if to_file is true,
  --           otherwise outfile does not matter;
  --   to_file if true, then all error messages are written to file,
  --           otherwise, all error messages and some extra information
  --           is written to screen;
  --   n       number of polynomials.

  -- ON RETURN :
  --   la      if not fail, la contains n strings;
  --   fail    is true if there no n semicolons on file.

  procedure Parse_Polynomials
              ( file : in file_type; to_file : in boolean;
                n : in natural32;
                a : in String_Splitters.Link_to_Array_of_Strings;
                p : in Link_to_Laur_Sys; fail : out integer32 );

  -- DESCRIPTION :
  --   Parses the strings given in a into polynomials.

  -- REQUIRED : a'range = p'range.

  -- ON ENTRY :
  --   file     to write error messages if to_file is true;
  --   to_file  if true, then error messages are written to file;
  --   n        number of variables;
  --   a        input strings are supposed to contain polynomials.

  -- ON RETURN :
  --   p        parsed polynomials if fail differs from zero;
  --   fail     index of the first string in a where the parsing
  --            raised an exception.

  procedure Create_Output_File
              ( file : in out file_type; name : in string;
                no_output_file : out boolean );

  -- DESCRIPTION :
  --   Prompts the user for the name of an output file if name is empty,
  --   otherwise attempts to create the file with the given name.

  -- ON ENTRY :
  --   name     name for the output file.

  -- ON RETURN :
  --   file     opened for output if no_output_file is false;
  --   no_output_file is true if the user does not want output to file,
  --            if false if the file points to a valid output file.

  procedure Write_Results
              ( file : in file_type;
                nq,nv : in natural32;
                a : in String_Splitters.Link_to_Array_of_Strings;
                p : in Link_to_Laur_Sys; fail_index : in integer32 );

  -- DESCRIPTION :
  --   Writes the results of the parsing to file,
  --   along with a time stamp.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nq       number of equations;
  --   nv       number of variables;
  --   a        strings read from file;
  --   p        parsed polynomials;
  --   fail_index is the index of the string in a that led to a parsing error,
  --            if 0, then all strings were parsed correctly.

  procedure Read_from_File ( inname,outname : in string );

  -- DESCRIPTION :
  --   Tries to open the input file with name in the string inname
  --   and checks whether the outname is empty or not.
  --   Then starts the parsing.

  procedure Read_Input_File_Name;

  -- DESCRIPTION :
  --   Prompts the user for an input file name
  --   and tries to open a file with that name.
  --   If the user does not give the name for an input file,
  --   then there is no output file eiter.

  procedure Main ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   This procedure checks whether the polynomial system given on input is
  --   in a valid input format, by parsing it as a Laurent polynomial system.
  --   This procedure gets called with the option '-g'.
  --   Since phc -b adds to the original input file, the phc -g input copy
  --   will create a copy of the system in the file input to the file copy,
  --   along with some parsing diagnostics and a time stamp.
  --   For systems that are parsed well, the output file starts with
  --   the string version submitted to phc (not the parsed polynomials).
  --   This copying and time stamp of phc -g is a good feature
  --   to keep a catalog of all systems submitted to PHCpack.

end Parse_Strings_to_Polynomials;
