with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with String_Splitters;                   use String_Splitters;

package String_System_Readers is

-- DESCRIPTION :
--   Provides a Wrapper around a routine to read a system from file,
--   as a pointer to an array of strings, needed in various drivers.
--   Reading a system as an array of strings allows to defer the
--   numerical evaluation of rational constants, e.g. 2/3,
--   until the working precision is set.

  procedure Read_System
             ( file : in out file_type; name : in string;
               n,m : out natural32; p : out Link_to_Array_of_Strings );

  -- DESCRIPTION :
  --   Tries to open the file with given name in the string filename
  --   with intent to read a polynomial system.  If this fails, then
  --   a message is written to the user and p will be null on return.
  --   Otherwise, the p on return will contain a polynomial system
  --   and the file on return is open for input.

  -- ON ENTRY :
  --   name     name of the input file, could be an empty string.
  --
  -- ON RETURN :
  --   file     input file opened for reading, could be a file
  --            with the given name, or a name entered by the user;
  --   n        number of equations;
  --   m        number of variables;
  --   p        array of n strings representing polynomials in m variables.

end String_System_Readers;
