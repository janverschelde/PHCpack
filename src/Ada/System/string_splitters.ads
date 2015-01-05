with text_io;                            use text_io;

package String_Splitters is

-- DESCRIPTION :
--   The routines in this package allow to read a polynomial system
--   into an array of strings, all terminated by the semicolon.

-- DATA STRUCTURES :

  type Link_to_String is access string;
  type Array_of_Strings is array ( integer range <> ) of Link_to_String;
  type Link_to_Array_of_Strings is access Array_of_Strings;

-- READ OPERATORS :

  function Read_String return string;

  -- DESCRIPTION :
  --   Reads a string from standard input and returns it to the caller.
  --   The length of the input string must be smaller than 256 characters.

  procedure Append ( s : in out Link_to_String; t : in string );

  -- DESCRIPTION :
  --   Appends the string t to the string s.

  function Read_till_Delimiter
             ( file : in file_type; d : in character )
             return Link_to_String;

  -- DESCRIPTION :
  --   Reads the file till the end or till d is encountered.
  --   The resulting string is returned.

  function Read_till_Delimiter
             ( file : in file_type; n : in natural; d : in character )
             return Array_of_Strings;

  -- DESCRIPTION :
  --   Reads n strings from file till end of file or d returned for 
  --   each of them.

  procedure get ( n,m : out natural; p : out Link_to_Array_of_Strings );
  procedure get ( file : in file_type; 
                  n,m : out natural; p : out Link_to_Array_of_Strings );

  -- DESCRIPTION :
  --   Reads a polynomial system as an array of strings, returned in p.

  -- ON RETURN :
  --   n         number of polynomials in p;
  --   m         number of variables in each polynomial;
  --   p         array of strings, each delimited by semicolon.

  function Count_Delimiters ( s : string; d : character ) return natural;

  -- DESCRIPTION :
  --   Returns the number of occurrences of the character d in the string s.

  function Split ( n : natural; s : string; d : character )
                 return Array_of_Strings;

  -- DESCRIPTION :
  --   Returns an array of n strings by splitting s
  --   along the delimiter d.

  -- ON ENTRY :
  --   n         number of occurrences of d in s;
  --   s         string representation of n polynomials, separated by d;
  --   d         delimiter symbol is typically a semicolon.

-- DESTRUCTORS :

  procedure Clear ( s : in out Link_to_String );
  procedure Clear ( s : in out Array_of_Strings );
  procedure Clear ( s : in out Link_to_Array_of_Strings );

  -- DESCRIPTION :
  --   Deallocation of the memory space occupied by s.

end String_Splitters;
