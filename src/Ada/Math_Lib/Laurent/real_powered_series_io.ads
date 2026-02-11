with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Real_Powered_Series_IO is

-- DESCRIPTION :
--   The input and output formats are defined for real powered series,
--   defined by a vector of coefficients of range 0..d,
--   and a vector of powers of range 1..d, where d is a positive integer.

-- WRITE OUTPUT :

  function to_string ( c : Standard_Complex_Vectors.Vector;
                       p : Standard_Floating_Vectors.Vector;
                       t : character := 't';
                       vrblvl : integer32 := 0 ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of a series given by
  --   coefficients c and powers p, using the symbol t as the variable.
  --   The verbose level is vrblvl.

  procedure put ( c : in Standard_Complex_Vectors.Vector;
                  p : in Standard_Floating_Vectors.Vector;
                  t : in character := 't' );
  procedure put ( file : in file_type;
                  c : in Standard_Complex_Vectors.Vector;
                  p : in Standard_Floating_Vectors.Vector;
                  t : in character := 't' );

  -- DESCRIPTION :
  --   Given in c are the coefficients and in p the powers of a series,
  --   writes the series in symbolic format, using t as the variable.

  procedure put_line ( c : in Standard_Complex_Vectors.Vector;
                       p : in Standard_Floating_Vectors.Vector;
                       t : in character := 't' );
  procedure put_line ( file : in file_type;
                       c : in Standard_Complex_Vectors.Vector;
                       p : in Standard_Floating_Vectors.Vector;
                       t : in character := 't' );

  -- DESCRIPTION :
  --   Given in c are the coefficients and in p the powers of a series,
  --   writes the series in symbolic format, using t as the variable.
  --   After each term, a new line symbol is written.

-- PARSE INPUT :

  function Size_Count
             ( s : string; t : character := 't' ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of occurrences of the character t in s,
  --   which equals the number of powers in the string representation
  --   of the series, with symbol t.

  procedure parse_string
              ( s : in string;
                c : out Standard_Complex_Vectors.Link_to_Vector;
                p : out Standard_Floating_Vectors.Link_to_Vector;
                t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Parses the string s for the coefficients c and powers p
  --   of the symbolic format of a series in t.
  --   The verbose level is vrblvl.

  procedure get ( c : out Standard_Complex_Vectors.Vector;
                  p : out Standard_Floating_Vectors.Vector;
                  t : in character := 't'; vrblvl : in integer32 := 0 );
  procedure get ( file : in file_type;
                  c : out Standard_Complex_Vectors.Vector;
                  p : out Standard_Floating_Vectors.Vector;
                  t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads the constant coefficient c(0) and then as many terms
  --   as the size of c and p, of a series with symbol t.
  --   The verbose level is given in vrblvl.

  procedure get ( size : in integer32;
                  c : out Standard_Complex_Vectors.Link_to_Vector;
                  p : out Standard_Floating_Vectors.Link_to_Vector;
                  t : in character := 't'; vrblvl : in integer32 := 0 );
  procedure get ( file : in file_type; size : in integer32;
                  c : out Standard_Complex_Vectors.Link_to_Vector;
                  p : out Standard_Floating_Vectors.Link_to_Vector;
                  t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Allocates space in c and p for a series of the given size,
  --   which is then read from standard input or from file.
  --   The verbose level is given in vrblvl.

end Real_Powered_Series_IO;
