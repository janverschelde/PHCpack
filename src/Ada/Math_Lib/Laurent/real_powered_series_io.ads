with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Real_Powered_Series_IO is

-- DESCRIPTION :
--   The input and output formats are defined for real powered series,
--   defined by a vector of coefficients of range 0..d,
--   and a vector of powers of range 1..d, where d is a positive integer.

  function to_string ( c : Standard_Complex_Vectors.Vector;
                       p : Standard_Floating_Vectors.Vector;
                       t : character := 't' ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of a series given by
  --   coefficients c and powers p, using the symbol t as the variable.

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

end Real_Powered_Series_IO;
