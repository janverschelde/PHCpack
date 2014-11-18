with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;

package Multprec_Complex_Numbers_io is

-- DESCRIPTION :
--   This package provides input/output routines for complex numbers.
--   A complex number is displayed as two floating numbers, representing
--   respectively the real and imaginary part of the complex number.

-- INPUT :

  procedure get ( x : in out Complex_Number );
  procedure get ( file : in file_type; x : in out Complex_Number );

  -- DESCRIPTION :
  --   Reads a complex number x from standard input or from file.

  procedure get ( s : in string; c : in out Complex_Number;
                  last : out integer );

  -- DESCRIPTION :
  --   Parses the string s for a complex number c.
  --   On return in last is that first nondigit character that follows c
  --   in the string s.

-- OUTPUT :

  procedure put ( x : in Complex_Number );
  procedure put ( file : in file_type; x : in Complex_Number );

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 );

  procedure put ( c : in Complex_Number; dp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; dp : in natural32 ); 

  function Character_Size ( c : Complex_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters the string representation
  --   of the complex floating-point number c as the sum of the
  --   character sizes of the real and imaginary parts of c,
  --   plus two for the two separating spaces in between.

  procedure put ( s : out string; c : in Complex_Number );

  -- DESCRIPTION :
  --   Writes the number f to a string.

  -- REQUIRED :
  --   s'last = Character_Size(c).

end Multprec_Complex_Numbers_io;
