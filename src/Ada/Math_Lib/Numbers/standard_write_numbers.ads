with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers; 
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_Write_Numbers is

-- DESCRIPTION :
--   This packages provides utilities to write polynomials with
--   complex coefficients in a more compact way in case the coefficients
--   are integers or real numbers.

  function Is_Imag ( c : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the real part of the number c equals 0.0.

  function Is_Real ( c : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the imaginary part of the number c equals 0.0.

  function Is_Integer ( f : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the error from rounding f to an integer is 0.0.

  function Length ( n : integer32 ) return natural32;
  function Length ( n : double_float ) return natural32;
  function Length ( n : Complex_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters needed to write n.
 
  procedure Write_Number ( file : in file_type; i : in integer32;
                           cnt : out natural32 );

  -- DESCRIPTION : 
  --   Writes the integer number with only one blank space before it.
  --   The cnt on return equals the number of characters written.

  procedure Write_Number ( file : in file_type; f : in double_float;
                           cnt : out natural32 );

  -- DESCRIPTION :
  --   If f is an integer, then f is written as an integer to file,
  --   otherwise the default scientific format is used.
  --   The cnt on return equals the number of characters written.

  procedure Write_Number ( file : in file_type; c : in Complex_Number;
                           cnt : out natural32 );

  -- DESCRIPTION :
  --   Writes the complex number c to file, checking first whether
  --   it is imaginary, real, or integer.
  --   The cnt on return equals the number of characters written.

  procedure Write_Coefficient ( file : in file_type; c : in Complex_Number;
                                cnt : out natural32 );

  -- DESCRIPTION :
  --   Writes the coefficient before a nonconstant term.
  --   In the special case where
  --     c is +1 : nothing is written;
  --     c is -1 : only the sign is written;
  --     c is +i : then i* is written;
  --     c is -i : then -i* is written;
  --   otherwise the Write_Number handles the writing.
  --   Except for c = +/-1, a '*' is written.

  procedure Write_Plus ( file : in file_type; c : in Complex_Number;
                         cnt : out natural32 );

  -- DESCRIPTION :
  --   Writes the plus sign if the complex number c is positive real,
  --   cnt equals 1 on return if + was written, otherwise cnt is 0.

end Standard_Write_Numbers;
