with text_io;                            use text_io;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;

package Multprec_Write_Numbers is

-- DESCRIPTION :
--   This package provides utilities to write multiprecision complex numbers
--   in a more compact way when they are integer or real numbers.

  function Is_Imag ( c : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the real part of the number c equals 0.0.

  function Is_Real ( c : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the imaginary part of the number c equals 0.0.

  function Is_Positive_Real ( c : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Is_Real(c) and REAL_PART(c) > 0.

  function Is_Positive_Imag ( c : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Is_Imag(c) and IMAG_PART(c) > 0.

  procedure Write_Number ( file : in file_type; c : in Complex_Number );

  -- DESCRIPTION :
  --   Checks first whether the number is real or purely imaginary
  --   before writing it in the standard format.

  procedure Write_Coefficient ( file : in file_type; c : in Complex_Number );

  -- DESCRIPTION :
  --   Writes the coefficient before a nonconstant term.
  --   In the special case where
  --     c is +1 : nothing is written;
  --     c is -1 : only the sign is written;
  --     c is +i : then i* is written;
  --     c is -i : then -i* is written;
  --   otherwise the Write_Number handles the writing.
  --   Except for c = +/-1, a '*' is written.

  procedure Write_Plus ( file : in file_type; c : in Complex_Number );

  -- DESCRIPTION :
  --   Writes the plus sign if the complex number c is positive real.

end Multprec_Write_Numbers;
