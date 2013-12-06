with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Floating64_Numbers;        use Multprec_Floating64_Numbers;

package Multprec_Floating64_Numbers_io is

-- DESCRIPTION :
--   Basic i/o routines for multi-precision floating numbers.

  procedure get ( f : in out Floating_Number );
  procedure get ( file : in file_type; f : in out Floating_Number );
  procedure get ( f : in out Floating_Number; c : in out character );
  procedure get ( file : in file_type;
                  f : in out Floating_Number; c : in out character );

  -- DESCRIPTION :
  --   Reading of a multiprecision floating-point number f.
  --   On input, the character c can be the leading character of f,
  --   otherwise, the default is a space.  On output, c is the last
  --   character read, which is the first nondigit immediately after f.

  procedure put ( f : in Floating_Number );
  procedure put ( file : in file_type; f : in Floating_Number );

  procedure put ( f : in Floating_Number; fore,aft,exp : in natural32 );
  procedure put ( file : in file_type;
                  f : in Floating_Number; fore,aft,exp : in natural32 );

  -- DESCRIPTION :
  --   Formatted output of a floating-point number.

  -- ON ENTRY :
  --   f          floating-point number;
  --   fore       number of places before the decimal point, including sign,
  --              additional spaces will be introducted;
  --   aft        number of places after the decimal point,
  --              if the fraction of f is longer, then it will be truncated,
  --              otherwise zeros are introduced;
  --   exp        number of places for the exponent, including sign,
  --              if the exponent of f is longer, then it will be printed
  --              in full, otherwise zeros are introducted.

  procedure put ( f : in Floating_Number; dp : in natural32 );
  procedure put ( file : in file_type;
                  f : in Floating_Number; dp : in natural32 );

  -- DESCRIPTION : put(f,dp) = put(f,dp,dp,dp).

  function Character_Size ( f : Floating_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of characters the string representation
  --   of the floating-point number f occupies.

  procedure put ( s : out string; f : in Floating_Number );

  -- DESCRIPTION :
  --   Writes the number f to a string.
  -- REQUIRED :
  --   s'last = Character_Size(f).


end Multprec_Floating64_Numbers_io;
