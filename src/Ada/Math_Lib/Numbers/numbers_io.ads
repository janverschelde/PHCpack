with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;

package Numbers_io is

-- DESCRIPTION :
--   This package provides user friendly input routines.
--   The exception handlers allow for repeated trials.

  procedure Read_Positive ( p : out positive );
  procedure Read_Natural ( n : out natural32 );
  procedure Read_Integer ( i : out integer32 );
  procedure Read_Single_Float ( f : out single_float );
  procedure Read_Double_Float ( f : out double_float );
  procedure Read_Double_Double ( f : out double_double );
  procedure Read_Quad_Double ( f : out quad_double );

  -- DESCRIPTION :
  --   Reads a number from standard input.
  --   As long as the value obtained is not of the right type,
  --   the user will be asked to try again.

  function Number_of_Decimal_Places ( n : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of decimal places n occupies.
  --   This is useful for formatted output.

  procedure Read_Positive_Float ( f : in out double_float );
  procedure Read_Positive_Double_Double ( f : in out double_double );
  procedure Read_Positive_Quad_Double ( f : in out quad_double );

  -- DESCRIPTION :
  --   Reads a float and forces the user to enter again
  --   in case the number entered is negative.

end Numbers_io;
