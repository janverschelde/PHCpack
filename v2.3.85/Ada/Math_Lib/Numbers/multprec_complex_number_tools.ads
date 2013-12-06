with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Multprec_Complex_Numbers;

package Multprec_Complex_Number_Tools is

-- DESCRIPTION :
--   Here are some tools for multi-precision numbers.

  function Round ( c : Multprec_Complex_Numbers.Complex_Number )
                 return Standard_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Rounds the multi-precision number to the nearest standard
  --   complex number.

  function Create ( f : double_float ) 
                  return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the complex number with real part equal to f.

  function Create ( re,im : double_float ) 
                  return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the complex number with real part equal to re
  --   and imaginary part equal to im.

  function Create ( c : Standard_Complex_Numbers.Complex_Number )
                  return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Converts the standard complex number into a multi-precision number.

  procedure Set_Size ( c : in out Multprec_Complex_Numbers.Complex_Number;
                       size : in natural32 );

  -- DESCRIPTION :
  --   Adjusts the size of the number to the given size.

end Multprec_Complex_Number_Tools;
