with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;

package Extended_Random_Numbers is

-- DESCRIPTION :
--   The functions in this package take a give number and extend it
--   with random digits until a number of the given size is reached.
--   The extension is such that the most significant decimal places
--   are preserved in the number on return.  If the size is smaller
--   than the size of n, then the number on return has the same leading
--   decimal places as n, but is truncated up to the given size.

  function Extended_Random
             ( n : Natural_Number; size : natural32 ) return Natural_Number;

  function Extended_Random
             ( i : Integer_Number; size : natural32 ) return Integer_Number;

  function Extended_Random
             ( f : Floating_Number; size : natural32 ) return Floating_Number;

  function Extended_Random
             ( c : Complex_Number; size : natural32 ) return Complex_Number;

end Extended_Random_Numbers;
