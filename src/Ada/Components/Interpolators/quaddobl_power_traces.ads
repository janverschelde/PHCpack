with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;

package QuadDobl_Power_Traces is

-- DESCRIPTION :
--   The traces t_k of a degree d polynomial are linked with the signed
--   sums s_k of the k-th powers of the roots by the Newton identities :
--
--            k-1
--            ---
--     s   +  \     s  t     + k t  = 0,  for k = 1,2,..,d.
--      k     /      i  k-i       k
--            ---
--            i=1
--
--   With these relations we convert power sums into traces, and vice versa.
--   The power sums are defined by
--     s_k = x_1^k + x_2^k + .. + x_d^k,  for k = 1,2,..d,
--   and the traces by
--     t_k = (-1)^k times the sum of all possible k-products of the roots.
--   This version is for double double complex numbers.

  function Power_Sums_to_Trace
             ( s,t : Vector; k : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value at the k-th trace, given the power sums in s(1..k)
  --   and the previous traces in t(1..k-1).

  function Traces_to_Power_Sum 
             ( t,s : Vector; k : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value at the k-th power sum, given the traces in t(1..k)
  --   and previous power sums in s(1..k-1).

  function Power_Sums_to_Traces ( s : Vector ) return Vector;

  -- DESCRIPTION :
  --   Given power sums in s, the vector on return has the same range as s
  --   and contains the traces.

  function Traces_to_Power_Sums ( t : Vector ) return Vector;

  -- DESCRIPTION :
  --   Given traces in t, the vector on return has the same range as t
  --   and contains the power sums.

end QuadDobl_Power_Traces;
