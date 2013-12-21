with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Quad_Double_Renormalizations is

-- DESCRIPTION :
--   This package provides basic functions for quad double arithmetic,
--   based on the file qd_inline.h of the QD-2.3.9 software library.
--   The code is based on QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  procedure quick_renorm ( c0,c1,c2,c3,c4 : in out double_float );

  -- DESCRIPTION :
  --   Does a quick renormalization of the number
  --   represented by the five given doubles.

  procedure renorm4 ( c0,c1,c2,c3 : in out double_float );

  -- DESCRIPTION :
  --   Does a renormalization of the number
  --   represented by the four given doubles.

  procedure renorm5 ( c0,c1,c2,c3,c4 : in out double_float );

  -- DESCRIPTION :
  --   Does a renormalization of the number
  --   represented by the five given doubles.

  procedure three_sum ( a,b,c : in out double_float );

  -- DESCRIPTION :
  --   Generalization of the two_sum of Double_Double_Basics.

  procedure three_sum2 ( a,b,c : in out double_float );

  -- DESCRIPTION :
  --   Less accurate version the three_sum from above.

  procedure quick_three_accum 
              ( a,b,s : in out double_float; c : in double_float );

  -- DESCRIPTION :
  --   Adds c to the double double pair (a,b).
  --   If the result does not fit in two doubles,
  --   then the sum is returned in s and (a,b) contains the remainder.
  --   Otherwise, the returned value in s is zero
  --   and on return (a,b) stores the sum.

end Quad_Double_Renormalizations;
