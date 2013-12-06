with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Double_Double_Basics is

-- DESCRIPTION :
--   This package provides basic functions for double double arithmetic,
--   based on the file inline.h of the QD-2.3.9 software library.
--   The code is based on QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  procedure quick_two_sum ( a,b : in double_float;
                            s,err : out double_float );

  -- DESCRIPTION :
  --   Assuming |a| >= |b|, returns a+b and in err the error.

  -- ON ENTRY :
  --   a,b      two doubles: |a| >= |b|.

  -- ON RETURN :
  --   s        returned sum of a and b.
  --   err      error value, b - (s - a).

  procedure quick_two_diff ( a,b : in double_float;
                             s,err : out double_float );

  -- DESCRIPTION :
  --   Assuming |a| >= |b|, returns a-b and in err the error.

  -- ON ENTRY :
  --   a,b      two doubles: |a| >= |b|.

  -- ON RETURN :
  --   s        returned a minus b.
  --   err      error value, (a - s) - b.

  procedure two_sum ( a,b : in double_float; s,err : out double_float );

  -- DESCRIPTION :
  --   Computes fl(a+b) and err(a+b).

  -- ON ENTRY :
  --   a,b      two doubles.

  -- ON RETURN :
  --   s        approximation for the sum of a and b is returned;
  --   err      error of a + b.

  procedure two_diff ( a,b : in double_float; s,err : out double_float );

  -- DESCRIPTION :
  --   Computes fl(a-b) and err(a-b).

  -- ON ENTRY :
  --   a,b      two doubles.

  -- ON RETURN :
  --   s        approximation for the difference of a with b is returned;
  --   err      error of a - b.

  procedure split ( a : in double_float; hi,lo : out double_float );

  -- DESCRIPTION :
  --   Computes high and low word of a.

  -- ON ENTRY :
  --   a        some double float.

  -- ON RETURN :
  --   hi       high word of a;
  --   lo       low word of a.

  procedure two_prod ( a,b : in double_float; p,err : out double_float );

  -- DESCRIPTION :
  --   Computes fl(a*b) and err(a*b).

  -- ON ENTRY :
  --   a,b      two doubles.

  -- ON RETURN :
  --   p        returned approximation for a*b;
  --   err      error on the approximated product.

  procedure two_sqr ( a : in double_float; q,err : out double_float );

  -- DESCRIPTION :
  --   Computes fl(a*a) and err(a*a) faster than two_prod.

  function nint ( d : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns nearest integer to d.

  function aint ( d : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the truncated integer.

end Double_Double_Basics;
