package body Double_Double_Basics is

  procedure quick_two_sum ( a,b : in double_float;
                            s,err : out double_float ) is
  begin
    s := a + b;
    err := b - (s - a);
  end quick_two_sum;

  procedure quick_two_diff ( a,b : in double_float;
                             s,err : out double_float ) is
  begin
    s := a - b;
    err := (a - s) - b;
  end quick_two_diff;

  procedure two_sum ( a,b : in double_float; s,err : out double_float ) is

    bb : double_float;

  begin
    s := a + b;
    bb := s - a;
    err := (a - (s - bb)) + (b - bb);
  end two_sum;

  procedure two_diff ( a,b : in double_float; s,err : out double_float ) is

    bb : double_float;

  begin
    s := a - b;
    bb := s - a;
    err := (a - (s - bb)) - (b + bb);
  end two_diff;

  procedure split ( a : in double_float; hi,lo : out double_float ) is

    QD_SPLITTER : constant double_float := 134217729.0; -- 2^27 + 1
    QD_SPLIT_THRESH : constant double_float := 6.69692879491417e+299; -- 2^996
    temp : double_float;
    aa : double_float := a;

  begin
    if ( a > QD_SPLIT_THRESH or a < -QD_SPLIT_THRESH ) then
      aa := aa*3.7252902984619140625E-09;  -- 2^-28
      temp := QD_SPLITTER * aa;
      hi := temp - (temp - aa);
      lo := a - hi;
      hi := hi*268435456.0;  -- 2^28
      lo := 268435456.0;     -- 2^28
    else
      temp := QD_SPLITTER * a;
      hi := temp - (temp - a);
      lo := a - hi;
    end if;
  end split;

  procedure two_prod ( a,b : in double_float; p,err : out double_float ) is

    a_hi,a_lo,b_hi,b_lo : double_float;

  begin
    p := a*b;
    split(a,a_hi,a_lo);
    split(b,b_hi,b_lo);
    err := ((a_hi*b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo;
  end two_prod;

  procedure two_sqr ( a : in double_float; q,err : out double_float ) is
 
    hi,lo : double_float;

  begin
    q := a*a;
    split(a,hi,lo);
    err := ((hi*hi - q) + 2.0*hi*lo) + lo*lo;
  end two_sqr;

  function nint ( d : double_float ) return double_float is

  -- Note: copied from "arm-05.txt" /usr/local/gnat/share/doc/gnat/txt:
  -- S'Floor denotes a function with the following specification:
  --   function S'Floor (X : T) return T
  -- The function yields the value Floor(X), i.e., the largest (most positive)
  -- integral value less than or equal to X.  When X is zero, the result has 
  -- the sign of X; a zero result otherwise has a positive sign.
  -- S'Ceiling denotes a function with the following specification:
  --   function S'Ceiling (X : T) return T
  -- The function yields the value Ceiling(X), i.e., the smallest
  -- (most negative) integral value greater than or equal to X.
  -- When X is zero, the result has the sign of X; a zero result otherwise
  -- has a negative sign when S'Signed_Zeros is True.

  begin
    if (d = double_float'floor(d)) 
     then return d;
     else return double_float'floor(d + 0.5);
    end if;
  end nint;

  function aint ( d : double_float ) return double_float is
  begin
    if d >= 0.0
     then return double_float'floor(d);
     else return double_float'ceiling(d);
    end if;
  end aint;

end Double_Double_Basics;
