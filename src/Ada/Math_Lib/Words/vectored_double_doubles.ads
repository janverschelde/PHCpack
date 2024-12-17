with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;

package Vectored_Double_Doubles is

-- DESCRIPTION :
--   The vectored sum of an array of double doubles postpones
--   the normalization of the sum to the very end.

  procedure Split ( v : in Double_Double_Vectors.Vector;
                    v0,v1,v2,v3 : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the doubles in double double numbers of v.

  -- ON ENTRY :
  --   v        a vector of double double numbers.

  -- ON RETURN :
  --   v0       high word of the high double of the numbers in v;
  --   v1       low word of the high double of the numbers in v;
  --   v2       high word of the low double of the numbers in v;
  --   v3       low word of the low double of the numbers in v.

  procedure Split ( v : in DoblDobl_Complex_Vectors.Vector;
                    v0re,v1re : out Standard_Floating_Vectors.Vector;
                    v2re,v3re : out Standard_Floating_Vectors.Vector;
                    v0im,v1im : out Standard_Floating_Vectors.Vector;
                    v2im,v3im : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the doubles of the real and imaginary parts of the
  --   complex double doubles in the vector v.

  -- ON ENTRY :
  --   v        a vector of complex double double numbers.

  -- ON RETURN :
  --   v0re     high word of the real high double of the numbers in v;
  --   v1re     low word of the real high double of the numbers in v;
  --   v2re     high word of the real low double of the numbers in v;
  --   v3re     low word of the real low double of the numbers in v;
  --   v0im     high word of the imaginary high double of the numbers in v;
  --   v1im     low word of the imaginary high double of the numbers in v;
  --   v2im     high word of the imaginary low double of the numbers in v;
  --   v3im     low word of the imaginary low double of the numbers in v.

  procedure Sum ( v0,v1,v2,v3 : in Standard_Floating_Vectors.Vector;
                  s0,s1,s2,s3 : out double_float );

  -- DESCRIPTION :
  --   Sums the numbers in v0, v1, v2, v3 into s0, s1, s2, s3.

  -- REQUIRED : v0'range = v1'range = v2'range = v3'range.

  procedure Sum ( v0re,v1re,v2re,v3re : in Standard_Floating_Vectors.Vector;
                  v0im,v1im,v2im,v3im : in Standard_Floating_Vectors.Vector;
                  s0re,s1re,s2re,s3re : out double_float;
                  s0im,s1im,s2im,s3im : out double_float );

  -- DESCRIPTION :
  --   Sums the numbers in v0re, v1re, v2re, v3re into s0re, s1re, s2re, s3re,
  --   and the numbers in v0im, v1im, v2im, v3im into s0im, s1im, s2im, s3im.

  -- REQUIRED :
  --   v0re'range = v1re'range = v2re'range = v3re'range,
  --   v0im'range = v1im'range = v2im'range = v3im'range.

  function to_Double_Double
             ( s0,s1,s2,s3 : double_float;
               verbose : boolean := true ) return double_double;
 
  -- DESCRIPTION :
  --   Given in s0, s1, s2, s3 the words of a sum,
  --   returns the double double number of the sum.
  --   If verbose, then the errors of the quick sum are shown.

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im : double_float;
               verbose : boolean := true ) return Complex_Number;
 
  -- DESCRIPTION :
  --   Given in s0re, s1re, s2re, s3re the real words of a sum,
  --   and in s0im, s1im, s2im, s3im the real words of a sum,
  --   returns the complex double double number of the sum.
  --   If verbose, then the errors of the quick sum are shown.

end Vectored_Double_Doubles;
