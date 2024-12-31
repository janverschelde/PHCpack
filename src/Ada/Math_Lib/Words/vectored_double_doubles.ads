with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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

-- BASIC PROCEDURES :

  procedure Split ( v : in Double_Double_Vectors.Vector;
                    v0,v1,v2,v3 : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the doubles in double double numbers of v.

  -- REQUIRED : v0'range = v1'range = v2'range = v4'range = v'range.

  -- ON ENTRY :
  --   v        a vector of double double numbers.

  -- ON RETURN :
  --   v0       high word of the high double of the numbers in v;
  --   v1       low word of the high double of the numbers in v;
  --   v2       high word of the low double of the numbers in v;
  --   v3       low word of the low double of the numbers in v.

  procedure Signed_Split
              ( v : in Double_Double_Vectors.Vector;
                p0,p1,p2,p3 : out Standard_Floating_Vectors.Vector;
                m0,m1,m2,m3 : out Standard_Floating_Vectors.Vector;
                np0,np1,np2,np3,nm0,nm1,nm2,nm3 : out integer32 );

  -- DESCRIPTION :
  --   Splits the doubles in double double numbers of v,
  --   taking into account their signs.

  -- REQUIRED :
  --   m0'range = m1'range = m2'range = m4'range can fit all negative numbers,
  --   and p0'range = p1'range = p2'range = p4'range for all other numbers.

  -- ON ENTRY :
  --   v        a vector of double double numbers.

  -- ON RETURN :
  --   p0       high word of the high double of the nonnegative numbers in v;
  --   p1       low word of the high double of the nonnegative numbers in v;
  --   p2       high word of the low double of the nonnegative numbers in v;
  --   p3       low word of the low double of the nonnegative numbers in v;
  --   m0       high word of the high double of the negative numbers in v;
  --   m1       low word of the high double of the negative numbers in v;
  --   m2       high word of the low double of the negative numbers in v;
  --   m3       low word of the low double of the negative numbers in v;
  --   np0      number of elements in p0;
  --   np1      number of elements in p1;
  --   np2      number of elements in p2;
  --   np3      number of elements in p3;
  --   nm0      number of elements in m0;
  --   nm1      number of elements in m1;
  --   nm2      number of elements in m2;
  --   nm3      number of elements in m3.

  procedure Split ( v : in DoblDobl_Complex_Vectors.Vector;
                    v0re,v1re : out Standard_Floating_Vectors.Vector;
                    v2re,v3re : out Standard_Floating_Vectors.Vector;
                    v0im,v1im : out Standard_Floating_Vectors.Vector;
                    v2im,v3im : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the doubles of the real and imaginary parts of the
  --   complex double doubles in the vector v.

  -- REQUIRED : 
  --   v0re'range = v1re'range = v2re'range = v3re'range = v'range, and
  --   v0im'range = v1im'range = v2im'range = v3im'range = v'range.

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

  procedure Quarter ( v : in Double_Double_Vectors.Vector;
                      v0,v1,v2,v3 : out Standard_Floating_Vectors.Vector;
                      v4,v5,v6,v7 : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the high and low doubles of the numbers in v
  --   into four equal parts, with v0, v1, v2, v3 the quarters of the
  --   high doubles and v4, v5, v6, v7 the quarters of the low doubles.

  -- REQUIRED : 
  --   v0'range = v1'range = v2'range = v3'range = v'range, and
  --   v4'range = v5'range = v6'range = v7'range = v'range.

  -- ON ENTRY :
  --   v        a vector of double double numbers.

  -- ON RETURN :
  --   v0,v1,v2,v3 are the words of the highest doubles in v;
  --   v4,v5,v6,v7 are the words of the lowest doubles in v.

  function Signs ( f0,f1,f2,f3,f4,f5,f6,f7 : double_float ) return string;

  -- DESCRIPTION :
  --   Returns a string encoding the signs of the 8 given doubles.

  procedure Signed_Convolutions ( sx : in string );

  -- DESCRIPTION :
  --   Given the arrays of strings for x, writes the sign patterns
  --   when making the convolution products of x with itself.

  procedure Signed_Convolutions ( sx,sy : in string );

  -- DESCRIPTION :
  --   Given the arrays of strings for x and y, writes the sign patterns
  --   when making the convolution products of x with y.

  procedure Signed_Quarter
              ( v : in Double_Double_Vectors.Vector;
                x0,x1,x2,x3 : out Standard_Floating_Vectors.Vector;
                x4,x5,x6,x7 : out Standard_Floating_Vectors.Vector;
                pm0,pm1,pm2,pm3 : out Standard_Floating_Vectors.Vector;
                pm4,pm5,pm6,pm7 : out Standard_Floating_Vectors.Vector;
                mp0,mp1,mp2,mp3 : out Standard_Floating_Vectors.Vector;
                mp4,mp5,mp6,mp7 : out Standard_Floating_Vectors.Vector;
                nbx,npm,nmp : out integer32;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Quarters the numbers in v, taking into account their sign.

  -- REQUIRED : 
  --   All output vectors have the same range as v.

  -- ON ENTRY :
  --   v        a vector of double double numbers.

  -- ON RETURN :
  --   x0, .., x7 : all quarters have the same sign;
  --   pm0, .., pm7 : first four quarters are nonnegative,
  --                  last four quarters are negative;
  --   mp0, .., mp7 : first four quarters are negative,
  --                  last four quarters are nonnegative;
  --   nbx      number of quarters of the same sign;
  --   npm      number of quarters of the type (+, -);
  --   nmp      number of quarters of the type (-, +).

  procedure Signed_Quarter
              ( x,y : in Double_Double_Vectors.Vector;
                xs0,xs1,xs2,xs3 : out Standard_Floating_Vectors.Vector;
                xs4,xs5,xs6,xs7 : out Standard_Floating_Vectors.Vector;
                ys0,ys1,ys2,ys3 : out Standard_Floating_Vectors.Vector;
                ys4,ys5,ys6,ys7 : out Standard_Floating_Vectors.Vector;
                xd0,xd1,xd2,xd3 : out Standard_Floating_Vectors.Vector;
                xd4,xd5,xd6,xd7 : out Standard_Floating_Vectors.Vector;
                yd0,yd1,yd2,yd3 : out Standard_Floating_Vectors.Vector;
                yd4,yd5,yd6,yd7 : out Standard_Floating_Vectors.Vector;
                ns,nd : out integer32; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Quarters the numbers in x and y, taking into account their sign,
  --   assuming all numbers in x and y are signed balanced.

  -- REQUIRED : 
  --   All output vectors have the same range as x and y.

  -- ON ENTRY :
  --   x        a vector of double double numbers;
  --   y        a vector of double double numbers.

  -- ON RETURN :
  --   xs*, ys* : both high and low have the same sign;
  --   xd*, yd* : both high and low have different sign;
  --   ns       number of pairs (x(i),y(i)) with the same sign;
  --   nd       number of pairs (x(i),y(i)) with different signs.

  procedure Balanced_Quarter_Product
              ( dim : in integer32;
                x0,x1,x2,x3 : in Standard_Floating_Vectors.Vector;
                x4,x5,x6,x7 : in Standard_Floating_Vectors.Vector;
                y0,y1,y2,y3 : in Standard_Floating_Vectors.Vector;
                y4,y5,y6,y7 : in Standard_Floating_Vectors.Vector;
                s0,s1,s2,s3,s4,s5,s6,s7 : out double_float );

  -- DESCRIPTION :
  --   Given balanced quarter double doubles in vectors of length dim,
  --   returns the subsums of their inner product.

  -- ON ENTRY :
  --   dim      all vectors have range 1..dim;
  --   x0, ..., x7 are the balanced quarter quad doubles of x;
  --   y0, ..., y7 are the balanced quarter quad doubles of y.

  -- ON RETURN :
  --   s0, ..., s7 are the subsums of the inner product of x and y.

  procedure Quarter ( v : in DoblDobl_Complex_Vectors.Vector;
                      v0re,v1re : out Standard_Floating_Vectors.Vector;
                      v2re,v3re : out Standard_Floating_Vectors.Vector;
                      v4re,v5re : out Standard_Floating_Vectors.Vector;
                      v6re,v7re : out Standard_Floating_Vectors.Vector;
                      v0im,v1im : out Standard_Floating_Vectors.Vector;
                      v2im,v3im : out Standard_Floating_Vectors.Vector;
                      v4im,v5im : out Standard_Floating_Vectors.Vector;
                      v6im,v7im : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the high and the low doubles of the real and imaginary parts
  --   of the numbers in v into four equal parts.

  -- REQUIRED : 
  --   v0re'range = v1re'range = v2re'range = v3re'range = v'range,
  --   v4re'range = v5re'range = v6re'range = v7re'range = v'range,
  --   v0im'range = v1im'range = v2im'range = v3im'range = v'range, and
  --   v4im'range = v5im'range = v6im'range = v7im'range = v'range.

  -- ON ENTRY :
  --   v        a vector of complex double double numbers.

  -- ON RETURN :
  --   v0re,v1re,v2re,v3re are the words of the highest doubles
  --            of the real parts in v;
  --   v4re,v5re,v6re,v7re are the words of the lowest doubles
  --            of the real parts in v;
  --   v0im,v1im,v2im,v3im are the words of the highest doubles
  --            of the imaginary parts in v;
  --   v4im,v5im,v6im,v7im are the words of the lowest doubles
  --            of the imaginary parts in v.

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
  --   v0re'range = v1re'range = v2re'range = v3re'range, and
  --   v0im'range = v1im'range = v2im'range = v3im'range.

  procedure Product ( x0,x1,x2,x3 : in Standard_Floating_Vectors.Vector;
                      x4,x5,x6,x7 : in Standard_Floating_Vectors.Vector;
                      y0,y1,y2,y3 : in Standard_Floating_Vectors.Vector;
                      y4,y5,y6,y7 : in Standard_Floating_Vectors.Vector;
                      s0,s1,s2,s3,s4,s5,s6,s7 : out double_float );

  -- DESCRIPTION :
  --   Makes the inner product of the vectors x and y,
  --   given by quartered components in x0, x1, x2, x3 for x,
  --   and in y0, y1, y2, y3 for y.
  --   The resulting sums in s0, s1, s2, s3, s4, s5, s6, and s7
  --   correspond to the successive words in the result.

  -- REQUIRED :
  --   x0'range = x1'range = x2'range = x3'range, and
  --   y0'range = y1'range = y2'range = y3'range.

  -- ON ENTRY :
  --   x0,x1,x2,x3   quarters of the high doubles in x;
  --   x4,x5,x6,x7   quarters of the low doubles in x;
  --   y0,y1,y2,y3   quarters of the high doubles in y;
  --   y4,y5,y6,y7   quarters of the low doubles in y.

  -- ON RETURN :
  --   s0            highest word of the inner product;
  --   s1            second highest word of the inner product;
  --   s2            third highest word of the inner product;
  --   s3            fourth highest word of the inner product;
  --   s4            fourth lowest word of the inner product;
  --   s5            third lowest word of the inner product;
  --   s6            second lowest word of the inner product;
  --   s7            lowest word of the inner product.

  procedure Product ( x0re,x1re : in Standard_Floating_Vectors.Vector;
                      x2re,x3re : in Standard_Floating_Vectors.Vector;
                      x4re,x5re : in Standard_Floating_Vectors.Vector;
                      x6re,x7re : in Standard_Floating_Vectors.Vector;
                      x0im,x1im : in Standard_Floating_Vectors.Vector;
                      x2im,x3im : in Standard_Floating_Vectors.Vector;
                      x4im,x5im : in Standard_Floating_Vectors.Vector;
                      x6im,x7im : in Standard_Floating_Vectors.Vector;
                      y0re,y1re : in Standard_Floating_Vectors.Vector;
                      y2re,y3re : in Standard_Floating_Vectors.Vector;
                      y4re,y5re : in Standard_Floating_Vectors.Vector;
                      y6re,y7re : in Standard_Floating_Vectors.Vector;
                      y0im,y1im : in Standard_Floating_Vectors.Vector;
                      y2im,y3im : in Standard_Floating_Vectors.Vector;
                      y4im,y5im : in Standard_Floating_Vectors.Vector;
                      y6im,y7im : in Standard_Floating_Vectors.Vector;
                      s0re,s1re,s2re,s3re : out double_float;
                      s4re,s5re,s6re,s7re : out double_float;
                      s0im,s1im,s2im,s3im : out double_float;
                      s4im,s5im,s6im,s7im : out double_float );

  -- DESCRIPTION :
  --   Makes the inner product of complex vectors x and y,
  --   given by quartered components of the high and low doubles
  --   in their real and imaginary parts.
  --   The sums in s0re, s1re, .., s7re are the successive words
  --   in the real part of the result, and the successive words in
  --   the imaginary part of the result are in s0im, s1im, .., s7im.
  --   Note that the product is not complex conjugated.

  -- REQUIRED :
  --   x0'range = x1'range = x2'range = x3'range, and
  --   y0'range = y1'range = y2'range = y3'range.

  -- ON ENTRY :
  --   x0re,x1re,x2re,x3re are the quarters of the high doubles
  --             in the real parts of x;
  --   x4re,x5re,x6re,x7re are the quarters of the low doubles
  --             in the real parts of x;
  --   x0im,x1im,x2im,x3im are the quarters of the high doubles
  --             in the imaginary parts of x;
  --   x4im,x5im,x6im,x7im are the quarters of the low doubles
  --             in the imaginary parts of x;
  --   y0re,y1re,y2re,y3re are the quarters of the high doubles
  --             in the real parts of y;
  --   y4re,y5re,y6re,y7re are the quarters of the low doubles
  --             in the real parts of y;
  --   y0im,y1im,y2im,y3im are the quarters of the high doubles
  --             in the imaginary parts of y;
  --   y4im,y5im,y6im,y7im are the quarters of the low doubles
  --             in the imaginary parts of y.

  -- ON RETURN :
  --   s0re      highest word of the real part of the product;
  --   s1re      second highest word of the real part of the product;
  --   s2re      third highest word of the real part of the product;
  --   s3re      fourth highest word of the real part of the product;
  --   s4re      fourth lowest word of the real part of the product;
  --   s5re      third lowest word of the real part of the product;
  --   s6re      second lowest word of the real part of the product;
  --   s7re      lowest word of the real part of the product;
  --   s0im      highest word of the imaginary part of the product;
  --   s1im      second highest word of the imaginary part of the product;
  --   s2im      third highest word of the imaginary part of the product;
  --   s3im      fourth highest word of the imaginary part of the product;
  --   s4im      fourth lowest word of the imaginary part of the product;
  --   s5im      third lowest word of the imaginary part of the product;
  --   s6im      second lowest word of the imaginary part of the product;
  --   s7im      lowest word of the imaginary part of the product.

  procedure Write_Subsums ( s0,s1,s2,s3 : in double_float );

  -- DESCRIPTION :
  --   Writes all subsums given in the four arguments.

  procedure Write_Subsums ( s0,s1,s2,s3,s4,s5,s6,s7 : in double_float );

  -- DESCRIPTION :
  --   Writes all subsums given in the eight arguments.

  function to_Double_Double
             ( s0,s1,s2,s3 : double_float;
               verbose : boolean := true ) return double_double;
 
  -- DESCRIPTION :
  --   Given in s0, s1, s2, s3 the words of a sum,
  --   returns the double double number of the sum.
  --   If verbose, then the errors of the quick sum are shown.

  function to_Double_Double
             ( s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
               verbose : boolean := true ) return double_double;
 
  -- DESCRIPTION :
  --   Given in s0, s1, s2, s3, s4, s5, s6, s7 the words of a sum,
  --   returns the double double number of the sum.
  --   If verbose, then the errors of the quick sum are shown.

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im : double_float;
               verbose : boolean := true ) return Complex_Number;
 
  -- DESCRIPTION :
  --   Given in s0re, s1re, s2re, s3re the real words of a sum,
  --   and in s0im, s1im, s2im, s3im the imaginary words of a sum,
  --   returns the complex double double number of the sum.
  --   If verbose, then the errors of the quick sum are shown.

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re : double_float;
               s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im : double_float;
               verbose : boolean := true ) return Complex_Number;
 
  -- DESCRIPTION :
  --   Given in s0re, s1re, .., s7re the real words of a sum,
  --   and in s0im, s1im, .., s7im the imaginary words of a sum,
  --   returns the complex double double number of the sum.
  --   If verbose, then the errors of the quick sum are shown.

-- SIGN AWARE WRAPPERS :

  function Sum ( v : Double_Double_Vectors.Vector;
                 verbose : boolean := true ) return double_double;

  -- DESCRIPTION :
  --   Splits the numbers and returns the sum as a double double.
  --   For better accuracy, takes into account the signs of the numbers.

  function Sum ( v : DoblDobl_Complex_Vectors.Vector;
                 verbose : boolean := true ) return Complex_Number;

  -- DESCRIPTION :
  --   Splits the numbers and returns the sum as a double double.
  --   For better accuracy, takes into account the signs of the numbers.

  function Squared_Norm
             ( x : Double_Double_Vectors.Vector;
               verbose : boolean := true ) return double_double;
 
  -- DESCRIPTION :
  --   Returns the product of x with itself, with sign aware quartering.

  function Squared_Norm
             ( x : DoblDobl_Complex_Vectors.Vector;
               verbose : boolean := true ) return double_double;
 
  -- DESCRIPTION :
  --   Returns the product of x with itself, with sign aware quartering,
  --   wrapping the squared norm for double double vectors.

  function Product ( x,y : Double_Double_Vectors.Vector;
                     verbose : boolean := true ) return double_double;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors of x and y,
  --   via sign aware quartering.

  -- REQUIRED : x'range = y'range.

  function Product ( x,y : DoblDobl_Complex_Vectors.Vector;
                     verbose : boolean := true ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors of x and y,
  --   via sign aware quartering.

  -- REQUIRED : x'range = y'range.

end Vectored_Double_Doubles;
