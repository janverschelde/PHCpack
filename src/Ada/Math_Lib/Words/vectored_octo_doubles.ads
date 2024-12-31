with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with OctoDobl_Complex_Numbers;           use OctoDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with Octo_Double_Vectors;
with OctoDobl_Complex_Vectors;

package Vectored_Octo_Doubles is

-- DESCRIPTION :
--   The vectored operation on arrays of octo doubles postpones
--   the normalization of the result to the very end.

-- BASIC PROCEDURES :

  function Sign_Balance
             ( x : Octo_Double_Vectors.Vector; verbose : boolean := true )
             return Octo_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector where the numbers are sign balanced, that is:
  --   all parts have the same sign.

  procedure Signed_Quarter
              ( x,y : in Octo_Double_Vectors.Vector;
                xs00,xs01,xs02,xs03 : out Standard_Floating_Vectors.Vector;
                xs04,xs05,xs06,xs07 : out Standard_Floating_Vectors.Vector;
                xs08,xs09,xs10,xs11 : out Standard_Floating_Vectors.Vector;
                xs12,xs13,xs14,xs15 : out Standard_Floating_Vectors.Vector;
                xs16,xs17,xs18,xs19 : out Standard_Floating_Vectors.Vector;
                xs20,xs21,xs22,xs23 : out Standard_Floating_Vectors.Vector;
                xs24,xs25,xs26,xs27 : out Standard_Floating_Vectors.Vector;
                xs28,xs29,xs30,xs31 : out Standard_Floating_Vectors.Vector;
                ys00,ys01,ys02,ys03 : out Standard_Floating_Vectors.Vector;
                ys04,ys05,ys06,ys07 : out Standard_Floating_Vectors.Vector;
                ys08,ys09,ys10,ys11 : out Standard_Floating_Vectors.Vector;
                ys12,ys13,ys14,ys15 : out Standard_Floating_Vectors.Vector;
                ys16,ys17,ys18,ys19 : out Standard_Floating_Vectors.Vector;
                ys20,ys21,ys22,ys23 : out Standard_Floating_Vectors.Vector;
                ys24,ys25,ys26,ys27 : out Standard_Floating_Vectors.Vector;
                ys28,ys29,ys30,ys31 : out Standard_Floating_Vectors.Vector;
                xd00,xd01,xd02,xd03 : out Standard_Floating_Vectors.Vector;
                xd04,xd05,xd06,xd07 : out Standard_Floating_Vectors.Vector;
                xd08,xd09,xd10,xd11 : out Standard_Floating_Vectors.Vector;
                xd12,xd13,xd14,xd15 : out Standard_Floating_Vectors.Vector;
                xd16,xd17,xd18,xd19 : out Standard_Floating_Vectors.Vector;
                xd20,xd21,xd22,xd23 : out Standard_Floating_Vectors.Vector;
                xd24,xd25,xd26,xd27 : out Standard_Floating_Vectors.Vector;
                xd28,xd29,xd30,xd31 : out Standard_Floating_Vectors.Vector;
                yd00,yd01,yd02,yd03 : out Standard_Floating_Vectors.Vector;
                yd04,yd05,yd06,yd07 : out Standard_Floating_Vectors.Vector;
                yd08,yd09,yd10,yd11 : out Standard_Floating_Vectors.Vector;
                yd12,yd13,yd14,yd15 : out Standard_Floating_Vectors.Vector;
                yd16,yd17,yd18,yd19 : out Standard_Floating_Vectors.Vector;
                yd20,yd21,yd22,yd23 : out Standard_Floating_Vectors.Vector;
                yd24,yd25,yd26,yd27 : out Standard_Floating_Vectors.Vector;
                yd28,yd29,yd30,yd31 : out Standard_Floating_Vectors.Vector;
                ns,nd : out integer32 );

  -- DESCRIPTION :
  --   Quarters the numbers in x and y, taking into account their sign,
  --   assuming all numbers in x and y are signed balanced.

  -- REQUIRED : 
  --   All output vectors have the same range as x and y.

  -- ON ENTRY :
  --   x        a vector of octo double numbers;
  --   y        a vector of octo double numbers.

  -- ON RETURN :
  --   xs*, ys* : both high and low have the same sign;
  --   xd*, yd* : both high and low have different sign;
  --   ns       number of pairs (x(i),y(i)) with the same sign;
  --   nd       number of pairs (x(i),y(i)) with different signs.

  procedure Balanced_Quarter_Product
              ( dim : in integer32;
                x00,x01,x02,x03 : in Standard_Floating_Vectors.Vector;
                x04,x05,x06,x07 : in Standard_Floating_Vectors.Vector;
                x08,x09,x10,x11 : in Standard_Floating_Vectors.Vector;
                x12,x13,x14,x15 : in Standard_Floating_Vectors.Vector;
                x16,x17,x18,x19 : in Standard_Floating_Vectors.Vector;
                x20,x21,x22,x23 : in Standard_Floating_Vectors.Vector;
                x24,x25,x26,x27 : in Standard_Floating_Vectors.Vector;
                x28,x29,x30,x31 : in Standard_Floating_Vectors.Vector;
                y00,y01,y02,y03 : in Standard_Floating_Vectors.Vector;
                y04,y05,y06,y07 : in Standard_Floating_Vectors.Vector;
                y08,y09,y10,y11 : in Standard_Floating_Vectors.Vector;
                y12,y13,y14,y15 : in Standard_Floating_Vectors.Vector;
                y16,y17,y18,y19 : in Standard_Floating_Vectors.Vector;
                y20,y21,y22,y23 : in Standard_Floating_Vectors.Vector;
                y24,y25,y26,y27 : in Standard_Floating_Vectors.Vector;
                y28,y29,y30,y31 : in Standard_Floating_Vectors.Vector;
                s00,s01,s02,s03,s04,s05,s06,s07 : out double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : out double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : out double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : out double_float );

  -- DESCRIPTION :
  --   Given balanced quarter octo doubles in vectors of length dim,
  --   returns the subsums of their inner product.

  -- ON ENTRY :
  --   dim      all vectors have range 1..dim;
  --   x00, ..., x31 are the balanced quarter octo doubles of x;
  --   y00, ..., y31 are the balanced quarter octo doubles of y.

  -- ON RETURN :
  --   s00, ..., s31 are the subsums of the inner product of x and y.


  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float );
  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float );
  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : in double_float );
  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : in double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : in double_float );

  -- DESCRIPTION :
  --   Writes all subsums of the convolutions of a product of two vectors.

  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                verbose : boolean := true ) return Octo_double;
  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                verbose : boolean := true ) return Octo_double;
  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
                verbose : boolean := true ) return Octo_double;
  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : double_float;
                verbose : boolean := true ) return Octo_double;

  -- DESCRPTION :
  --   Adds the given doubles to a octo double.
  --   If verbose, then the subsums are written.

-- SIGN AWARE WRAPPERS :

  function Product ( x,y : Octo_Double_Vectors.Vector;
                     verbose : boolean := true ) return octo_double;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors of x and y,
  --   via sign aware quartering.

  -- REQUIRED : x'range = y'range.

  function Product ( x,y : OctoDobl_Complex_Vectors.Vector;
                     verbose : boolean := true ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors of x and y,
  --   via sign aware quartering.

  -- REQUIRED : x'range = y'range.

end Vectored_Octo_Doubles;
