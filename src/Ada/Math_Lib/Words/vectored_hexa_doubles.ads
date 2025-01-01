with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Vectored_Hexa_Doubles is

-- DESCRIPTION :
--   The vectored operation on arrays of hexa doubles postpones
--   the normalization of the result to the very end.

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
                x32,x33,x34,x35 : in Standard_Floating_Vectors.Vector;
                x36,x37,x38,x39 : in Standard_Floating_Vectors.Vector;
                x40,x41,x42,x43 : in Standard_Floating_Vectors.Vector;
                x44,x45,x46,x47 : in Standard_Floating_Vectors.Vector;
                x48,x49,x50,x51 : in Standard_Floating_Vectors.Vector;
                x52,x53,x54,x55 : in Standard_Floating_Vectors.Vector;
                x56,x57,x58,x59 : in Standard_Floating_Vectors.Vector;
                x60,x61,x62,x63 : in Standard_Floating_Vectors.Vector;
                y00,y01,y02,y03 : in Standard_Floating_Vectors.Vector;
                y04,y05,y06,y07 : in Standard_Floating_Vectors.Vector;
                y08,y09,y10,y11 : in Standard_Floating_Vectors.Vector;
                y12,y13,y14,y15 : in Standard_Floating_Vectors.Vector;
                y16,y17,y18,y19 : in Standard_Floating_Vectors.Vector;
                y20,y21,y22,y23 : in Standard_Floating_Vectors.Vector;
                y24,y25,y26,y27 : in Standard_Floating_Vectors.Vector;
                y28,y29,y30,y31 : in Standard_Floating_Vectors.Vector;
                y32,y33,y34,y35 : in Standard_Floating_Vectors.Vector;
                y36,y37,y38,y39 : in Standard_Floating_Vectors.Vector;
                y40,y41,y42,y43 : in Standard_Floating_Vectors.Vector;
                y44,y45,y46,y47 : in Standard_Floating_Vectors.Vector;
                y48,y49,y50,y51 : in Standard_Floating_Vectors.Vector;
                y52,y53,y54,y55 : in Standard_Floating_Vectors.Vector;
                y56,y57,y58,y59 : in Standard_Floating_Vectors.Vector;
                y60,y61,y62,y63 : in Standard_Floating_Vectors.Vector;
                s00,s01,s02,s03,s04,s05,s06,s07 : out double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : out double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : out double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : out double_float;
                s32,s33,s34,s35,s36,s37,s38,s39 : out double_float;
                s40,s41,s42,s43,s44,s45,s46,s47 : out double_float;
                s48,s49,s50,s51,s52,s53,s54,s55 : out double_float;
                s56,s57,s58,s59,s60,s61,s62,s63 : out double_float );

  -- DESCRIPTION :
  --   Given balanced quarter hexa doubles in vectors of length dim,
  --   returns the subsums of their inner product.

  -- ON ENTRY :
  --   dim      all vectors have range 1..dim;
  --   x00, ..., x63 are the balanced quarter hexa doubles of x;
  --   y00, ..., y63 are the balanced quarter hexa doubles of y.

  -- ON RETURN :
  --   s00, ..., s63 are the subsums of the inner product of x and y.

  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : in double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : in double_float;
                s32,s33,s34,s35,s36,s37,s38,s39 : in double_float;
                s40,s41,s42,s43,s44,s45,s46,s47 : in double_float;
                s48,s49,s50,s51,s52,s53,s54,s55 : in double_float;
                s56,s57,s58,s59,s60,s61,s62,s63 : in double_float );

  -- DESCRIPTION :
  --   Writes all subsums of the convolutions of a product of two vectors.

  function to_hexa_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : double_float;
                s32,s33,s34,s35,s36,s37,s38,s39 : double_float;
                s40,s41,s42,s43,s44,s45,s46,s47 : double_float;
                s48,s49,s50,s51,s52,s53,s54,s55 : double_float;
                s56,s57,s58,s59,s60,s61,s62,s63 : double_float;
                verbose : boolean := true ) return hexa_double;

  -- DESCRPTION :
  --   Adds the given doubles to a hexa double.
  --   If verbose, then the subsums are written.

end Vectored_Hexa_Doubles;
