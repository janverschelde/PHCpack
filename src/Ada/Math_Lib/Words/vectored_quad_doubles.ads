with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;

package Vectored_Quad_Doubles is

-- DESCRIPTION :
--   The vectored operation on arrays of quad doubles postpones
--   the normalization of the result to the very end.

-- BASIC PROCEDURES :

  function Sign_Balance
             ( x : Quad_Double_Vectors.Vector; verbose : boolean := true )
             return Quad_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector where the numbers are sign balanced, that is:
  --   all parts have the same sign.

  procedure Signed_Quarter
              ( x,y : in Quad_Double_Vectors.Vector;
                xs0,xs1,xs2,xs3 : out Standard_Floating_Vectors.Vector;
                xs4,xs5,xs6,xs7 : out Standard_Floating_Vectors.Vector;
                xs8,xs9,xsA,xsB : out Standard_Floating_Vectors.Vector;
                xsC,xsD,xsE,xsF : out Standard_Floating_Vectors.Vector;
                ys0,ys1,ys2,ys3 : out Standard_Floating_Vectors.Vector;
                ys4,ys5,ys6,ys7 : out Standard_Floating_Vectors.Vector;
                ys8,ys9,ysA,ysB : out Standard_Floating_Vectors.Vector;
                ysC,ysD,ysE,ysF : out Standard_Floating_Vectors.Vector;
                xd0,xd1,xd2,xd3 : out Standard_Floating_Vectors.Vector;
                xd4,xd5,xd6,xd7 : out Standard_Floating_Vectors.Vector;
                xd8,xd9,xdA,xdB : out Standard_Floating_Vectors.Vector;
                xdC,xdD,xdE,xdF : out Standard_Floating_Vectors.Vector;
                yd0,yd1,yd2,yd3 : out Standard_Floating_Vectors.Vector;
                yd4,yd5,yd6,yd7 : out Standard_Floating_Vectors.Vector;
                yd8,yd9,ydA,ydB : out Standard_Floating_Vectors.Vector;
                ydC,ydD,ydE,ydF : out Standard_Floating_Vectors.Vector;
                ns,nd : out integer32 );

  -- DESCRIPTION :
  --   Quarters the numbers in x and y, taking into account their sign,
  --   assuming all numbers in x and y are signed balanced.

  -- REQUIRED : 
  --   All output vectors have the same range as x and y.

  -- ON ENTRY :
  --   x        a vector of quad double numbers;
  --   y        a vector of quad double numbers.

  -- ON RETURN :
  --   xs*, ys* : both high and low have the same sign;
  --   xd*, yd* : both high and low have different sign;
  --   ns       number of pairs (x(i),y(i)) with the same sign;
  --   nd       number of pairs (x(i),y(i)) with different signs.

  procedure Balanced_Quarter_Product
              ( dim : in integer32;
                x0,x1,x2,x3 : in Standard_Floating_Vectors.Vector;
                x4,x5,x6,x7 : in Standard_Floating_Vectors.Vector;
                x8,x9,xA,xB : in Standard_Floating_Vectors.Vector;
                xC,xD,xE,xF : in Standard_Floating_Vectors.Vector;
                y0,y1,y2,y3 : in Standard_Floating_Vectors.Vector;
                y4,y5,y6,y7 : in Standard_Floating_Vectors.Vector;
                y8,y9,yA,yB : in Standard_Floating_Vectors.Vector;
                yC,yD,yE,yF : in Standard_Floating_Vectors.Vector;
                s0,s1,s2,s3,s4,s5,s6,s7 : out double_float;
                s8,s9,sA,sB,sC,sD,sE,sF : out double_float );

  -- DESCRIPTION :
  --   Given balanced quarter quad doubles in vectors of length dim,
  --   returns the subsums of their inner product.

  -- ON ENTRY :
  --   dim      all vectors have range 1..dim;
  --   x0, ..., xF are the balanced quarter quad doubles of x;
  --   y0, ..., yF are the balanced quarter quad doubles of y.

  -- ON RETURN :
  --   s0, ..., sF are the subsums of the inner product of x and y.

  procedure Write_Subsums
              ( s0,s1,s2,s3,s4,s5,s6,s7 : in double_float );
  procedure Write_Subsums
              ( s0,s1,s2,s3,s4,s5,s6,s7 : in double_float;
                s8,s9,sA,sB,sC,sD,sE,sF : in double_float );

  -- DESCRIPTION :
  --   Writes all subsums of the convolutions of a product of two vectors.

  function to_quad_double
              ( s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
                verbose : boolean := true ) return quad_double;
  function to_quad_double
              ( s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
                s8,s9,sA,sB,sC,sD,sE,sF : double_float;
                verbose : boolean := true ) return quad_double;

  -- DESCRPTION :
  --   Adds the given doubles to a quad double.
  --   If verbose, then the subsums are written.

-- SIGN AWARE WRAPPERS :

  function Product ( x,y : Quad_Double_Vectors.Vector;
                     verbose : boolean := true ) return quad_double;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors of x and y,
  --   via sign aware quartering.

  -- REQUIRED : x'range = y'range.

  function Product ( x,y : QuadDobl_Complex_Vectors.Vector;
                     verbose : boolean := true ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors of x and y,
  --   via sign aware quartering.

  -- REQUIRED : x'range = y'range.

end Vectored_Quad_Doubles;
