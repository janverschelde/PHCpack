with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Double_Integer_Arithmetic is

-- DESCRIPTION :
--   A double integer number is represented as an unevaluated sum
--   of two 64-bit integers, high and low, with value high*base + low.
--   The arithmetical operations return a carry over.
--   If the carry over is nonzero, then the double integer number
--   can not store the result of the operation correctly.
--   For the multiplications, the base should be 2^60,
--   that is: all 64-bit integers should have size no larger than 60.

  procedure Add ( xhi,xlo,yhi,ylo : in integer64; 
                  zhi,zlo,carry : out integer64;
                  verbose : in boolean := true );

  -- DESCRIPTION :
  --   Adds two double integer numbers.

  -- ON ENTRY :
  --   xhi       high word of the first integer x;
  --   xlo       low word of the first integer x;
  --   yhi       high word of the second integer y;
  --   ylo       low word of the second integer y;
  --   verbose   if verbose, then prints intermediate results.

  -- ON RETURN :
  --   zhi       high word of the sum x + y;
  --   zlo       low word of the sum x + y;
  --   carry     carry over, if nonzero, then a double integer number
  --             can no longer store the sum correctly.

  procedure Mul ( x,y : in integer64; 
                  zhi,zlo,carry : out integer64;
                  verbose : in boolean := true );

  -- DESCRIPTION :
  --   Multiplies two integer numbers into a double integer number.
  --   If both x and y do not exceed 60 bits, then there is no carry
  --   and the result are correct.

  -- ON ENTRY :
  --   x         first integer;
  --   y         second integer;
  --   verbose   if verbose, then prints intermediate results.

  -- ON RETURN :
  --   zhi       high word of the product x * y;
  --   zlo       low word of the product x * y;
  --   carry     carry over, if nonzero, then a double integer number
  --             can no longer store the product correctly.

  procedure Dbl_Mul ( xhi,xlo,yhi,ylo : in integer64; 
                      zhihi,zlohi,zhilo,zlolo,carry : out integer64;
                      verbose : in boolean := true );

  -- DESCRIPTION :
  --   Multiplies two double integer numbers into 
  --   a quadruple integer number.
  --   If the size of all xhi, xlo, yhi, ylo is not larger than 60 bits,
  --   then the result is correct.

  -- ON ENTRY :
  --   xhi       high word of the first integer;
  --   xlo       low word of the first integer;
  --   yhi       high word of the second integer;
  --   ylo       low word of the second integer;
  --   verbose   if verbose, then prints intermediate results.

  -- ON RETURN :
  --   zhihi     highest word of the product x * y;
  --   zlohi     second highest word of the product x * y;
  --   zhilo     second lowest word of the product x * y;
  --   zlolo     lowest word of the product x * y;
  --   carry     carry over, if nonzero, then a quadruple integer number
  --             can no longer store the product correctly.

end Double_Integer_Arithmetic;
