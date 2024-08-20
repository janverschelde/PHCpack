with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Double_Exponential_Arithmetic is

-- DESCRIPTION :
--   Provides basic arithmetic operations on complex exponential series,
--   in double precision.  An exponential series is represented by
--   (1) a complex vector of coefficients; and
--   (2) a corresponding vector of real exponents.

  function Inverse ( cff : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the inverse series.

  -- REQUIRED : cff(0) is nonzero.

  function Convolute ( a,b : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the product of the series with coefficients
  --   with the same exponents and same length.

  -- REQUIRED : a'last = b'last.

  procedure Add ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  tol : in double_float := 1.0e-14 );

  -- DESCRIPTION :
  --   Computes the sum of two exponential series,
  --   truncated at the same degree.
  --   If ccf'last equals twice the truncation degree,
  --   then the result contains all terms of the sum,
  --   and (a + b) - b - a = (a + b) - a - b = 0. 

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series;
  --   tol        tolerance used to omit zero coefficients.

  -- ON RETURN :
  --   ccf        coefficients of the sum;
  --   cxp        exponents of the sum.

  procedure Sub ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  tol : in double_float := 1.0e-14 );

  -- DESCRIPTION :
  --   Computes the difference of two exponential series,
  --   truncated at the same degree.
  --   If ccf'last equals twice the truncation degree,
  --   then the result contains all terms of the difference,
  --   and (a + b) - b - a = (a + b) - a - b = 0. 

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series;
  --   tol        tolerance used to omit zero coefficients.

  -- ON RETURN :
  --   ccf        coefficients of the difference;
  --   cxp        exponents of the difference.

  procedure Mul ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  prdcf,wrkcf : in out Standard_Complex_Vectors.Vector;
                  prdxp,wrkxp : in out Standard_Floating_Vectors.Vector;
                  tol : in double_float := 1.0e-14 );

  -- DESCRIPTION :
  --   Computes the product of two exponential series,
  --   truncated at the same degree.
  --   If ccf'last is (deg+1)*deg, where deg is the truncation degree,
  --   then all terms of the product can be stored.

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series;
  --   tol        tolerance to omit zero coefficients.

  -- ON RETURN :
  --   ccf        coefficients of the product;
  --   cxp        exponents of the product.
  --   prdcf      work space coefficients for product term;
  --   prdxp      work space exponents for product term;
  --   wrkcf      work space coefficients for increment;
  --   wrkxp      work space exponents for increment.

  procedure Div ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  invbcf,prdcf,wrkcf : in out Standard_Complex_Vectors.Vector;
                  prdxp,wrkxp : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the quotient of two exponential series,
  --   truncated at the same degree.

  -- REQUIRED : bcf(0) is nonzero.

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series.

  -- ON RETURN :
  --   ccf        coefficients of the quotient;
  --   cxp        exponents of the quotient;
  --   invbcf     coefficients of inverse of second series;
  --   prdcf      work space coefficients for product term;
  --   prdxp      work space exponents for product term;
  --   wrkcf      work space coefficients for increment;
  --   wrkxp      work space exponents for increment.

end Double_Exponential_Arithmetic;
