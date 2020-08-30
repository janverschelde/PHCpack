with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with DoblDobl_Complex_Vectors;

package Test_DoblDobl_Performance is

-- DESCRIPTION :
--   Tests the performance on double double complex vector computations.

  function Inner_Product
             ( x,y : DoblDobl_Complex_Vectors.Vector )
             return Complex_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors x and y,
  --   without taking complex conjugates.

  -- REQUIRED : x'range = y'range.

  procedure Inner_Product_Inlined
                ( zrehi,zimhi,zrelo,zimlo : out double_float;
                  xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                  xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                  yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  yimlo : in Standard_Floating_Vectors.Link_to_Vector);

  -- DESCRIPTION :
  --   Computes the inner product of two double double complex vectors
  --   x and y to form the result z, in 4-vector representation.
  --   No external procedures are called.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrehi      high parts of the real parts for the numbers in x;
  --   ximhi      high parts of the imaginary parts or the numbers in x;
  --   xrelo      low parts of the real parts for the numbers in x;
  --   ximlo      low parts of the imaginary parts for the numbers in x;
  --   yrehi      high parts of the real parts for the numbers in y;
  --   yimhi      high parts of the imaginary parts or the numbers in y;
  --   yrelo      low parts of the real parts for the numbers in y;
  --   yimlo      low parts of the imaginary parts for the numbers in y.

  -- ON RETURN :
  --   zrehi      high parts of the real parts of the inner product;
  --   zimhi      high parts of the imaginary parts of the inner product.
  --   zrelo      low parts of the real parts of the inner product;
  --   zimlo      low parts of the imaginary parts of the inner product.

  procedure Multiply
              ( first,second : in DoblDobl_Complex_Vectors.Link_to_Vector;
                product : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the coefficients of first with second
  --   and stores the results in the product.

  -- REQUIRED :
  --   All vectors have the same range.

  procedure Test_Add ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and tests their sum.
  --   The dimension of the vectors is given in dim.

  procedure Test_Update ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and tests the update.
  --   The dimension of the vectors is given in dim.

  procedure Test_Inner_Product ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and tests their inner product.
  --   The dimension of the vectors is given in dim.

  procedure Test_Multiply ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and tests their convolution.
  --   The dimension of the vectors is given in dim.

  procedure Timing_Add ( dim,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates random vectors and times their sum,
  --   for a frequency equal to frq.

  procedure Timing_Update ( dim,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates random vectors and times the update,
  --   for a frequency equal to frq.

  procedure Timing_Multiply ( dim,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates random vectors and times their convolution,
  --   for a frequency equal to frq.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_DoblDobl_Performance;
