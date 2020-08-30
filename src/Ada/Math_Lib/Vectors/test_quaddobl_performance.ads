with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;

package Test_QuadDobl_Performance is

-- DESCRIPTION :
--   Development of better performing computations with vectors of
--   complex numbers in quad double precision.

  function Inner_Product
             ( x,y : QuadDobl_Complex_Vectors.Vector )
             return Complex_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors x and y,
  --   without taking complex conjugates.

  -- REQUIRED : x'range = y'range.

  procedure Multiply
              ( first,second : in QuadDobl_Complex_Vectors.Link_to_Vector;
                product : in QuadDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the coefficients of first with second
  --   and stores the results in the product.

  -- REQUIRED :
  --   All vectors have the same range.

  procedure Test_Add ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and tests their sum.
  --   The dimension of the vectors is given in dim.

  procedure Test_Two_Add ( dim : in integer32 );

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

  procedure Test_Multiply ( deg : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and tests their convolution.
  --   The degree of the coefficient vectors is given in deg.

  procedure Timing_Add ( dim,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and times the sum,
  --   for a frequency equal to the value of frq.
  --   The dimension of the vectors is given in dim.

  procedure Timing_Two_Add ( dim,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and times the sum,
  --   in the 2-vector representation, for a frequency frq.
  --   The dimension of the vectors is given in dim.

  procedure Timing_Update ( dim,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors of dimension dim 
  --   and times the update, for the frequency frq.

  procedure Timing_Inner_Product ( dim,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and times the inner product,
  --   in the 2-vector representation.
  --   The dimension of the vectors is in dim
  --   and the frequency of the inner products is in frq.

  procedure Timing_Multiply ( deg,frq : in integer32 );

  -- DESCRIPTION :
  --   Generates two random vectors and times the convolution
  --   in the 2-vector representation.
  --   The degree of the coefficient vectors is in deg
  --   and the frequency of the convolutions in frq.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_QuadDobl_Performance;
