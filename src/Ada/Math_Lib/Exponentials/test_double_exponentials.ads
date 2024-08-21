with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Test_Double_Exponentials is

-- DESCRIPTION :
--   Tests arithmetic operations on complex exponential series.

  procedure Make_Random_Exponentials
              ( deg : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                exp : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates deg+1 random complex coefficients and
  --   corresponding exponents for the terms in the series.

  procedure Write_Exponential_Series
              ( file : in file_type;
                cff : in Standard_Complex_Vectors.Vector;
                exp : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the exponential series to file.

  procedure Test_Inverse ( deg : in integer32 );

  -- DESCRIPTION :
  --   Makes a random series truncated at degree deg,
  --   computes its inverse and the product to check.

  procedure Test_Sum ( adeg,bdeg : in integer32 );

  -- DESCRIPTION :
  --   Makes two random series truncated at degree adeg and bdeg,
  --   computes their sum and the difference to check.

  procedure Test_Product ( deg : in integer32 );

  -- DESCRIPTION :
  --   Makes two random series truncated at degree deg,
  --   computes their product and the quotient to check.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for the number of terms,
  --   generates random exponential series,
  --   and tests the arithmetical operations.

end Test_Double_Exponentials;
