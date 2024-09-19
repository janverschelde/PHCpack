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
                sxp : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates deg+1 random complex coefficients cff and
  --   corresponding exponents sxp for the terms in the series.
  --   The leading coefficient has exponent zero,
  --   all other exponents are in [1,2], sorted in increasing order.

  function Is_Sorted ( xp : Standard_Floating_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the sequence xp is sorted in increasing order,
  --   false otherwise.

  procedure Write_Exponential_Series
              ( file : in file_type;
                cff : in Standard_Complex_Vectors.Vector;
                sxp : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the exponential series given by coefficients cff
  --   and corresponding exponents sxp to file.

  procedure Test_Inverse
              ( deg : in integer32;
                cff : in Standard_Complex_Vectors.Vector;
                sxp : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes and tests the inverse for the series defined
  --   by the coefficients in cff and exponents in sxp.
  --   The deg on input is the number of independent exponents.

  procedure Test_Inverse ( deg : in integer32 );

  -- DESCRIPTION :
  --   Makes a random series truncated at degree deg,
  --   computes its inverse and the product to check.

  procedure Test_Sum ( adeg,bdeg : in integer32 );

  -- DESCRIPTION :
  --   Makes two random series truncated at degree adeg and bdeg,
  --   computes their sum and the difference to check.

  procedure Test_Multiplicative_Commutativity
              ( adeg,bdeg : in integer32;
                acf,bcf : in Standard_Complex_Vectors.Vector;
                axp,bxp : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests whether multiplying the first series with the second
  --   is the same as the second series with the first.

  procedure Test_Product
              ( adeg,bdeg : in integer32;
                acf,bcf : in Standard_Complex_Vectors.Vector;
                axp,bxp : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests on the series of degrees adeg and bdeg,
  --   with coefficients in acf and bcf,
  --   and corresponding exponents in axp and bxp,
  --   after extending the second series.

  procedure Test_Product ( adeg,bdeg : in integer32 );

  -- DESCRIPTION :
  --   Makes two random series truncated at degree adeg and bdeg,
  --   computes their product and the quotient to check.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for the number of terms,
  --   generates random exponential series,
  --   and tests the arithmetical operations.

end Test_Double_Exponentials;
