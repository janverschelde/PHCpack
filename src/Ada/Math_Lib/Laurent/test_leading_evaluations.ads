with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Test_Leading_Evaluations is

-- DESCRIPTION :
--   Tests the evaluation of Laurent monomials at leading terms
--   of series with real powers.

  procedure Test_Monomial_Derivative
              ( deg : in Standard_Integer_Vectors.Vector;
                pwr : in Standard_Floating_Vectors.Vector;
                cff : in Standard_Complex_Vectors.Vector;
                idx : in integer32; err : out double_float );

  -- DESCRIPTION :
  --   Tests the derivative of the monomial with exponents in deg
  --   at the series with leading powers in pwr and coefficients in cff,
  --   with respect to the variable index idx.
  --   On return in err is the magnitude of the error of a random point test.

  procedure Test_Indexed_Derivative
              ( deg : in Standard_Integer_Vectors.Vector;
                cff : in Standard_Complex_Vectors.Vector;
                err : out double_float );

  -- DESCRPITION :
  --   Computes all derivatives up to the third order
  --   and compares the outcome of the indexed derivative function
  --   with the values of the dedicated derivative functions.
  --   Returns the sum of the errors.

  procedure Show_Indices ( dim,idxsum : in integer32 );

  -- DESCRIPTION :
  --   Writes all index vectors of range 1..dim and
  --   with sum of the indices equal to idxsum.

  procedure Test_Indexed_Monomial_Derivatives ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Laurent monomial in dim many variables
  --   and tests the indexed derivatives.

  procedure Test_Indexed_Derivatives;

  -- DESCRIPTION :
  --   Prompts for the number of variables and then tests the
  --   indexed derivatives of a random Laurent monomial.

  procedure Show_Numbers ( dim,nbr : in integer32 );

  -- DESCRIPTION :
  --   Writes all numbers of size dim with digits in 0..nbr.

  procedure Test_Number_Enumeration;

  -- DESCRIPTION :
  --   Prompts for the dimension and a base number
  --   and then all numbers in this base.

  procedure Test_Monomial ( dim : in integer32 );
               
  -- DESCRIPTION :
  --   Tests monomial evaluation and differentation in dim many variables
  --   at the leading terms of a series with real positive powers.

  procedure Test_Polynomial ( nbr,dim : in integer32 );

  -- DESCRIPTION :
  --   Tests polynomial evaluation and differentiation in dim many variables,
  --   of a polynomial with nbr many terms, at the leading terms of a series
  --   with real positive powers.

  procedure Test_System
              ( nbp,dim : in integer32;
                nbm : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests polynomial evaluation in dim many variables,
  --   of a system of nbp polynomials, where the i-th polynomial
  --   has nbm(i) many terms, at the leading terms of a series with 
  --   real positive powers.

  procedure Test_Homotopy
              ( dim : in integer32;
                nbm,nbt : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests generation of a Laurent homotopy of dim polynomials
  --   in dim variables, with number of monomials in nbm,
  --   with a generated power series solution, where the number
  --   of terms in each series is defined by nbt.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and then launches a test.

end Test_Leading_Evaluations;
