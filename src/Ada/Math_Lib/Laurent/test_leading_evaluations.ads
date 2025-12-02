with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Test_Leading_Evaluations is

-- DESCRIPTION :
--   Tests the evaluation of Laurent monomials at leading terms
--   of series with real powers.

  function Random_Monomial
             ( dim,low,upp : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of exponents of dimension dim,
  --   with values randomly generated between low and upp.

  function Random_Polynomial
             ( nbr,dim,low,upp : integer32 )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns as many as nbr random exponents in a vector of range 1..nbr,
  --   of dimension dim, with values between low and upp.

  function Random_Leading_Powers
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..dim with random powers of series
  --   with real positive powers. 

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

  procedure Test_Monomial ( dim : in integer32 );
               
  -- DESCRIPTION :
  --   Tests monomial evaluation and differentation in dim many variables
  --   at the leading terms of a series with real positive powers.

  procedure Test_Polynomial ( nbr,dim : in integer32 );

  -- DESCRIPTION :
  --   Tests polynomial evaluation and differentiation in dim many variables,
  --   of a polynomial with nbr many terms, at the leading terms of a series
  --   with real positive powers.

  procedure Test_System ( nbp,nbr,dim : in integer32 );

  -- DESCRIPTION :
  --   Tests polynomial evaluation and differentiation in dim many variables,
  --   of a system of nbp polynomials, each with nbr many terms,
  --   at the leading terms of a series with real positive powers.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and then launches a test.

end Test_Leading_Evaluations;
