with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;

package Test_Newton_Puiseux is

-- DESCRIPTION :
--   Tests the Newton-Puiseux algorithm on a random Laurent homotopy.

  procedure Define_Product_Homotopy
              ( dim : in integer32;
                nbm,nbt : in Standard_Integer_Vectors.Vector;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                hct : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given the dimension, number of monomials and number of terms,
  --   defines a Laurent homotopy with a power series solution.

  -- ON ENTRY :
  --   dim      number of equations and variables in the homotopy;
  --   nbm      nbm(i) equals the number of monomials in the system;
  --   nbt      nbt(i) equals the number of terms in the i-th series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cff      coefficients of the power series solution;
  --   pwr      real powers of the series solution;
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial.

  procedure Define_Binomial_Homotopy
              ( dim : in integer32;
                nbm : in Standard_Integer_Vectors.Vector;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                hct : out Standard_Floating_VecVecs.VecVec;
                intpow : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given the dimension, number of monomials and number of terms,
  --   defines a Laurent homotopy starting at canonical binomials.

  -- ON ENTRY :
  --   dim      number of equations and variables in the homotopy;
  --   nbm      nbm(i) equals the number of monomials in the system
  --            used to augment the canonical binomial system;
  --   intpow   if true, then converts the real powers of t to integers;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial.

  procedure Test_Product_Homotopy ( dim : in integer32 );

  -- DESCRIPTION :
  --   Runs a test on a Laurent homotopy of dimension dim,
  --   on a random system with a random series solution,
  --   which occurs as a factor of the homotopy.

  procedure Test_Binomial_Homotopy ( dim : in integer32 );

  -- DESCRIPTION :
  --   Runs a test on a Laurent homotopy of dimension dim,
  --   on a canonical binomial system augmented with a random Laurent system.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the the dimension and then launches the test.

end Test_Newton_Puiseux;
