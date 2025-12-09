with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;

package Test_Newton_Puiseux is

-- DESCRIPTION :
--   Tests the Newton-Puiseux algorithm on a random Laurent homotopy.

  procedure Evaluate_and_Differentiate
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                zt0 : in Standard_Complex_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                cjm : out Standard_Complex_Matrices.Matrix;
                ejm : out Standard_Floating_Matrices.Matrix;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates the Laurent homotopy at the leading constants
  --   and computes the Jacobian matrix.
 
  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   zt0      leading, constant coefficients of the series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      leading coefficients of the value;
  --   ydg      leading exponents of the value;
  --   cjm      leading coefficients of the Jacobian matrix;
  --   ejm      leading exponents of the Jacobian matrix.

  function Positive_Minimum
             ( v : Standard_Floating_Vectors.Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the smallest positive number in v.

  procedure Run_Newton_Step
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs one step with Newton's method on a Laurent homotopy,
  --   starting at the leading coefficients of a power series.

  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of a power series solution;
  --   pwr      exponents of a power series solution;
  --   vrblvl   is the verbose level.

  procedure Scale_Homotopy_Powers
              ( hct : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Subtracts the minimum power from all other powers in hct.

  procedure Define_Homotopy
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

  procedure Test ( dim : in integer32 );

  -- DESCRIPTION :
  --   Runs a test on a Laurent homotopy of dimension dim,
  --   on a random system with a random series solution.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the the dimension and then launches the test.

end Test_Newton_Puiseux;
